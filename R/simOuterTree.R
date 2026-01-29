#' sim.outer.tree
#'
#' Given the events generated under the model by forward-time simulation in 
#' simulate.dynamics(), we determine the subset of events that are involved 
#' in the transmission tree relating the sampled hosts.  This will generally
#' be a much smaller subset of events
#' 
#' @param mod:  R6 object of class Model
#' @param eventlog:  character or data frame of event log from forward-time 
#'                   simulation of population trajectories
#' @return data frame, containing subset of events that comprise
#'         the transmission tree
#' @examples
#' sim <- simulate.dynamics(mod)
#' outer <- sim.outer.tree(mod, sim)
#' @export
sim.outer.tree <- function(mod, eventlog, chunk.size=100) {
  # container for recording outer events
  outer.log <- data.frame(
    time=numeric(),
    event=character(),
    host.anc=character(),
    host.des=character(),
    comp.anc=character(),
    comp.des=character()
  )
  
  cnames <- mod$get.compartments()
  k <- length(cnames)
  
  if (is.character(eventlog)) {
    eventlog <- read.csv(eventlog)  # load logfile into memory
  }
  counts <- get.counts(eventlog, mod)
  stopifnot(all(names(counts)[2:ncol(counts)] == cnames))
  stopifnot(all(counts$time == c(0, eventlog$time)))

  # append counts to event log
  for (cn in cnames) {
    eventlog[[cn]] <- counts[[cn]][2:nrow(counts)]
  }
  
  sampling <- mod$get.sampling()
  if (sampling$mode != "compartment") {
    stop("Unsupported sampling mode", sampling$mode)
  }
  
  # confirm sufficient numbers in sampled compartments
  targets <- sampling$targets
  for (cn in names(targets)) {
    size <- sum(eventlog$event=="migration" & eventlog$dest==cn)
    target <- targets[[cn]]
    if (target < size) {
      stop("Insufficient number in sampling compartment ", 
           cn, " (", size , "<", target, ")")
    }
  }
  nsamples <- sum(unlist(targets))
  sampled <- setNames(rep(0, length(targets)), names(targets))
  
  active <- HostSet$new()  # store active Host lineages
  
  # iterate through events in reverse (start with most recent)
  for (row in seq(nrow(eventlog), 1, -1)) {
    # check stopping criterion
    if (sum(sampled)==nsamples &  # sampled all hosts
        active$count.type() == 1) {  # reached root of transmission tree
      break  
    }
    
    e <- eventlog[row, ]  # retrieve the current event
    e.prev <- ifelse(row>1, eventlog[row-1,], NULL)  # previous event
    
    if (e$event == "migration") {
      # transition of a Host from one compartment to another
      
      if (e$dest %in% names(targets)) {
        # this was a sampling event
        
        if (sum(sampled) < nsamples) {
          host <- Host$new(
            name=paste0(e$dest),
            compartment=e$src, 
            sampling.time=e$time
            )
          active$add.host(host)
          sampled[[e$dest]] <- sampled[[e$dest]]+1
        }

      } else {  
        # this was NOT a sampling event
        
        # check if an active lineage underwent a migration
        active.in.comp <- active$count_type(e$dest)
        if (active.in.comp > 0) {
          tot <- e[[e$dest]]  # size of compartment AFTER migration
          prob <- active.in.comp / tot
          if (runif(1) < prob) {
            
            # choose an eligible Host at random to migrate
            host <- active$sample.host(e$dest)
            host$set.compartment(e$src)
            
            # record this event
            outer.log[nrow(outer.log)+1, ] <- 
              c(time=e$time, event=e$event, host.dec=host$get.name(), 
                comp.anc=e$src, comp.dec=e$dest)
          }
          # otherwise ignore this migration
        }
      }
      
    } else if (e$event == "transmission") {
      
      # what were the compartment sizes at transmission?
      n.recip <- ifelse(is.null(e.prev), e[[e$dest]]-1, e.prev[[e$dest]]) 
      stopifnot(n.recip > 0)
      
      # number of Hosts carrying Pathogens in destination compartment
      n.active.recip <- active$count_type(e$dest)
      if (e$dest == e$src) {
        # recipient Host cannot be its own source
        n.active.recip <- n.active.recip - 1
      }
      stopifnot(n.active.recip <= n.recip)
      
      # is recipient an active Host? otherwise ignore
      if (sample(1:n.recip, 1) <= n.active.recip) {
        # TODO: check if Pathogen sampling date in recipient Host precedes
        # TODO: this transmission event
        
        # if the recipient's Compartment is infected, this is a superinfection
        is.superinfect <- mod$get.infected[[e$dest]]
        
        # if this is not a superinfection, then retire the recipient Host
        recip <- active$sample.host(e$dest, remove=!is.superinfect)
        recip$set.transmission.time(e$time)
        
        # next identify source Host
        n.source <- e.prev[[e$src]]
        stopifnot(n.source > 0)
        
        # is source an active Host?
        n.active.source <- active$count.type(e$src)
        stopifnot(n.active.source <= n.source)
        
        if (sample(1:n.source, 1) <= n.active.source) {
          # source is an active host
          source <- active$sample.host(e$src)
        } else {
          # source is an unsampled host
          source <- Host$new(compartment=e$src, unsampled=TRUE)
          active$add.host(source)  # HostSet will assign a label
        }
        recip$set.source(source)  # record source
        
        # record transmission event in outer log
        outer.log[nrow(outer.log)+1, ] <- 
          c(time=e$time, event=e$event, 
            host.anc=source$get.name(), host.dec=recip$get.name(), 
            comp.anc=e$src, comp.dec=e$dest)
      }
      
    }
    # birth events do not affect infected compartments, so they cannot 
    # participate in the transmission tree
    
    # because we are tracing infected host lineages back in time, we can 
    # ignore death events
  }
  return(outer.log)
}
