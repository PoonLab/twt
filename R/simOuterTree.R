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
  stopifnot(all(names(counts)[2:ncol(counts)] == cnames))
  
  if (is.character(eventlog) & file.exists(eventlog)) {
    # load logfile into memory
    eventlog <- read.csv(eventlog)
  }
  counts <- get.counts(eventlog, mod)
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
  for (cn in names(sampling$targets)) {
    size <- sum(eventlog$event=="migration" & eventlog$dest==cn)
    target <- sampling$targets[[cn]]
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
        sum(active) == 1) {  # reached root of transmission tree
      break  
    }
    
    # retrieve the current event
    e <- eventlog[row, ]
    
    # lookahead to previous event
    e.prev <- ifelse(row>1, eventlog[row-1,], NULL)
    
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
      
      # what are the compartment sizes at transmission?
      n.source <- e.prev[[e$src]]
      stopifnot(n.source > 0)
      n.recip <- e.prev[[e$dest]]
      stopifnot(n.recip > 0)
      
      # is recipient an active Host? otherwise ignore
      n.active.recip <- active$count_type(e$dest)
      stopifnot(n.active.recip <= n.recip)
      
      if (sample(1:n.recip, 1) <= n.active.recip) {
        # TODO: check if Pathogen sampling date in recipient Host precedes
        # TODO: this transmission event
        
        recip <- active$sample.host(e$dest, remove=TRUE)
        
        # is source an active Host?
        n.active.source <- active$count.type(e$src)
        stopifnot(n.active.source <= n.source)
        
        if (sample(1:n.source, 1) <= n.active.source) {
          # source is an active host
          source <- active$sample.host(e$src)
        } else {
          # source is an unsampled host
          source <- Host$new(compartment=e$src)
        }
      }
    }
  }
}