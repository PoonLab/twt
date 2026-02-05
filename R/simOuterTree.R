#' sim.outer.tree
#'
#' Given the events generated under the model by forward-time simulation in 
#' simulate.dynamics(), we determine the subset of events that are involved 
#' in the transmission tree relating the sampled hosts.  This will generally
#' be a much smaller subset of events.
#' 
#' @param mod:  R6 object of class Model
#' @param eventlog:  character or data frame of event log from forward-time 
#'                   simulation of population trajectories
#' @return data frame, containing subset of events that comprise
#'         the transmission tree
#' @examples
#' require(twt)
#' mod <- Model$new(yaml.load_file("examples/SIRS_serial.yaml"))
#' eventlog <- sim.dynamics(mod)
#' outer <- sim.outer.tree(mod, eventlog)
#' @export
sim.outer.tree <- function(mod, eventlog, chunk.size=100) {
  cnames <- mod$get.compartments()
  k <- length(cnames)
  
  # append counts to event log
  if (is.character(eventlog)) {
    eventlog <- read.csv(eventlog)  # load logfile into memory
  }
  counts <- get.counts(eventlog, mod)
  stopifnot(all(names(counts)[2:ncol(counts)] == cnames))
  stopifnot(all(counts$time == c(0, eventlog$time)))
  for (cn in cnames) {
    eventlog[[cn]] <- counts[[cn]][2:nrow(counts)]
  }
  
  # confirm sufficient numbers in sampled compartments 
  sampling <- mod$get.sampling()
  if (sampling$mode != "compartment") {
    stop("Unsupported sampling mode", sampling$mode)
  }
  targets <- sampling$targets
  for (cn in names(targets)) {
    size <- sum(eventlog$event=="migration" & eventlog$dest==cn)
    target <- targets[[cn]]
    if (target < size) {
      stop("Insufficient number in sampling compartment ", 
           cn, " (", size , "<", target, ")")
    }
  }
  
  outer <- OuterTree$new(mod=mod)  # recording outer events
  active <- outer$get.active()
  
  # iterate through events in reverse (start with most recent)
  for (row in seq(nrow(eventlog), 1, -1)) {
    # check stopping criterion
    if (outer$nsamples() == sum(unlist(targets)) &  # sampled all hosts
        active$count.type() == 1) {  # reached root of transmission tree
      break  
    }
    
    e <- eventlog[row, ]  # retrieve the current event
    if (row > 1) {
      e.prev <- eventlog[row-1, ]  # previous event
    } else {
      e.prev <- e  # set counts to initial sizes
      e.prev$time <- 0
      e.prev$event <- NA
      e.prev$from.comp <- NA
      e.prev$to.comp <- NA
      e.prev$source <- NA
      for (cn in cnames) {
        e.prev[[cn]] <- mod$get.init.sizes()[[cn]]
      }
    }
    
    if (e$event == "migration") {
      # transition of a Host from one compartment to another
      .do.migration(e, outer)
    } else if (e$event == "transmission") {
      .do.transmission(e, e.prev, outer)
    }
    
    # birth events do not affect infected compartments, so they cannot 
    #  participate in the transmission tree
    
    # death events are ignored because we are tracing infected host lineages 
    #  back in time
  }
  return(outer)
}


#' Handle migration event
#' @param e:  row from event log including counts
#' @param active:  R6 object of class HostSet
#' @param outer:  R6 object of class OuterTree
#' @keywords internal
.do.migration <- function(e, outer) {
  targets <- unlist(outer$get.targets())
  from.comp <- e[['from.comp']]
  to.comp <- e[['to.comp']]
  active <- outer$get.active()
  
  if (to.comp %in% names(targets)) {
    # this was a sampling event
    
    if (outer$nsamples() < sum(targets)) {
      # have not sampled all hosts yet
      host <- Host$new(
        compartment=from.comp, # we are going back in time 
        sampling.time=e[['time']]
      )
      active$add.host(host)
      outer$add.sample(host)
    }
    
  } else {  
    # this was NOT a sampling event
    
    active.in.comp <- active$count.type(to.comp)
    if (active.in.comp > 0) {
      # check if migration affected active lineage
      
      tot <- e[[to.comp]]  # size of compartment AFTER migration
      prob <- active.in.comp / tot
      if (runif(1) < prob) {
        # choose an eligible Host at random to migrate
        host <- active$sample.host(to.comp)
        host$set.compartment(from.comp)
        
        event <- list(time=e[['time']], event=e[['event']], 
                   from.comp=from.comp, to.comp=to.comp,
                   from.host=NA, to.host=NA)
        outer$add.event(event)
      }
      # otherwise ignore this migration
    }
  }
  
  print(paste('migration', active$get.names()))
}



#' Handle transmission event
#' 
#' There is a lot of potential for confusion here.
#'   source: the Host that transmits its infection to the recipient Host
#'   recip:  the recipient Host to which the infection is transmitted
#'   from.comp:  the Compartment that the recipient Host was a member of 
#'               *before* becoming infected, e.g., 'S' susceptibles
#'   to.comp:  the Compartment that the recipient Host becomes a member of 
#'             *after* becoming infected, e.g., 'I' infected
#' 
#' @param e:  row from event log including counts
#' @param e.prev:  preceding row from event log
#' @param outer:  R6 object of class OuterTree
#' @keywords internal
.do.transmission <- function(e, e.prev, outer) {
  
  # how many recipients are in the destination compartment after transmission?
  n.recip <- e[[e$to.comp]]
  if (e$source == e$to.comp) {
    n.recip <- n.recip - 1  # recipient cannot be source
  }
  stopifnot(n.recip > 0)
  
  # how many of these hosts are active, i.e., carrying Pathogens?
  active <- outer$get.active()
  n.active.recip <- active$count.type(e$to.comp)
  if (n.active.recip > n.recip) {
    if (e$source == e$to.comp) {
     n.active.recip <- n.active.recip - 1 
    } else {
      stop("There are more active recipients (", n.active.recip, ") than ",
           "total recipients (", n.recip, ") at event ", e)      
    }
  }
  
  # is recipient an active Host? otherwise ignore
  if (sample(1:n.recip, 1) <= n.active.recip) {
    # TODO: check if Pathogen sampling date in recipient Host precedes
    #       this transmission event
    
    # is the recipient's previous Compartment an infected state?
    is.superinfect <- outer$get.infected(e$from.comp)
    recip <- active$sample.host(
      e$to.comp,  # pick an active member of the destination compartment
      remove=!is.superinfect  # if superinfection, host remains active
      )
    recip$set.transmission.time(e$time)  # can be appended
    
    # next identify source Host
    n.source <- e[[e$source]]
    stopifnot(n.source > 0)
    n.active.source <- active$count.type(e$source)
    stopifnot(n.active.source <= n.source)
    if (e$to.comp == e$source) { 
      n.source <- n.source - 1 
      }
    
    if (sample(1:n.source, 1) <= n.active.source) {
      # source is an active host
      source <- active$sample.host(e$source)
      
    } else {
      # source is an inactive (not directly sampled) host
      source <- Host$new(compartment=e$source, unsampled=TRUE)
      active$add.host(source)  # HostSet will assign a label
    }
    
    recip$set.source(source)  # record source (can be appended)
    outer$add.retired(recip)  # transfer Host from active to retired Set
    
    # record transmission event in outer log
    event <- list(
      time=e$time, event=e$event, 
      from.comp=e$from.comp, to.comp=e$to.comp,
      from.host=source$get.name(), to.host=recip$get.name()
      )
    print(paste('transmission', active$get.names()))
    outer$add.event(event)
  }
  
}
