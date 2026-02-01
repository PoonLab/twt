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
  active <- HostSet$new()  # store active Host lineages
  
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
      e.prev <- NULL  # at start of event log
    }
    
    if (e$event == "migration") {
      # transition of a Host from one compartment to another
      .do.migration(e, active, outer.log)
    } else if (e$event == "transmission") {
      .do.transmission(e, e.prev, active, outer.log)
    }
    
    # birth events do not affect infected compartments, so they cannot 
    #  participate in the transmission tree
    
    # death events are ignored because we are tracing infected host lineages 
    #  back in time
  }
  return(outer.log)
}


#' Handle migration event
#' @param e:  row from event log
#' @param active:  R6 object of class HostSet
#' @param outer:  R6 object of class OuterTree
.do.migration <- function(e, active, outer) {
  targets <- outer$get.targets()
  
  if (e$dest %in% names(targets)) {
    # this was a sampling event
    
    if (outer$nsamples() < sum(targets)) {
      # have not sampled all hosts yet
      host <- Host$new(
        name=paste0(e$dest),
        compartment=e$src, 
        sampling.time=e$time
      )
      active$add.host(host)
      outer$add.sample(e$dest)
    }
    
  } else {  
    # this was NOT a sampling event
    
    active.in.comp <- active$count.type(e$dest)
    if (active.in.comp > 0) {
      # check if migration affected active lineage
      
      tot <- e[[e$dest]]  # size of compartment AFTER migration
      prob <- active.in.comp / tot
      if (runif(1) < prob) {
        # choose an eligible Host at random to migrate
        host <- active$sample.host(e$dest)
        host$set.compartment(e$src)
        
        event <- c(time=e$time, event=e$event, host.dec=host$get.name(), 
                   comp.anc=e$src, comp.dec=e$dest)
        outer$add.event(event)
      }
      # otherwise ignore this migration
    }
  }
}



#' Handle transmission event
#' 
#' There is a lot of potential for confusion here.
#'   source: the Host that transmits its infection to another
#'   recip:  the recipient Host to which the infection is transmitted
#'   src:  the Compartment that the recipient Host was a member of *before*
#'         becoming infected, e.g., 'S' susceptibles
#'   dest:  the Compartment that the recipient Host becomes a member of 
#'          *after* becoming infected, e.g., 'I' infected
#' 
#' @param e:  row from event log
#' @param e.prev:  preceding row from event log
#' @param active:  R6 object of class HostSet
#' @param outer:  R6 object of class OuterTree
.do.transmission <- function(e, e.prev, active, outer) {
  
  # how many recipients are in the destination compartment after transmission?
  n.recip <- e[[e$dest]]
  stopifnot(n.recip > 0)
  
  # number of Hosts carrying sampled Pathogens in destination compartment
  n.active.recip <- active$count.type(e$dest)
  if (n.active.recip > n.recip) {
    stop("There are more active recipients (", n.active.recip, ") than ",
         "total recipients (", n.recip, ") at event ", e)
  }
  
  # is recipient an active Host? otherwise ignore
  if (sample(1:n.recip, 1) <= n.active.recip) {
    # TODO: check if Pathogen sampling date in recipient Host precedes
    # TODO: this transmission event
    
    # if the recipient's Compartment is infected, this is a superinfection
    is.superinfect <- outer$get.infected()[[e$src]]
    
    # if this is not a superinfection, then retire the recipient Host
    recip <- active$sample.host(e$src, remove=!is.superinfect)
    recip$set.transmission.time(e$time)  # can be appended
    
    # next identify source Host - note the recipient enters the same 
    # compartment as the source: the *destination* compartment!
    n.source <- e.prev[[e$dest]]
    stopifnot(n.source > 0)
    
    # is source an active Host?
    n.active.source <- active$count.type(e$src)
    stopifnot(n.active.source <= n.source)
    
    if (sample(1:n.source, 1) <= n.active.source) {
      # source is an active host
      source <- active$sample.host(e$dest)
      
    } else {
      # source is an unsampled host
      source <- Host$new(compartment=e$src, unsampled=TRUE)
      active$add.host(source)  # HostSet will assign a label
    }
    
    recip$set.source(source)  # record source
    
    # record transmission event in outer log
    event <- c(time=e$time, event=e$event, 
               host.anc=source$get.name(), host.dec=recip$get.name(), 
               comp.anc=e$src, comp.dec=e$dest)
    outer$add.event(event)
  }
  
}
