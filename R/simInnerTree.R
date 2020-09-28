#' sim.inner.tree
#' 
#' Simulate the coalescence of Lineages within Compartments, and resolve
#' migration events that may involve sampled Lineages.
#' 
#' @param run:  R6 object of class Run.  The Run object tracks the locations of 
#'        Lineage objects among Compartments that are modified by transmission 
#'        and migration; and tracks the presence/absence of Lineages that change 
#'        with sampling, coalescence and bottleneck events.
#'        
#' @return R6 object of class EventLogger
#' 
#' @examples 
#' # load model
#' path <- system.file('extdata', 'SI.yaml', package='twt')
#' 
#' # load file and parse to construct MODEL object
#' settings <- yaml.load_file(path)
#' mod <- Model$new(settings)
#' 
#' # simulate outer tree - returns a Run object carrying EventLogger
#' run <- sim.outer.tree(mod)
#' 
#' # simulate inner tree
#' tree <- sim.inner.tree(run)
#' plot(tree)  # converts to a Phylo object
#' 
#' @export
sim.inner.tree <- function(run) {
  eventlog <- run$get.eventlog()
  events <- eventlog$get.all.events()
  
  # retrieve all Compartments in outer events
  comps <- c(run$get.compartments(), run$get.unsampled.hosts())
  if (any(!is.element(names(comps), union(events$compartment1, events$compartment2)))) {
    stop("Mismatch between compartment names in Run and EventLogger objects.")
  }
  types <- run$get.types()
  
  # revert Compartments to their derived ("recipient") Types
  # Types determine coalescent rates during simulation
  transitions <- events[events$event.type=='transition', ]
  for (i in order(transitions$time, decreasing=TRUE)) {
    e <- transitions[i,]
    comp <- comps[[e$compartment1]]
    dev.type <- types[[e$type1]]
    comp$set.type(dev.type)
  }
  
  # initialize simulation at most recent sampling time
  current.time <- 0.
  
  #extant.lineages <- run$get.extant.lineages(current.time)
  n.extant <- run$get.num.extant(current.time)
  if (n.extant == 0) {
    stop ('There must be at least one lineage sampled at time t=0.')
  }
  
  
  # iterate through events in reverse time (most recent first)
  events <- events[order(events$time), ]
  row <- 1
  while (row <= nrow(events)) {
    
    if ( n.extant == 1 & 
         run$get.num.extant(max(events$time)) == 1 ) {
         #length(run$get.extant.lineages(max(events$time)))==1 ) {
      # coalesced to final ancestral Lineage
      break
    }
    #print(n.extant)
    
    e <- events[row, ]  # retrieve outer event
    
    # draw waiting times for coalescence of extant lineages (named vector)
    c.times <- sample.coalescents(run, current.time)
    
    if (length(c.times) > 0) {
      c.time <- min(c.times)
      t.delta <- e$time - current.time  # time left until next event
      
      if (c.time < t.delta) {
        comp.name <- names(c.times)[which.min(c.times)]
        comp <- comps[[comp.name]]  # retrieve Compartment obj from list
        
        current.time <- current.time + c.time
        .resolve.coalescent(run, comp, time=current.time)
        
        #print('coalescent')
        next  # try another coalescent event (no change to row)
      }
    }
    
    # resolve the next event
    if (e$event.type == 'transmission') {
      .resolve.transmission(run, e)
    }
    else if (e$event.type == 'migration') {
      .resolve.migration(run, e)
    } 
    else if (e$event.type == 'transition') {
      .resolve.transition(run, e)
    } 
    else {
      stop("Error in simInnerTree: unknown event type ", e$event.type)
    }
    #print(e$event.type)
    
    # move to next event
    row <- row+1
    current.time <- e$time
    #extant.lineages <- run$get.extant.lineages(current.time)
    n.extant <- run$get.num.extant(current.time)
  }
 
  # coalesce residual Lineages in root Compartment
  while (n.extant > 1) {
    c.times <- sample.coalescents(run, current.time)
    if (length(c.times) == 0) {
      break
    }
    c.time <- min(c.times)
    comp.name <- names(c.times)[which.min(c.times)]
    comp <- comps[[comp.name]]
    
    current.time <- current.time + c.time
    .resolve.coalescent(run, comp, current.time)
    
    n.extant <- run$get.num.extant(current.time) #length(run$get.extant.lineages(current.time))
  }
  
  return(eventlog)
}


#' .resolve.coalescent
#' 
#' Helper function records the coalescence of two extant Lineages into an
#' ancestral Lineage.
#' 
#' @param run:  an R6 object of class Run
#' @param comp:  an R6 object of class Compartment
#' @param time:  double; time of coalescent event
#' @param is.bottleneck:  if TRUE, mark event as 'bottleneck' in log
#' 
#' @keywords internal
.resolve.coalescent <- function(run, comp, time, is.bottleneck=FALSE) {
  # retrieve extant lineages
  lineages <- run$get.extant.lineages(time, comp)
  if (length(lineages) < 2) {
    stop("Error in simInnerTree: <2 extant lineages in Compartment ", 
         comp$get.name())
  }
  pair <- sample(lineages, 2)
  line1 <- pair[[1]]
  line2 <- pair[[2]]
  
  # create ancestral Lineage
  ancestor <- Lineage$new(
    name=run$get.node.ident(),  # label internal nodes as "Node1", etc.
    sampling.time=time,
    location=comp
    )
  run$add.lineage(ancestor)
  
  # remove coalesced lineages
  run$remove.lineage(line1)
  run$remove.lineage(line2)
  
  # update event log
  event.type <- ifelse(is.bottleneck, 'bottleneck', 'coalescent')
  eventlog <- run$get.eventlog()
  eventlog$add.event(type=event.type, time=time, line1=line1$get.name(),
                     line2=ancestor$get.name(), comp1=comp$get.name())
  eventlog$add.event(type=event.type, time=time, line1=line2$get.name(),
                     line2=ancestor$get.name(), comp1=comp$get.name())
}


#' .resolve.transmission
#' 
#' Helper function records the transmission of Lineages from source
#' to recipient Compartments and updates the state of the Run object.
#' 
#' @param run:  R6 object of class Run
#' @param e:  row from EventLogger dataframe
#' 
#' @keywords internal
.resolve.transmission <- function(run, e) {
  if (e$event.type != 'transmission') {
    stop("resolve.transmission called on event of type '", e$event.type,
         "', expecting 'transmission'")
  }

  # retrieve Compartment objects for transmission event
  comps <- c(run$get.compartments(), run$get.unsampled.hosts())
  
  recipient <- comps[[e$compartment1]]
  if (is.null(recipient)) {
    stop("resolve.transmission() failed to locate recipient Compartment ",
         e$compartment1)
  }
  
  source <- comps[[e$compartment2]]
  if (is.null(source)) {
    stop("resolve.transmission() failed to locate source Compartment ", 
         e$compartment2)
  }
  
  # apply bottleneck in recipient
  survivors <- .resolve.bottleneck(run, recipient)
  for (lineage in survivors) {
    # move "surviving" Lineage to source Compartment
    recipient$remove.lineage(lineage)
    lineage$set.location(source)
    run$move.lineage(lineage, source)
    source$add.lineage(lineage)
  }
  
  # update eventlog
  eventlog <- run$get.eventlog()
  eventlog$record.transmission(recipient, survivors)
}


#' .resolve.bottleneck
#' 
#' Helper function records the bottleneck of extant Lineages, with a
#' random sample being transferred to the source Compartment.
#' 
#' @param run:  R6 object of class Run
#' @param comp:  R6 object of class Compartment
#' 
#' @return list of Lineage objects to pass to source Compartment
#' @keywords internal
.resolve.bottleneck <- function(run, comp) {
  
  time <- comp$get.branching.time()
  if (!is.numeric(time)) {
    stop("Error in .resolve.bottleneck: Compartment ", comp$get.name(), 
         " has no assigned branching time.")
  }
  
  mean.size <- comp$get.type()$get.bottleneck.size()
  theta <- comp$get.type()$get.bottleneck.theta()
  if (theta == 0) {
    bottleneck.size <- mean.size
  } else {
    bottleneck.size <- rztnbinom(n=1, mu=mean.size, theta=theta)  
  }
  
  if (!is.numeric(bottleneck.size)) {
    stop("Error in .resolve.bottleneck: Compartment ", comp$get.name(),
         " has no bottleneck size specified.")
  }
  
  # retrieve extant Lineages in compartment
  lineages <- run$get.extant.lineages(time, comp)
  if (length(lineages) == 0) {
    return(c())
  }
  
  while (length(lineages) > bottleneck.size) {
    pair <- sample(lineages, 2)  # replace defaults to FALSE
    .resolve.coalescent(run, comp, time, is.bottleneck=TRUE)
    # update lineages
    lineages <- run$get.extant.lineages(time, comp)
  }
  
  return(lineages)
}



#' .resolve.migration
#' 
#' Helper function records the migration of extant Lineages between
#' Compartments.
#' 
#' @param run:  R6 object of class Run
#' @param e:  a row from EventLogger$get.all.events() dataframe
#' 
#' @keywords internal
.resolve.migration <- function(run, e) {
  if (e$event.type != 'migration') {
    stop("resolve.migration called on event of type '", e$event.type,
         "', expecting 'migration'")
  }
  
  # retrieve Compartment objects for migration event
  comps <- c(run$get.compartments(), run$get.unsampled.hosts())
  
  recipient <- comps[[e$compartment1]]
  if (is.null(recipient)) {
    stop("resolve.migration() failed to locate recipient Compartment ",
         e$compartment1)
  }
  
  source <- comps[[e$compartment2]]
  if (is.null(source)) {
    stop("resolve.migration() failed to locate source Compartment ", 
         e$compartment2)
  }
  
  # retrieve extant lineages in recipient Compartment
  lineages <- run$get.extant.lineages(e$time, recipient)
  n.extant <- length(lineages)
  
  # determine how many extant Lineages in recipient (if any) were 
  # transferred by migration (secondary contact) from source
  eff.size <- recipient$get.type()$get.effective.size()
  bottleneck.size <- recipient$get.type()$get.bottleneck.size()
  
  # sampling without replacement (hypergeometric distribution)
  count <- rhyper(nn=1, m=n.extant, n=eff.size-n.extant,
                  k=bottleneck.size)
  
  if (count > 0) {
    migrants <- sample(lineages, count)
    
    for (line in migrants) {
      recipient$remove.lineage(line)
      line$set.location(source)
      source$add.lineage(line)
      
      # update Run object
      run$move.lineage(line, source)
    }
    
    # update eventlog
    eventlog <- run$get.eventlog()
    eventlog$record.migration(recipient, source, e$time, migrants)
  }

}


#' .resolve.transition
#' 
#' Helper function records the transition of a Compartment from one Type
#' to another.  Note this REVERTS a Compartment at time 0 (most recent time)
#' from derived type1 to ancestral type2.
#' 
#' @param run:  R6 object of class Run
#' @param e:  a row from an EventLogger dataframe
#' 
#' @keywords internal
.resolve.transition <- function(run, e) {
  # retrieve Compartment object by name
  comps <- c(run$get.compartments(), run$get.unsampled.hosts())
  comp <- comps[[e$compartment1]]
  
  # retrieve CompartmentTypes
  types <- run$get.types()
  type1 <- types[[e$type1]]  # derived
  type2 <- types[[e$type2]]  # ancestral
  
  # apply reverse transition (derived -> ancestral)
  comp$set.type(type2)
  
  # locate the original row
  eventlog <- run$get.eventlog()
  idx <- which(eventlog$event.type == 'transition' && 
                 eventlog$time == e$time && 
                 eventlog$compartment1 == e$compartment1)
}


#' sample.coalescents
#' 
#' Draw waiting times for all compartments that have two or more extant lineages.
#' 
#' @param run: an R6 object of class Run'
#' @param current.time: current simulation time for inner tree
#' @return Named numeric vector of waiting times per extant compartment
#' 
#' @examples 
#' # this model specifies 10 compartments with 3 lineages each,
#' # and constant coalescent rate 1.0
#' path <- system.file('extdata', 'SI.yaml', package='twt')
#' settings <- yaml.load_file(path)
#' mod <- Model$new(settings)
#' run <- Run$new(mod)
#' 
#' # these should average 0.333
#' sample.coalescents(run, 0)
#' 
#' @export
sample.coalescents <- function(run, current.time){
  # retrieve all compartments
  comps <- c(run$get.compartments(), run$get.unsampled.hosts())
  compnames <- names(comps)
  
  # retrieves compartment names for extant lineages
  ext.lineages.compnames <- run$get.locations()
  
  # retrieves compartments with multiple extant lineages
  counts <- table(unlist(ext.lineages.compnames))
  ext.comps <- comps[names(counts)[counts >= 2]]
  
  # calculate waiting times per Compartment
  sapply(ext.comps, function(comp) {
    k <- run$get.num.extant(current.time, comp$get.name())
    .rexp.coal(k, comp, current.time)
  })
}


#' rexp.coal
#' 
#' Draw exponential waiting time to next coalescent event of Lineages 
#' within a given Compartment.  Note that if the waiting time exceeds 
#' this Compartment's transmission event, it will be discarded by 
#' sim.inner.tree().
#' 
#' @param k:  integer, number of extant lineages
#' @param comp:  R6 object of class Compartment
#' @param time:  double, current simulation (reverse) time
#' 
#' @return  Random variate from waiting time distribution.
#' 
#' @keywords internal
.rexp.coal <- function(k, comp, time) {
  if (k < 2) {
    stop("Error in rexp.coal: called on Compartment with <2 lineages.")
  }
  
  b.time <- comp$get.branching.time()
  
  ctype <- comp$get.type()
  eff.size <- ctype$get.effective.size()
  gen.time <- ctype$get.generation.time()
  
  # FIXME: this is time-consuming
  pieces <- as.data.frame(ctype$get.popn.growth.dynamics())
  
  #k <- length(run$get.extant.lineages(time=time, comp=comp))
  #

  if ( is.null(b.time) || is.na(as.logical(b.time)) ) {
    # no branching time, handle index case
    if (is.null(pieces) || nrow(pieces) == 0) {
      if (is.null(eff.size)) {
        stop("Error in rexp.coal: CompartmentType ", ctype$get.name(), 
             " has no defined effective size nor growth model.")
      }
      return(gen.time*rexp(n=1, rate=choose(k,2) / eff.size))
    }
    else {
      # assume we are at last stage of piecewise growth
      # FIXME: should use infection time of index case
      eff.size <- pieces$startPopn[which.max(pieces$startTime)]
      return(gen.time*rexp(n=1, rate=choose(k,2) / eff.size))
    }
  }
  
  
  if (time > b.time) {
    stop("Error in rexp.coal: time ", time, " is further back than ",
         "infection time of Compartment ", comp$get.name())
  }
  
  
  if (is.null(pieces) || nrow(pieces) == 0) {
    if (is.null(eff.size)) {
      stop("Error in rexp.coal: CompartmentType ", ctype$get.name(), 
           " has no defined effective size nor growth model.")
    }
    else {
      # constant coalescent rate
      return(gen.time*rexp(n=1, rate=choose(k,2) / eff.size))
    }
  }
  else {
    # use piecewise linear growth model (overrides coalescent.rate)
    
    # get (forward) time since infection of Compartment
    f.time <- comp$get.branching.time() - time
    current.time <- f.time
    
    # iterate through pieces, starting with most recent (max startTime)
    pieces <- pieces[order(pieces$startTime, decreasing=TRUE), ]
    
    for (i in 1:nrow(pieces)) {
      piece <- pieces[i, ]
      if (piece$startTime > current.time) {
        # this interval is in the future!
        next
      }
      
      # configure parameters
      delta.t <- ifelse(is.na(piece$endTime), 
                        0, 
                        piece$endTime - current.time)
      beta <- piece$slope
      if (is.na(beta) & piece$startPopn == piece$endPopn) {
        beta <- 0
      }
      n.0 <- piece$endPopn
      n.eff <- n.0 + beta * delta.t
      
      # for derivation see https://github.com/PoonLab/twt/issues/90
      wait <- gen.time * 
        ifelse( beta == 0, 
               -n.eff/choose(k,2) * log(runif(1)),
               n.eff/beta * (runif(1)^(-beta/choose(k, 2)) - 1) )
      result <- time + (f.time - (current.time-wait))

      if (i == nrow(pieces) ||  # last interval
          (current.time - wait) > piece$startTime  # event within this interval
          ) {
        return(result)
      }
      
      # else go to next interval
      current.time <- piece$startTime  # should equal next piece's endTime
    }
  }
}

