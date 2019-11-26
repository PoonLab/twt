#' sim.inner.tree
#' 
#' Simulate the coalescence of Lineages within Compartments, and resolve
#' migration events that may involve sampled Lineages.
#' 
#' @param obj:  R6 object of class Model or Run.  If user provides a Model 
#'        object, then an EventLogger object (such as produced from a Newick tree
#'        string by `eventlog.from.tree`) must be provided.
#'        The supplied or derived Run object tracks the locations of Lineage objects
#'        among Compartments that are modified by transmission and migration; and 
#'        tracks the presence/absence of Lineages that change with sampling, 
#'        coalescence and bottleneck events.  
#' @param e: (optional)  R6 object of class EventLogger
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
sim.inner.tree <- function(obj, e=NA) {
  
  # handle arguments
  if ( is.element('Run', class(obj)) ) {
    run <- obj
    if ( nrow(run$get.eventlog()$get.all.events()) == 0 ) {
      stop("Error in sim.inner.tree(): empty EventLogger in Run object. ",
           "Did you run sim.outer.tree() before calling this function?")
    }
  }
  else if ( is.element('Model', class(obj)) ) {
    if (is.environment(e)) {
      run <- Run$new(obj)
      # don't modify the original EventLogger object
      run$set.eventlog(e$clone())
    } 
    else {
      stop("You must provide an EventLogger to sim.inner.tree() if `mod` is a Model object")
    }
  }
  else {
    stop("sim.inner.tree(): argument `mod` must be R6 class Model or Run.")
  }
  
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
  for (comp in comps) {
    temp <- transitions[transitions$comp1 == comp$get.name()]
    if (nrow(temp) == 0) {
      next  # this Compartment did not undergo any transitions
    }
    last.type <- temp$type1[which.min(temp$time)]
    ctype <- types[[last.type]]
    comp$set.type(ctype)
  }
  
  # initialize simulation at most recent sampling time
  current.time <- 0.
  extant.lineages <- run$get.extant.lineages(current.time)
  n.extant <- length(extant.lineages)
  if (n.extant == 0) {
    stop ('There must be at least one lineage sampled at time t=0.')
  }
  
  
  # iterate through events in reverse time (most recent first)
  events <- events[order(events$time), ]
  row <- 1
  while (row <= nrow(events)) {
    
    if (n.extant == 1) {
      # coalesced to final ancestral Lineage
      break
    }
    
    e <- events[row,]  # retrieve outer event
    
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
      
    # move to next event
    row <- row+1
    current.time <- e$time
    extant.lineages <- run$get.extant.lineages(current.time)
    n.extant <- length(extant.lineages)
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
    
    n.extant <- length(run$get.extant.lineages(current.time))
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
    stop("Error in generate.bottleneck: Compartment ", comp$get.name(), 
         " has no assigned branching time.")
  }
  
  bottleneck.size <- comp$get.type()$get.bottleneck.size()
  if (!is.numeric(bottleneck.size)) {
    stop("Error in generate.bottleneck: Compartment ", comp$get.name(),
         " has no bottleneck size specified.")
  }
  
  # retrieve extant Lineages in compartment
  lineages <- run$get.extant.lineages(time, comp)
  if (length(lineages) == 0) {
    stop("Error in generate.bottleneck: Compartment ", comp$get.name(),
         " has no Lineages to bottleneck!")
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
  eff.size <- 1/recipient$get.type()$get.coalescent.rate()
  bottleneck.size <- recipient$get.type()$get.bottleneck.size()
  
  # sampling without replacement (hypergeometric distribution)
  count <- rhyper(nn=1, m=n.extant, n=eff.size-n.extant,
                  k=bottleneck.size)
  migrants <- sample(lineages, count)
  
  for (line in migrants) {
    recipient$remove.lineage(line)
    line$set.location(source)
    source$add.lineage(line)
  }
  
  # update eventlog
  eventlog <- run$get.eventlog()
  eventlog$record.migration(recipient, source, e$time, migrants)
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
#' mod <- MODEL$new(settings)
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
  ext.lineages.compnames <- sapply(
    run$get.extant.lineages(current.time),
    function(x) {
      x$get.location()$get.name()
    }
  )
  
  # counts the number of extant lineages in one compartment
  # calculated for parameter `k` in function `wait.time`
  num.ext.lineages <- function(x) {
    length(which(ext.lineages.compnames==x$get.name()))
  }
  
  # retrieves compartments with multiple extant lineages
  counts <- sapply(comps, function(x) num.ext.lineages(x))
  ext.comps <- comps[counts >= 2]
  
  # calculate waiting times per Compartment
  sapply(ext.comps, function(comp) {
    .rexp.coal(run, comp, current.time)
  })
}


#' rexp.coal
#' 
#' Draw exponential waiting time to next coalescent event of Lineages 
#' within a given Compartment.  Note that if the waiting time exceeds 
#' this Compartment's transmission event, it will be discarded by 
#' sim.inner.tree().
#' 
#' @param run:  R6 object of class Run
#' @param comp:  R6 object of class Compartment
#' @param time:  double, current simulation (reverse) time
#' 
#' @return  Random variate from waiting time distribution.
#' 
#' @keywords internal
.rexp.coal <- function(run, comp, time) {
  b.time <- comp$get.branching.time()
  
  ctype <- comp$get.type()
  c.rate <- ctype$get.coalescent.rate()
  pieces <- as.data.frame(ctype$get.popn.growth.dynamics())
  
  k <- length(run$get.extant.lineages(time=time, comp=comp))
  if (k < 2) {
    stop("Error in rexp.coal: called on Compartment with <2 lineages.")
  }
  
  if (is.null(b.time) || is.na(b.time)) {
    # no branching time, handle index case
    if (is.null(pieces)) {
      if (is.null(c.rate)) {
        stop("Error in rexp.coal: CompartmentType ", ctype$get.name(), 
             " has no defined coalescent rate nor growth model.")
      }
      return(rexp(n=1, rate=choose(k,2) * c.rate))
    }
    else {
      # assume we are at last stage of piecewise growth
      # TODO: it would be more accurate to use sampling time distribution
      #       to determine infection time of index case (but what if index
      #       case is unsampled?)
      c.rate <- 1. / pieces$startPopn[which.max(pieces$startTime)]
      return(rexp(n=1, rate=choose(k,2) * c.rate))
    }
  }
  
  
  if (time > b.time) {
    stop("Error in rexp.coal: time ", time, " is further back than ",
         "infection time of Compartment ", comp$get.name())
  }
  
  
  if (is.null(pieces)) {
    if (is.null(c.rate)) {
      stop("Error in rexp.coal: CompartmentType ", ctype$get.name(), 
           " has no defined coalescent rate nor growth model.")
    }
    else {
      # constant coalescent rate
      return(rexp(n=1, rate=choose(k,2)*c.rate))
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
      delta.t <- piece$endTime - current.time
      beta <- piece$slope
      n.0 <- piece$endPopn
      n.eff <- n.0 + beta * delta.t
      
      # for derivation see https://github.com/PoonLab/twt/issues/90
      wait <- ifelse(beta==0, 
                     -n.eff/choose(k,2) * log(runif(1)),
                     n.eff/beta * (runif(1)^(-beta/choose(k, 2)) - 1)
                     )
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

