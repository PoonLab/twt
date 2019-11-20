#' sim.inner.tree
#' 
#' Simulate the coalescence of Lineages within Compartments, and resolve
#' migration events that may involve sampled Lineages.
#' 
#' @param model:  R6 object of class Model or Run.  If user provides a Model 
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
sim.inner.tree <- function(model, e=NA) {
  
  # handle arguments
  if ( is.element('Run', class(mod)) ) {
    run <- mod
    if ( nrow(run$get.eventlog()$get.all.events()) ==0 ) {
      stop("Error in sim.inner.tree(): empty EventLogger in Run object. ",
           "Did you run sim.outer.tree() before calling this function?")
    }
  }
  else if ( is.element('Model', class(mod)) ) {
    if (is.environment(e)) {
      run <- Run$new(mod)
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
  if (any(!is.element(comps, union(events$compartment1, events$compartment2)))) {
    stop(paste("Mismatch between compartment names in Run and EventLogger objects."))
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
    c.times <- calc.coal.wait.times(run, current.time)
    
    if (length(c.times) > 0) {
      c.time <- min(c.times)
      t.delta <- e$time - current.time  # time left until next event
      
      if (c.time < t.delta) {
        comp.name <- names(c.times)[which.min(c.times)]
        comp <- comps[[comp.name]]  # retrieve Compartment obj from list
        
        current.time <- current.time + c.time
        resolve.coalescent(run, comp, time=current.time)
        
        next  # try another coalescent event (no change to row)
      }
    }
    
    # resolve the next event
    if (e$event.type == 'transmission') {
      resolve.transmission(run, e)
    }
    else if (e$event.type == 'migration') {
      resolve.migration(run, e)
    } 
    else if (e$event.type == 'transition') {
      resolve.transition(run, e)
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
  lineages <- run$get.extant.lineages(current.time, comp)
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
  if (!is.numeric(current.time)) {
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
    pair <- sample(linegaes, 2)  # replace defaults to FALSE
    resolve.coalecsent(run, comp, time, is.bottleneck=TRUE)
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
#' to another.  Note this reverts a Compartment at time 0 (most recent time)
#' from derived type2 to ancestral type1.
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
  type1 <- types[[e$type1]]
  type2 <- types[[e$type2]]
  
  # apply transition
  comp$set.type()
  
  # locate the original row
  eventlog <- run$get.eventlog()
  idx <- which(eventlog$event.type == 'transition' && 
                 eventlog$time == e$time && 
                 eventlog$compartment1 == e$compartment1)
}




#' rcoal.linear
#' 
#' Draw waiting time to next coalescent event for the current interval
#' of a piecewise linear model describing the population 
#' growth dynamics over time, where N(t) = alpha + beta*t
#' 
#' Inverse cumulative distribution function from Romero-Severson et al., 2017 
#' https://doi.org/10.1534/genetics.117.300284
#' 
#' @param k: number of extant lineages that could potentially coalesce
#' @param t1: most recent time point on piecewise interval -- origin of 
#'        coalescent (reverse-time) process.
#' @param alpha: intercept, population size at t=0 (start of infection)
#' @param beta: slope of linear growth interval
#' 
#' @return  Random variate from waiting time distribution.
#' 
#' @examples 
#' rcoal.linear()
#' 
#' @export
rcoal.linear <- function(k, t1, alpha, beta=0){
  u <- runif(1)  # random uniform deviate
  (1-(1-u)^(beta/choose(k,2)))*(alpha+beta*t1)/beta
}


#' calc.coal.wait.times
#' 
#' Draw waiting times for all compartments that have two or more extant lineages.
#' 
#' @param run: an R6 object of class Run'
#' @param current.time: current simulation time for inner tree
#' @return Named numeric vector of waiting times per compartment
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
#' calc.coal.wait.times(run, 0)
#' 
#' @export
calc.coal.wait.times <- function(run, current.time, dynamic=FALSE){
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
  waiting.times <- vector()
  
  for (comp in ext.comps) {
    k <- num.ext.lineages(comp)  # >=2 at this point
    this.type <- comp$get.type()
    
    if (is.null(this.type$get.death.rate.distr()) || 
        is.null(this.type$get.popn.growth.dynamics())) {
      
      # assume constant rate of coalescence
      c.rate <- this.type$get.coalescent.rate()
      waiting.times <- c(waiting.times, rexp(n=1, rate=choose(k,2)*c.rate))
      names(waiting.times)[length(waiting.times)] <- comp$get.name()
    }
    
    else {
      # population growth dynamics with user-specified piecewise linear model
      this.name <- comp$get.name()
      popn.growth <- this.type$get.popn.growth.dynamics()  # the model!
      
      # time is measured in reverse relative to start of simulation at t=0
      if (is.null(comp$get.branchingtime())) {
        # index case, branching time not determined by transmission from other case
        # arbitrarily set br. time to maximum time at limit of population growth model
        infection.time <- popn.growth[nrow(popn.growth), 'startTime']
      } else {
        infection.time <- comp$get.branching.time()
      }
      
      # amount of time until we reach start of infection (backwards)
      delta.t = overall.t <- infection.time - current.time
      
      # index to pieces within the remaining time interval
      piece.rows <- which(popn.growth[, 'startTime'] < delta.t)
      if (length(piece.rows) == 1) {
        # only one interval left - cast row vector as matrix
        pieces <- as.matrix(t(popn.growth[piece.rows, ]))
        row.names(pieces) <- 1
      } else {
        pieces <- popn.growth[piece.rows, ]
      }
      
      for (i in nrow(pieces):1) {
        piece <- pieces[i, ]
        wait <- wait.time(num.ext.lineages(this.name), piece['intercept'], piece['slope'])
        delta.t <- delta.t - wait
        
        if (delta.t < piece['startTime']) {
          # waiting time exceeds limit of this piece
          delta.t <- piece['startTime']
          if (delta.t == 0) {
            # past limit of initial piece - either a different within-host event occurs or 
            # we apply this compartment's bottleneck
            waiting.times <- c(waiting.times, overall.t)
            names(waiting.times)[[length(waiting.times)]] <- this.name
            break
          }
          # otherwise, go to next piece
        } else {
          # waiting time within this piece
          cumul.wait.time <- overall.t - delta.t  # relative to current time
          waiting.times <- c(waiting.times, cumul.wait.time)
          names(waiting.times)[[length(waiting.times)]] <- this.name
          break
        }
      }
      
    }  # end else
    
    # go to next compartment
  }
  
  return(waiting.times)
}

