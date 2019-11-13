# issue 48: Fix implementation of migration rates
# Reason: difficult to conceive and parameterize a migration process in reverse time
# Solution: Annotate outer transmission tree with migration events (going forward in time)
#           When we encounter a migration event during inner tree simulation in reverse time,
#           we check the probability that one of the sampled lineages we are following back
#           in time was a result of the migration.
#           This also means that a sampled compartment can become empty of sampled lineages
#           before we reach its transmission event. (It must contain unsampled infected lineages)

#### NOTE: There are many similarities between this code and `sim.outer.tree()` function ####


#' sim.migrations
#' 
#' Annotate and update the EventLogger object with migration events.
#' 
#' @param model: R6 object from MODEL$new()
#' @param eventlog: R6 Eventlogger object populated by sim.outer.tree()
#' 
#' @examples 
#' path <- system.file('extdata', 'structSI.yaml', package='twt')
#' settings <- yaml.load_file(path)
#' mod <- Model$new(settings)
#' 
#' run <- sim.outer.tree(mod)
#' sim.migrations(run)
#' tree$get.migration.events()  # display private data frame
#' tree2 <- sim.inner.tree(mod, tree)
#' 
#' @export
sim.migrations <- function(run) {
  if (!is.element('Run', class(run))) {
    stop("Argument to sim.migrations ")
  }
  
  # unpack Run object
  comps <- run$get.compartments()
  lineages <- run$get.lineages()
  types <- run$get.types()
  eventlog <- run$get.eventlog()
  
  # CompartmentType name for each sampled compartment
  indiv.types <- sapply(unlist(comps), function(a) {a$get.type()$get.name()})
  
  # record population totals and migration rates for all Types
  init.conds <- run$get.initial.conds()
  popn.totals <- init.conds$size
  popn.migration.rates <- .record.migration.rates(types, indiv.types)
  
  if (any(popn.migration.rates > 0)) {
    # record max sampling time of lineages (the time of the most recent sample)
    max.sampling.time <- min(sapply(lineages, function(x) { x$get.sampling.time() }))
    
    # record possible source types for each Compartment type
    possible.source.types <- .record.possible.source.types(types, lineages)
    
    # extract transmission events that occurred from `sim.outer.tree()`
    transmission.events <- eventlog$get.events('transmission')
    
    # generate migration events based on population dynamics up to the max sampling time
    m_events <- .calc.migration.events(
      popn.totals, popn.migration.rates, max.sampling.time, possible.source.types, 
      transmission.events
      )
    
    # store generated migration events into EventLogger, to be resolved later in `sim.inner.tree()`
    eventlog$store.migration.events(m_events)  
  }
}


#' .record.migration.rates
#'
#' Helper function stores migration rates for each CompartmentType 
#' specified by the user.
#'
#' @param types: list of CompartmentType objects
#' @param indiv.types: list of CompartmentType name of type `character` for each 
#' individual Compartment object
#' @return popn.migration.rates = matrix of population migration rates for each 
#' Type to Type comparison (for all possible migration events)
#' @keywords internal
.record.migration.rates <- function(types, indiv.types) {
  popn.migration.rates <- matrix(nrow=length(types),
                                 ncol=length(types),
                                 dimnames=list(names(types), names(types))
                                 )
  for (x in types) {
    for (y in names(types)) {
      mig.rate <- x$get.migration.rate(y)
      if (is.null(mig.rate)) {
        # assume missing rate is zero
        mig.rate <- 0
      }
      popn.migration.rates[x$get.name(), y] <- mig.rate
    }
  }
  
  popn.migration.rates
}


#' .record.possible.source.types
#' 
#' Stores list of possible Types of `source` with recipient types for a 
#' migration event as type character.
#'
#' @param types: list of CompartmentType objects
#' @param lineages: list of Lineage objects
#' @return  named list of possible Types of `source` with 
#' recipient types of type character
#' @keywords internal
.record.possible.source.types <- function(types, lineages) {
  # generate dictionary of different types of source that each recipient Type 
  # could possibly receive a migrating lineage from
  lapply(types, function(t) {
    rates <- t$get.migration.rates()
    names(rates)[rates>0]
  })
}


#' .calc.migration.events
#' 
#' Generates migration events only, based on population dynamics of the MODEL.
#'
#' @param init.conds: initial conditions specified in model settings, i.e., 
#'        originTime and size (per CompartmentType)
#' @param popn.migration.rates: rates of migration between different 
#'        CompartmentTypes
#' @param max.sampling.time: time of the most recent sample as a endpoint for 
#'        this function (type numeric)
#' @param possible.source.types: list of possible Sources that each recipient 
#'        type can receive a migrating lineage from
#' @param transmission.events: 
#'        
#' @return data frame of migration events, each made up of: time, 
#'         recipient Type, and source Type
#' @keywords internal
.calc.migration.events <- function(init.conds, popn.migration.rates, max.sampling.time, 
                                   possible.source.types, transmission.events) {
  # prepare outcome container
  m_events <- data.frame(time=numeric(), r_type=character(), s_type=character(), 
                         v_type=character(), stringsAsFactors = FALSE)
  
  total.outrate <- apply(popn.migration.rates, 1, sum)
  r.types <- colnames(popn.migration.rates)[total.outrate > 0]

  current.time <- init.conds$originTime
  starting.infection <- 1
  
  # check if most recent sample is earlier than the epidemic start time
  if (max.sampling.time > current.time) {
    stop ('Not possible to have Lineage sampling time(s) precede the start time of the "', v.name, '"epidemic. Please set the start time of the epidemic further back in time.')
  }
  
  for (i in transmission.events$time) {
    
    temp.time <- i
    while (temp.time > max.sampling.time) {
      
      # sample source and recipient types for this generated migration evnet
      r <- sample(r.types, 1)
      s <- sample(possible.source.types[[r]], 1)
      
      # calculate waiting time until the next migration event
      num.transmissions.occurred <- length(which(transmission.events$time >= temp.time))
      rate <- popn.migration.rates[s,r] * (num.transmissions.occurred + starting.infection - 1)    # substract by 1, can't receive a migration from itself
      if (rate <= 0) {
        break
      } else {
        delta.t <- rexp(n=1, rate=rate)
        migration.time <- temp.time - delta.t
      }
      
      if (migration.time < max.sampling.time) { 
        break
      } else {
        # store time, source and recipient types of migration event, and virus type
        m_events <- rbind(m_events, list(time=migration.time, r_type=r, s_type=s, v_type=v.name), stringsAsFactors=F)
        temp.time <- migration.time
      }
      
    }
    
  }
    
  
  
  m_events
} 
