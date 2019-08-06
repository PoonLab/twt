# issue 48: Fix implementation of migration rates
# Reason: difficult to conceive and parameterize a migration process in reverse time
# Solution: Annotate outer transmission tree with migration events (going forward in time)
#           When we encounter a migration event during inner tree simulation in reverse time,
#           we check the probability that one of the sampled lineages we are following back
#           in time was a result of the migration.
#           This also means that a sampled compartment can become empty of sampled lineages
#           before we reach its transmission event. (It must contain unsampled infected lineages)

#### NOTE: There are many similarities between this code and `sim.outer.tree()` function #######


sim.migrations <- function(model, eventlog) {
  # Annotate and update the EventLogger object with migration events
  # 
  # @param model: R6 object from MODEL$new()
  # @param eventlog: R6 Eventlogger object populated by sim.outer.tree()
  
  comps <- model$get.compartments()
  compnames <- model$get.names(comps)
  
  types <- model$get.types()
  indiv.types <- sapply(unlist(comps), function(a) {a$get.type()$get.name()})
  
  # record population totals and migration rates for all Types
  popn.totals <- model$get.origin.times()
  popn.migration.rates <- .record.migration.rates(types, indiv.types)

  # record max sampling time of lineages (the time of the most recent sample)
  max.sampling.time <- min(sapply(model$get.lineages(), function(x) {x$get.sampling.time()}))
  
  # record possible source types for each Compartment type
  possible.source.types <- .record.possible.source.types(types, model$get.lineages())
  
  # extract transmission events that occurred from `sim.outer.tree()`
  transmission.events <- eventlog$get.events('transmission')
  
  # generate migration events based on population dynamics up to the max sampling time
  m_events <- .calc.migration.events(popn.totals, popn.migration.rates, max.sampling.time, possible.source.types, transmission.events)
  
  # store generated migration events into EventLogger, to be resolved later in `sim.inner.tree()`
  eventlog$store.migration.events(m_events)
}



.record.migration.rates <- function(types, indiv.types) {
  # helper function stores migration rates for each CompartmentType specified by the user
  #
  # @param types: list of CompartmentType objects
  # @param indiv.types: list of CompartmentType name of type `character` for each individual Compartment object
  # @return popn.migration.rates = matrix of population migration rates for each Type to Type comparison (for all possible migration events)

  popn.migration.rates <- matrix(nrow=length(types),
                                 ncol=length(types),
                                 dimnames=list(names(types), names(types)))
  
  for (x in types) {
    for (y in names(types)) {
      mig.rate <- x$get.migration.rate(y)
      popn.migration.rates[x$get.name(), y] <- mig.rate
    }
  }
  
  popn.migration.rates
}



.record.possible.source.types <- function(types, lineages) {
  # stores list of possible Types of `source` with recipient types for a migration event as type character
  #
  # @param types: list of CompartmentType objects
  # @param lineages: list of Lineage objects
  # @return possible.source.types: names list of possible Types of `source` with recipient types of type character
  
  # generate dictionary of different types of source that each recipient Type could possible receive a migrating lineage from
  typenames <- sapply(types, function(x) {x$get.name()})
  possible.source.types <- as.list(rep(NA, length(typenames)))
  names(possible.source.types) <- typenames
  
  for (t in types) {
    mig.rates <- t$get.migration.rates()
    r.types <- names(mig.rates)
    for (r in 1:length(mig.rates)) {
      if (mig.rates[r] == 0) next
      else possible.source.types[[r.types[r]]] <- t$get.name()
    }
  }
  
  # cleanup
  possible.source.types <- possible.source.types[!is.na(possible.source.types)]
  possible.source.types
  
}



.calc.migration.events <- function(popn.totals, popn.migration.rates, max.sampling.time, possible.source.types, transmission.events) {
  # generates migration events only, based on population dynamics of the MODEL
  #
  # @param popn.totals: totals at time `t=0` of susceptible and infected specific to each CompartmentType
  # @param popn.migration.rates: rates of migration between different CompartmentTypes
  # @param max.sampling.time: time of the most recent sample as a endpoint for this function (type numeric)
  # @param possible.source.types: list of possible Sources that each recipient type can receive a migrating lineage from
  # @return m_events = data frame of migration events, each made up of: time, recipient Type, and source Type
  
  m_events <- data.frame(time=numeric(), r_type=character(), s_type=character(), v_type=character(), stringsAsFactors = FALSE)
  
  for (v in 1:nrow(popn.totals)) {
    
    v.name <- rownames(popn.totals)[v]
    virus <- popn.totals[v,]
    r.types <- names(possible.source.types)
    current.time <- as.numeric(virus['start'])     # current time starts at user given time for each epidemic
    starting.infection <- 1                              # starting infection
    
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
    
  }
  
  m_events
} 
