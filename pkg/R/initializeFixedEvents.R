init.fixed.samplings <- function(model) {
  # treeswithintrees/Wiki/Simulation Pseudocode step 2
  # retrieve sampling time and populate tip labels / times in ape::phylo object (building it tips up)
  # @param model = MODEL object

  # add lineage sampling events from Lineage objects
  lineages <- model$get.lineages()

  list(
    # store label w/ corresponding tip height in new ape::phylo object (not casted into `phylo` yet)
    tip.label = sapply(lineages, function(x) x$get.name()),

    # only used for calculating edge length
    tip.height = sapply(lineages, function(x) x$get.sampling.time())
  )
}



# after the objects are generated from user inputs, we need to initialize the list of fixed events
# Case 1 : User provides a host transmission tree
.to.eventlog <- function(newick) {
  # function converts an ape::phylo tree object into transmission events stored in a NEW EventLogger object
  # @param newick = Newick string in mode character
  # @return e = EventLogger object initialized with list of fixed transmission events
  
  e <- EventLogger$new()
  phy <- read.tree(text=newick)
  
  if (is.null(phy$node.label)) {
    # in the future, may be able to generate unique name for internal node
    stop('Node labels must be present in host transmission tree.')
  } 
  
  sapply(1:nrow(phy$edge), function(x) {
    s_ind <- phy$edge[x,1]
    r_ind <- phy$edge[x,2]
    
    sourceLabel <- phy$node.label[s_ind]
    
    if (r_ind <= length(phy$tip.label)) {
      recipientLabel <- phy$tip.label[r_ind]
    } else {
      recipientLabel <- phy$node.label[r_ind]
    }
    
    branching.time <- phy$edge.length[x]
    e$add.event('transmission', branching.time, NA, NA, recipientLabel, sourceLabel)
  })
  
  e
}



# Case 2 : User manually inputs a host transmission tree into YAML format under Compartments header
init.branching.events <- function(model, eventlog) {
  # @param model = MODEL object
  # @param eventlog = EventLogger object
  # @return eventlog = EventLogger object initialized with list of fixed transmission events

  # if the user input includes a tree (host tree) then add transmission events
  comps <- model$get.compartments()
  lineages <- model$get.lineages()

  transmissions <- sapply(comps, function(x) {
    branching.time <- x$get.branching.time()
      
    if (is.R6(x$get.source())) {
      source <- x$get.source()$get.name()
      xLin <- sapply(lineages, function(y){which(y$get.location()$get.name() == x$get.name())})
      lineage <- lineages[[ which(xLin == 1) ]]$get.name()
    } else {
      source <- x$get.source()
      lineage <- NA
    }

    # add transmission event to EventLogger object
    eventlog$add.event('transmission',  branching.time, lineage, NA, x$get.name(), source)
  })
  
  eventlog
}



# Case 3: no host tree provided, transmission events need to be generated
generate.transmission.events <- function(model, eventlog) {
  # simulate transmission events and fix them to the timeline of lineage sampled events
  # @param model = MODEL object
  # @param eventlog = EventLogger object
  # @return eventlog = EventLogger object populated with generated transmission events
  
  comps <- model$get.compartments()           
  compnames <- model$get.names(comps)
  sources <- c(comps, model$get.unsampled.hosts())
  sources.names <- model$get.names(sources)
  
  types <- model$get.types()
  indiv.types <- sapply(unlist(comps), function(a){a$get.type()$get.name()})
  storage <- .calc.popn.totals.rates(types, indiv.types)
  
  popn.totals <- storage$totals
  popn.rates <- storage$rates
  
  
  # record max sampling times of lineages for each Compartment
  storage <- .store.initial.samplings(sources, types, model$get.lineages())
  possible.source.types <- storage$s.types
  time.bands <- storage$initial.times
  
  t_events <- .calc.transmission.events(popn.totals, popn.rates, time.bands, possible.source.types)
  
  sapply(types, function(x) {
    r.indices <- which( sapply(sources, function(y) y$get.type()$get.name() == x$get.name()) )
    r.comps <- sources[r.indices]
    r.init.samplings <- time.bands[ which(names(time.bands) %in% sources.names[r.indices]) ]
    r.events <- t_events[ which(t_events$r_type == x$get.name()), ]
    
    .assign.transmission.times(r.comps, r.events, r.init.samplings, x)
  })

  
  # after transmission times are matched with infected Compartments as recipients, now have to assign source Compartments
  # order infection times from most recent to furthest back in time
  comps <- comps[ order(sapply(comps, function(x) x$get.branching.time())) ] 
  
  # note that each matched transmission time is associated with an event, which determines what TYPE the source is, just not which in particular
  # will need to separate source populations into lists that are as many as the number of distinct Types in the model
  source.popns = source.popns.names <- setNames(vector(length(types), mode="list"), names(types))
  for (x in 1:length(names(source.popns))) {
    specific.type <- names(source.popns)[x]
    s.pop.by.type <- sapply(sources, function(y) {
      if (y$get.type()$get.name() == specific.type) {
        if (is.na(y$get.branching.time())) {y}    # root no branching time
        else if (y$get.branching.time() > comps[[1]]$get.branching.time()) {y}
      } 
    })
    s.pop.by.type[sapply(s.pop.by.type, is.null)] <- NULL    # cleanup
    s.pop.by.type.names <- sapply(s.pop.by.type, function(z) z$get.name())
    source.popns[[specific.type]] <- s.pop.by.type
    source.popns.names[[specific.type]] <- s.pop.by.type.names
  }
  
  numActive <- length(comps)
  while (numActive > 1) {
  
    # start at first recipient (will give largest list of `sources`) --> for efficiency
    r_ind_comps <- 1                                        
    recipient <- comps[[r_ind_comps]]
    r_name <- recipient$get.name()
    r_type <- recipient$get.type()$get.name()
    
    # remove chosen recipient from relevant lists
    comps[[ r_ind_comps ]] <- NULL
    compnames <- compnames[-r_ind_comps]
    ind_source_popn <- which(source.popns.names[[r_type]] == r_name)
    if (length(ind_source_popn) != 0) {
      source.popns[[r_type]][[ind_source_popn]] <- NULL
      source.popns.names[[r_type]] <- source.popns.names[[r_type]][-ind_source_popn]
    }
    
    # list of possible sources is based off of the `s_type` recorded in the event associated with the transmission time in master copy `t_events`
    if (is.na(recipient$get.branching.time())) {
      s_name <- NA
      eventlog$add.event('transmission', recipient$get.branching.time(),NA, NA, r_name, s_name)
      break
      
    } else {
      s_type <- t_events[ which(t_events$time == recipient$get.branching.time()), 's_type']
      list.sources <- source.popns.names[[s_type]]
      
      s_ind_s_popn <- sample.int(length(list.sources), 1)
      s_name <- list.sources[s_ind_s_popn]
      source <- source.popns[[s_type]][[s_ind_s_popn]]
      
      eventlog$add.event('transmission', recipient$get.branching.time(),NA, NA, r_name, s_name)
    }
    
    # if source is an unsampled infected Compartment, now holds a sampled lineage we care about (promote us_comp)
    if (source$is.unsampled()) {
      # add us_comp to list `comps` (once first a source, can now be a recipient)
      comps[[length(comps)+1]] <- source
      compnames[[length(compnames)+1]] <- s_name
    }
    
    # update recipient object `source` attr
    recipient$set.source(source)
    
    # update number of active compartments (excludes recipients, includes promoted us_comps)
    comps <- comps[ order(sapply(comps, function(x) x$get.branching.time())) ]
    numActive <- length(comps)
  }

  eventlog
}




.calc.popn.totals.rates <- function(types, indiv.types) {
  # stores population totals and rates for each CompartmentType specified by the user
  # @param types = all CompartmentTypes
  # @param indiv.types = list of individual Types for each Compartment
  # @return list of population totals matrix and population transmission rates matrix
  
  popn.totals <- matrix(nrow=length(types), ncol=2,
                        dimnames=list(names(types), c('I', 'S')))
  popn.rates <- matrix(nrow=length(types),                                # source types
                       ncol=length(types),                                # recipient types
                       dimnames=list(names(types), names(types)))
  
  for (x in types) {                                                      # for each CompartmentType:
    I <- length(which(indiv.types == x$get.name())) + x$get.unsampled()   # 2. store number of active sampled compartments AND unsampled infected hosts at time t=0 (I)
    S <- x$get.susceptible()                                              # 3. store number of susceptibles at time t=0 (I)
    popn.totals[x$get.name(),] <- c(I,S)
    
    for (y in names(types)) {
      rate <- x$get.branching.rate(y)                                     # store instrinsic transmission rates for all typeA -> typeB pairs
      popn.rates[x$get.name(), y] <- rate
    }
  }
  
  list(totals=popn.totals, rates=popn.rates)
}




.store.initial.samplings <- function(infected, types, lineages) {
  # stores first sampling time of a Lineage for each sampled infected Compartment
  # @param infected = list of infected Compartment objects
  # @param types = list of CompartmentType objects
  # @param lineages = list of Lineage objects
  # @return list of possible Types of `source` and list of first Lineage sampling times for each Compartment
  
  possibleSourceTypes <- list()  # list of different types of source that each recipient could possibly receive a transmission from
  time.bands <- vector()         # vector of maximum sampling times for each Compartment
  
  lineage.locations <- sapply(lineages, function(x) {x$get.location()$get.name()})
  
  for (infected in infected) {
    # for loop generates a comprehensive dictionary of all possibilities of a source with a given recipient Compartment
    # filtered first within the for loop to remove all source -> recipient pairs with branching rates of 0
    # filtered again at while loop to skip unsampled compartments as potential recipients unless they were first a source
    recipientType <- infected$get.type()$get.name()
    recipientRates <- sapply(types, function(a) {
      if (a$get.branching.rate(recipientType) == 0) { NULL } 
      else { a$get.branching.rate(recipientType) }
    })
    recipientRates[sapply(recipientRates, is.null)] <- NULL      # remove source -> recipient pairs w/ branching rates of 0
    
    if (length(recipientRates) == 0) {
      # means that this compartment will never be a recipient (ie. example3.yaml 'blood' compartment)
      next
    } else {
      # check if infected is a us_comp (us_comps have no sampled lineages native to their compartment)
      infected.native.lineages <- which(lineage.locations == infected$get.name())
      if (length(infected.native.lineages) == 0) {
        # us_comp is not allowed to be a recipient without first being a source
        # placeholder 'NA' to be filtered out later in while loop
        single.infected.sampling.times <- NA   
      } else {
        single.infected.lineages <- lineages[ infected.native.lineages ]
        single.infected.sampling.times <- sapply(single.infected.lineages, function(b) {
          b$get.sampling.time()
        })
      }
      
      time.bands <- c(time.bands, max(single.infected.sampling.times))
      names(time.bands) <- c(names(time.bands)[nzchar(x=names(time.bands))], infected$get.name())
      possibleSourceTypes <- c(possibleSourceTypes, list(names(recipientRates)))
      names(possibleSourceTypes)[[length(possibleSourceTypes)]] <- recipientType
    }
  }
  
  list(s.types=possibleSourceTypes, initial.times=time.bands)
}




.calc.transmission.events <- function(popn.totals, popn.rates, init.samplings, possible.sources) {
  # creates transmission events only, based on population dynamics of the MODEL
  # @param popn.totals = totals at time `t=0` of susceptible and infected specific to each CompartmentType
  # @param popn.rates = rates of transmission between different CompartmentTypes
  # @param init.samplings = list of first sampling times for each Compartment
  # @param possible.sources = list of possible Sources that a Type can receive a transmission from 
  # @return t_events = data frame of transmission events, each made up of: time, recipient Type, and source Type 
  
  # the start time of the simulation is the time where at least 2 sampled infected compartments are active
  current.time <- as.numeric(init.samplings[order(init.samplings)[2]])
  t_events <- data.frame(time=numeric(), r_type=character(), s_type=character())
  
  while (sum(popn.totals[,'I']) > 1) {
    # filter population to determine which of these active compartments can be a recipient
    # based on if all possible intrinsic rates are 0, and the max sampling times of the active compartments
    qualified.sampled.recipients <- which(init.samplings <= current.time)
    
    # retrieve transmission types from dictionary of possible.sources for each qualified sampled recipient
    possible.transm <- possible.sources[ qualified.sampled.recipients ]
    
    # calc transmission rates among all source-recipient pairings of CompartmentTypes
    indiv.rates <- sapply(names(types), function(x) {
      sapply(names(types), function(y) {
        pairRate <- popn.rates[x,y] * (popn.totals[y,'S'] + 1) * (popn.totals[x,'I'])
        
        qualified.r <- which(names(possible.transm) == y)
        if (length(qualified.r) ==0) {
          nPairs <- 0
        } else {
          qualified.sr <- which(sapply(qualified.r, function(z) {
            x %in% possible.transm[z]
          }))
          if (length(qualified.sr) == 0) {
            nPairs <- 0
          } else {
            nPairs <- length(qualified.sr)                   # the number of pairs that have this transmission type
          }
        }
        
        pairRate * nPairs
      })
    })
    
    # total rate of ANY transmission event occurring is the weightes sum of these rates in the dictionary
    total.rate <- sum(indiv.rates)
    # sample waiting time & update the current time to the new upper bound in time (waiting time)
    waiting.time <- current.time + rexp(n=1, rate=total.rate)
    current.time <- waiting.time
    
    if (length(which(init.samplings <= current.time)) > length(qualified.sampled.recipients)) {
      # check if the waiting time exceeds any sampling time within the sampled compartments previously not qualifying as a recipient
      # re-start the filtering to include new qualified sampled infected recipients
      next
    } else {
      # randomly choose a source and recipient pair, uniformly distributed
      ind <- sample(1:length(possible.transm), 1)
      r_type <- names(possible.transm)[[ind]]
      s_type <- sample(possible.transm[[ind]], 1)
      
      # store the time, and the source and recipient types of transmission
      t_events <- rbind(t_events, list(time=current.time, r_type=r_type, s_type=s_type), stringsAsFactors=F)
      
      # remove r_type w/ all it's possible source types from list `possible.transm` (this Infected can no longer be a recipient)
      possible.transm[[ind]] <- NULL
      
      # update counts of total population
      popn.totals[r_type, 'I'] <- popn.totals[r_type, 'I'] - 1
      popn.totals[r_type, 'S'] <- popn.totals[r_type, 'S'] + 1
    }
  }
  
  t_events
}




.assign.transmission.times <- function(infected, events, initial.samplings, type) {
  # assignment of transmission times specific to each type
  # @param infected = list of Compartment objects to be assigned transmission times
  # @param events = list of possible transmission events for the infected
  # @param initial.samplings = list of respective initial sampling times for each Compartment
  # @param type = CompartmentType object
  
  infected.names <- sapply(infected, function(x) x$get.name())
  i.times <- initial.samplings[!is.na(initial.samplings)]                     # separate sampled infected comps from 
  u.times <- initial.samplings[is.na(initial.samplings)]                      # unsampled infected comps
  
  while (length(i.times) > 1) {
    # for sampled infected Compartments, assignments based off of Type-specific events and Type-specific waiting time distribution
    
    earliest.time <- max(i.times, na.rm=T)                                    # pick Compartment with earliest sampling time (furthest back in time)
    recipients.names <- names(i.times)[which(i.times==earliest.time)]         # group all infected with the same sampling time together
    
    recipients.inds <- sapply(recipients.names, function(x) which(infected.names == x))
    recipients <- infected[recipients.inds]
    
    # retrieve transmission times with a recipient type equal to the general type of `recipients`
    # weight each transmission time according to the `wait.time.distr` provided for this type
    possibleTimes <- events$time[ which(events$r_type == type$get.name()) ]
    weights <- sapply(possibleTimes, function(x) eval(parse(text=type$get.wait.time.distr())) )
    
    sampledTimes <- sample(possibleTimes, size=length(recipients), prob=weights, replace=F)
    
    for (x in recipients) {
      # now assign `sampledTimes` to `recipients`, uniformly distributed
      ind <- sample.int(length(sampledTimes), 1)
      t.time <- sampledTimes[ind]
      x$set.branching.time(t.time)
      
      sampledTimes <- sampledTimes[-ind]                                     # remove transmission time from sampledTimes
      events <- events[ -which(events$time == t.time), ]                     # remove transmission event from events
    }
    
    i.times <- i.times[ -which(i.times==earliest.time) ] 
  }
  
  while (length(u.times) > 1) {
    # for unsampled infected Compartments, assignments based off of Type-specific events and uniform transmission time distribution
    u.ind <- sample.int(length(u.times), 1)
    u.name <- names(u.times)[u.ind]
    u.comp <- infected[[ which(infected.names == u.name) ]]
    
    u.event <- sample(nrow(events), 1)
    t.time <- events[u.event, 'time']
    u.comp$set.branching.time(t.time)
    
    u.times <- u.times[-u.ind]
    events <- events[-u.event,]
  }
  
}

