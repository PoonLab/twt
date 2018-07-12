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
sim.outer.tree <- function(model, eventlog) {
  # simulate transmission events and fix them to the timeline of lineage sampled events
  # @param model = MODEL object
  # @param eventlog = EventLogger object
  # @return eventlog = EventLogger object populated with generated transmission events
  
  # ptm <- proc.time()   # benchmark start time
  
  comps <- model$get.compartments()           
  compnames <- model$get.names(comps)
  
  types <- model$get.types()
  indiv.types <- sapply(unlist(comps), function(a){a$get.type()$get.name()})
  
  # record population totals and transmission rates for all Types
  popn.totals <- model$get.origin.times()
  popn.rates <- .calc.popn.rates(types, indiv.types)
  
  # record max sampling times of lineages for each Compartment
  storage <- .store.initial.samplings(comps, types, model$get.lineages())
  time.bands <- storage$initial.times
  possible.source.types <- storage$s.types
  
  # generate transmission events based on population dynamics and Compartments' initial sampling times
  t_events <- .calc.transmission.events(popn.totals, popn.rates, time.bands, possible.source.types)
  
  # generate unsampled hosts
  last.indiv <- T
  for (t in types) {
    num.unsampled <- length(which(t_events$r_type == t$get.name())) - length(which(indiv.types == t$get.name()))
    if (num.unsampled > 0) {
      if (last.indiv) {
        model$generate.unsampled(num.unsampled+1, t)
        last.indiv <- F
      } else {
        model$generate.unsampled(num.unsampled, t)
      }
    }
  }
  
  # update `time.bands` to include unsampled hosts
  time.bands <- c(time.bands, rep(NA, length(model$get.unsampled.hosts())))
  names(time.bands) <- c(names(time.bands)[nzchar(x=names(time.bands))], model$get.names(model$get.unsampled.hosts()))
  
  # unsampled infected now calculated and generated, set source population
  sources <- c(comps, model$get.unsampled.hosts())
  sources.names <- model$get.names(sources)
  
  sapply(types, function(x) {
    # for each Type, assign transmission times to all infected compartments
    r.indices <- which( sapply(sources, function(y) y$get.type()$get.name() == x$get.name()) )
    r.comps <- sources[r.indices]
    r.init.samplings <- time.bands[ which(names(time.bands) %in% sources.names[r.indices]) ]
    r.events <- t_events[ which(t_events$r_type == x$get.name()), ]
    
    if (nrow(r.events) != 0) {            # no transmission events for case of ie. blood compartment
      .assign.transmission.times(r.comps, r.events, r.init.samplings, x)
    } else {
      sapply(r.comps, function(x) x$set.branching.time(NA))
    }
    
  })
  
  
  # mark1 <- proc.time() - ptm            # benchmark 1 (generation of transmission events)
  # cat("GTE:", round(mark1[['elapsed']],5), "s\n")
  # ptm1 <- proc.time()
  
  # after transmission times are matched with infected Compartments as recipients, now have to assign source Compartments
  # order infection times from most recent to furthest back in time
  new.order <- order(sapply(comps, function(x) x$get.branching.time()))
  comps <- comps[ new.order ]
  compnames <- compnames[ new.order ]
  
  # each assigned transmission time is associated with an event, which determines what TYPE its source is, just not which Compartment in particular
  
  source.popns = source.popns.names <- setNames(vector(length(types), mode="list"), names(types))
  # separate source populations into lists that are as many as the number of distinct Types in the model
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
  
  # assign source compartments for all sampled infected and promoted unsampled infected Compartments
  numActive <- length(comps)
  while (numActive > 1) {
    r_ind_comps <- 1                            # first recipient == most recent infection time == largest list of `sources` --> efficiency                             
    recipient <- comps[[r_ind_comps]]
    r_name <- recipient$get.name()
    r_type <- recipient$get.type()$get.name()
    
    # remove chosen recipient from relevant lists (`comps` and possibly `source.popns`)
    comps[[ r_ind_comps ]] <- NULL
    compnames <- compnames[-r_ind_comps]
    ind_source_popn <- which(source.popns.names[[r_type]] == r_name)
    if (length(ind_source_popn) != 0) {         # TRUE for all iterations except initial iteration case (bc starting recipient already removed)
      source.popns[[r_type]][[ind_source_popn]] <- NULL
      source.popns.names[[r_type]] <- source.popns.names[[r_type]][-ind_source_popn]
    }
    
    # list of possible sources is based off of the `s_type` recorded in the event associated with the transmission time in master copy `t_events`
    if (is.na(recipient$get.branching.time())) {
      # root case, "final" resolved transmission event (aka first recorded transmission, furthest back in time)
      break
      
    } else {
      s_type <- t_events[ which(t_events$time == recipient$get.branching.time()), 's_type']
      # refine list of sources previously separated by Type also by earlier branching times than recipient's branching time
      refined.list.sources <- sapply(source.popns[[s_type]], function(s) {
        if (is.na(s$get.branching.time())) s
        else if (s$get.branching.time() > recipient$get.branching.time()) s
        else NULL
      })
      refined.list.sources[sapply(refined.list.sources, is.null)] <- NULL   # cleanup
      
      # select source from a list of sources previously separated by Type
      s_ind_s_popn <- sample.int(length(refined.list.sources), 1)
      source <- refined.list.sources[[s_ind_s_popn]]
      s_name <- source$get.name()
      
      eventlog$add.event('transmission', recipient$get.branching.time(),NA, NA, r_name, s_name)
    }
    
    # if source is an unsampled infected Compartment, now holds a sampled lineage we care about (promote us_comp)
    if (source$is.unsampled() && source$get.name() %in% compnames == FALSE) {
      # add us_comp to list `comps` (once first a source, can now be a recipient)
      comps[[length(comps)+1]] <- source
      compnames[[length(compnames)+1]] <- s_name
    }
    
    # update recipient object `source` attr
    recipient$set.source(source)
    
    # update number of active compartments (excludes recipients, includes promoted us_comps)
    new.order <- order(sapply(comps, function(x) x$get.branching.time()))
    comps <- comps[ new.order ]
    compnames <- compnames[ new.order ]
    numActive <- length(comps)
  }
  
  
  # mark2 <- proc.time() - ptm1               # benchmark 2 (assignment of sources)
  # cat("AST:", round(mark2[['elapsed']],5), "s\n")
  # total <- proc.time() - ptm
  # cat("Total:", round(total[['elapsed']],5), "s\n")
    
  eventlog
}




.calc.popn.rates <- function(types, indiv.types) {
  # stores population rates for each CompartmentType specified by the user
  # @param types = all CompartmentTypes
  # @param indiv.types = list of individual Types for each Compartment
  # @return population transmission rates matrix
  
  popn.rates <- matrix(nrow=length(types),                                # source types
                       ncol=length(types),                                # recipient types
                       dimnames=list(names(types), names(types)))
  
  for (x in types) {                                                      # for each CompartmentType:
    # I <- length(which(indiv.types == x$get.name())) + x$get.unsampled()   # 2. store number of active sampled compartments AND unsampled infected hosts at time t=0 (I)
    # S <- x$get.susceptible()                                              # 3. store number of susceptibles at time t=0 (I)
    # popn.totals[x$get.name(),] <- c(I,S)
    
    for (y in names(types)) {
      rate <- x$get.branching.rate(y)                                     # store instrinsic transmission rates for all typeA -> typeB pairs
      popn.rates[x$get.name(), y] <- rate
    }
  }
  
  popn.rates
}




.store.initial.samplings <- function(infected, types, lineages) {
  # stores first sampling time of a Lineage for each sampled infected Compartment
  # @param infected = list of infected Compartment objects
  # @param types = list of CompartmentType objects
  # @param lineages = list of Lineage objects
  # @return list of possible Types of `source` and list of first Lineage sampling times for each Compartment
  
  # generate dictionary of different types of source that each recipient Type could possibly receive a transmission from
  typenames <- sapply(types, function(x) x$get.name())
  possibleSourceTypes <- as.list(rep(NA, length(typenames)))  
  names(possibleSourceTypes) <- typenames
  for (t in types) {
    rates <- t$get.branching.rates()
    r.types <- names(rates)
    for (r in 1:length(rates)) {
      if (rates[r] == 0) next
      else possibleSourceTypes[[r.types[r]]] <- t$get.name()
    }
  }
  possibleSourceTypes <- possibleSourceTypes[!is.na(possibleSourceTypes)] # cleanup

  
  time.bands <- vector()             # vector of initial (earliest) sampling times for each Compartment
  lineage.locations <- sapply(lineages, function(x) {x$get.location()$get.name()})
  
  for (i in infected) {
    # filtered to remove all source -> recipient pairs with branching rates of 0
    recipientType <- i$get.type()$get.name()
    
    if (recipientType %in% names(possibleSourceTypes)) {
      i.native.lineages <- which(lineage.locations == i$get.name())
      single.i.lineages <- lineages[ i.native.lineages ]
      single.i.sampling.times <- sapply(single.i.lineages, function(b) {
        b$get.sampling.time()
      })
      
      # store initial sampling time for this recipient Compartment
      time.bands <- c(time.bands, max(single.i.sampling.times))
      names(time.bands) <- c(names(time.bands)[nzchar(x=names(time.bands))], i$get.name())
      
    } else {
      # this compartment will never be a recipient (ie. example3.yaml 'blood' compartment)
      next
      
    }
  }
  
  list(s.types=possibleSourceTypes, initial.times=time.bands)
}




.calc.transmission.events <- function(popn.totals, popn.rates, init.samplings, possible.sources) {
  # generates transmission events only, based on population dynamics of the MODEL
  # @param popn.totals = totals at time `t=0` of susceptible and infected specific to each CompartmentType
  # @param popn.rates = rates of transmission between different CompartmentTypes
  # @param init.samplings = list of first sampling times for each Compartment
  # @param possible.sources = list of possible Sources that each recipient Type can receive a transmission from 
  # @return t_events = data frame of transmission events, each made up of: time, recipient Type, and source Type 
  
  t_events <- data.frame(time=numeric(), r_type=character(), s_type=character(), v_type=character())
  maxAttempts <- 5 
  
  for (attempt in 1:maxAttempts) {
    
    for (v in 1:nrow(popn.totals)) {
      
      v.name <- rownames(popn.totals)[v]
      virus <- popn.totals[v,]
      r.types <- names(virus)[-1]
      current.time <- as.numeric(virus['start'])          # current time starts at user given time for each epidemic
      num.infected <- 1                                   # starting infection
      
      # check sampling times with epidemic start time
      if (any(init.samplings > current.time)) {
        stop ('Not possible to have Compartment initial sampling time(s) precede the start time of the "', v.name, '" epidemic. Please set the start time of the epidemic further back in time.')
      }
      
      while (current.time > min(init.samplings) && all(virus != 1)) {
        # calculate total waiting time
        r <- sample(r.types, 1)
        s <- sample(possible.sources[[r]], 1)
        
        # total waiting time = exp(-beta * (S-1) * I)
        rate <- popn.rates[r,s] * (virus[s]-1) * num.infected
        delta.t <- rexp(n=1, rate=rate)
        current.time <- current.time - delta.t
        
        if (current.time < min(init.samplings)) { break }
        
        # store time, source and recipient types of transmission, and virus type
        t_events <- rbind(t_events, list(time=current.time, r_type=r, s_type=s, v_type=v.name), stringsAsFactors=F)
        
        # update counts of total population
        virus[s] <- virus[s] - 1
        num.infected <- num.infected + 1
        
      }
    
    }
    
    # check if viable # of transmission times @ each unique user-given Compartment `sampling.time`
    checks <- sapply(unique(init.samplings), function(x) {
      if (length(init.samplings==x) > length(t_events$time > x)) {
        # not enough viable transmission times at this `sampling.time`
        # regenerate transmission times up to maxAttempts
        FALSE
      } else {
        TRUE
      }
    })
    
    if (any(checks) == FALSE) {
      attempt <- attempt + 1
    } else if (attempt > maxAttempts) {
      stop ('Transmission times generated invalid matches to given `sampling.times` of Compartments. Please change origin time of the "', v.name, '" epidemic or modify transmission rates.')
    } else {
      break
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
    
    if (length(possibleTimes) < length(recipients)) {  # for case where one of these recipients will be labeled as the root node (ie. no unsampled hosts)
      sampledTimes <- sample(possibleTimes, size=length(recipients)-1, prob=weights, replace=F)
    } else {
      sampledTimes <- sample(possibleTimes, size=length(recipients), prob=weights, replace=F)
    }
    
    for (x in 1:length(sampledTimes)) {
      # now assign `sampledTimes` to `recipients`, uniformly distributed
      ind <- sample.int(length(recipients), 1)
      t.time <- sampledTimes[x]
      recipients[[ind]]$set.branching.time(t.time)
      
      recipients <- recipients[-ind]
      events <- events[ -which(events$time == t.time), ]                     # remove transmission event from events
    }
    
    if (length(recipients) != 0) sapply(recipients, function(x) x$set.branching.time(NA))
    
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
