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
  source.popn <- c(comps, model$get.unsampled.hosts())
  source.popn.names <- model$get.names(source.popn)
  types <- model$get.types()
  typenames <- model$get.names(types)
  
  indiv.types <- sapply(unlist(comps), function(a){a$get.type()$get.name()})
  popn.totals <- matrix(nrow=length(types),
                        ncol=2,
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
  
  
  
  # record max sampling times of lineages for each Compartment
  lineage.locations <- sapply(model$get.lineages(), function(x) {x$get.location()$get.name()})
  possibleSourceTypes <- list()  # list of different types of source that each recipient could possibly receive a transmission from
  time.bands <- vector()         # vector of maximum sampling times for each Compartment
  
  for (comp in source.popn) {
    # for loop generates a comprehensive dictionary of all possibilities of a source with a given recipient Compartment
    # filtered first within the for loop to remove all source -> recipient pairs with branching rates of 0
    # filtered again at while loop to skip unsampled compartments as potential recipients unless they were first a source
    recipientType <- comp$get.type()$get.name()
    recipientRates <- sapply(model$get.types(), function(a) {
      if (a$get.branching.rate(recipientType) == 0) { NULL } 
      else { a$get.branching.rate(recipientType) }
    })
    recipientRates[sapply(recipientRates, is.null)] <- NULL      # remove source -> recipient pairs w/ branching rates of 0
    
    if (length(recipientRates) == 0) {
      # means that this compartment will never be a recipient (ie. example3.yaml 'blood' compartment)
      next
    } else {
      # check if comp is a us_comp (us_comps have no sampled lineages native to their compartment)
      comp.native.lineages <- which(lineage.locations == comp$get.name())
      if (length(comp.native.lineages) == 0) {
        # us_comp is not allowed to be a recipient without first being a source
        # placeholder 'NA' to be filtered out later in while loop
        single.comp.sampling.times <- NA   
      } else {
        single.comp.lineages <- model$get.lineages()[ comp.native.lineages ]
        single.comp.sampling.times <- sapply(single.comp.lineages, function(b) {
          b$get.sampling.time()
        })
      }
      
      time.bands <- c(time.bands, max(single.comp.sampling.times))
      names(time.bands) <- c(names(time.bands)[nzchar(x=names(time.bands))], comp$get.name())
      possibleSourceTypes <- c(possibleSourceTypes, list(names(recipientRates)))
      names(possibleSourceTypes)[[length(possibleSourceTypes)]] <- recipientType
    }
  }
  
  
  
  t_events <- data.frame(time=numeric(), r_type=character(), s_type=character())
  
  # the start time of the simulation is the time where at least 2 sampled infected compartments are active
  current.time <- as.numeric(time.bands[order(time.bands)[2]])
  
  while (sum(popn.totals[,'I']) > 1) {
    # filter population to determine which of these active compartments can be a recipient
    # based on if all possible intrinsic rates are 0, and the max sampling times of the active compartments
    qualified.sampled.recipients <- which(time.bands <= current.time)
    
    # retrieve transmission types from dictionary of possibleSourceTypes for each qualified sampled recipient
    possibleTransmissions <- possibleSourceTypes[ qualified.sampled.recipients ]
    
    # calc transmission rates among all source-recipient pairings of CompartmentTypes
    indiv.rates <- sapply(typenames, function(x) {
      sapply(typenames, function(y) {
        pairRate <- popn.rates[x,y] * (popn.totals[y,'S'] + 1) * (popn.totals[x,'I'])
        
        qualified.r <- which(names(possibleTransmissions) == y)
        if (length(qualified.r) ==0) {
          nPairs <- 0
        } else {
          qualified.sr <- which(sapply(qualified.r, function(z) {
            x %in% possibleTransmissions[z]
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
    
    if (length(which(time.bands <= current.time)) > length(qualified.sampled.recipients)) {
      # check if the waiting time exceeds any sampling time within the sampled compartments previously not qualifying as a recipient
      # re-start the filtering to include new qualified sampled infected recipients
      next
    } else {
      # randomly choose a source and recipient pair, uniformly distributed
      ind <- sample(1:length(possibleTransmissions), 1)
      r_type <- names(possibleTransmissions)[[ind]]
      s_type <- sample(possibleTransmissions[[ind]], 1)
      
      # store the time, and the source and recipient types of transmission
      t_events <- rbind(t_events, list(time=current.time, r_type=r_type, s_type=s_type), stringsAsFactors=F)
      
      # remove r_type w/ all it's possible source types from list `possibleTransmissions` (this Infected can no longer be a recipient)
      possibleTransmissions[[ind]] <- NULL
      
      # update counts of total population
      popn.totals[r_type, 'I'] <- popn.totals[r_type, 'I'] - 1
      popn.totals[r_type, 'S'] <- popn.totals[r_type, 'S'] + 1
    }
  }
  

    
  sapply(types, function(x) {
    r.indices <- which( sapply(source.popn, function(y) y$get.type()$get.name() == x$get.name()) )
    r.comps <- source.popn[r.indices]
    initial.samplings <- time.bands[ which(names(time.bands) %in% source.popn.names[r.indices]) ]
    
    r.events <- t_events[ which(t_events$r_type == x$get.name()), ]
    
    .assign.transmission.times(r.comps, r.events, initial.samplings, x)
  })

  
  
  # after transmission times are matched with infected Compartments as recipients, now have to assign source Compartments
 
  # order infection times from most recent to furthest back in time
  comps <- comps[ order(sapply(comps, function(x) x$get.branching.time())) ]
  
  # note that each matched transmission time is associated with an event, which determines what TYPE the source is, just not which in particular
  # will need to separate source populations into lists that are as many as the number of distinct Types in the model
  source.popns = source.popns.names <- setNames(vector(length(types), mode="list"), typenames)
  for (x in names(source.popns)) {
    s.pop.by.type <- sapply(source.popn, function(y) {
      if (y$get.type()$get.name() == x) {
        if (is.na(y$get.branching.time())) {y}    # root no branching time
        else if (y$get.branching.time() > comps[[1]]$get.branching.time()) {y}
      } 
    })
    s.pop.by.type[sapply(s.pop.by.type, is.null)] <- NULL    # cleanup
    s.pop.by.type.names <- sapply(s.pop.by.type, function(z) z$get.name())
    source.popns[[x]] <- s.pop.by.type
    source.popns.names[[x]] <- s.pop.by.type.names
  }
  
  numActive <- length(comps)
  while (numActive > 1) {
    # start at first recipient (will give largest `source.popn`)
    # `source.popns` can be updated efficiently starting at the largest `source.popns` initially, then cut down little by little each time
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
      source.popns.names <- source.popns.names[[r_type]][-ind_source_popn]
    }
    
    # list of possible sources is based off of the `s_type` recorded in the event associated with the transmission time in master copy `t_events`
    s_type <- t_events[ which(t_events$time == recipient$get.branching.time()), 's_type']
    list.sources <- source.popns.names[[s_type]]
    
    s_ind_s_popn <- sample.int(length(list.sources), 1)
    s_name <- list.sources[s_ind_s_popn]
    source <- source.popns[[s_type]][[s_ind_s_popn]]

    
    eventlog$add.event('transmission', recipient$get.branching.time(),NA, NA, r_name, s_name)
   

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




.assign.transmission.times <- function(compartments, eventsList, initialSamplings, compartmentType) {
  # this function is to be applied to each specific CompartmentType
  # for sampled infected compartments, based off of Type-specific eventsList and Type-specific waiting time distribution
  # for unsampled infected compartments, based off of Type-specific eventsList and uniform transmission time distribution
  
  compartments.names <- sapply(compartments, function(x) x$get.name())
  i.times <- initialSamplings[!is.na(initialSamplings)]                     # separate sampled infected comps from 
  u.times <- initialSamplings[is.na(initialSamplings)]                      # unsampled infected comps
  
  while (length(i.times) > 1) {
    earliest.time <- max(i.times, na.rm=T)                                  # pick Compartment with earliest sampling time (furthest back in time)
    recipients.names <- names(i.times)[which(i.times==earliest.time)]       # group all Compartments with the same sampling time together
    
    recipients.inds <- sapply(recipients.names, function(x) which(compnames == x))
    recipients <- comps[recipients.inds]
    
    # retrieve transmission times with a recipient type equal to the general CompartmentType of `recipients`
    # weight each transmission time according to the `wait.time.distr` provided for this CompartmentType
    possibleTimes <- eventsList$time[ which(eventsList$r_type == compartmentType$get.name()) ]
    weights <- sapply(possibleTimes, function(x) eval(parse(text=compartmentType$get.wait.time.distr())) )
    
    sampledTimes <- sample(possibleTimes, size=length(recipients), prob=weights, replace=F)
    
    # now assign `sampledTimes` to `recipients`, uniformly distributed
    for (x in recipients) {
      ind <- sample.int(length(sampledTimes), 1)
      t.time <- sampledTimes[ind]
      x$set.branching.time(t.time)
      
      sampledTimes <- sampledTimes[-ind]                                     # remove transmission time from sampledTimes
      eventsList <- eventsList[ -which(eventsList$time == t.time), ]         # remove transmission event from eventsList
    }
    
    i.times <- i.times[ -which(i.times==earliest.time) ] 
  }
  
  # for unsampled infected recipients, randomly assign remaining transmission times, uniformly distributed
  while (length(u.times) > 1) {
    u.ind <- sample.int(length(u.times), 1)
    u.name <- names(u.times)[u.ind]
    u.comp <- compartments[[ which(compartments.names == u.name) ]]
    
    u.event <- sample(nrow(eventsList), 1)
    t.time <- eventsList[u.event, 'time']
    u.comp$set.branching.time(t.time)
    
    u.times <- u.times[-u.ind]
    eventsList <- eventsList[-u.event,]
  }
  
}

  


.to.transmission.tree <- function(eventlog) {
  # function converts the transmission events stored in an EventLogger object into a transmission tree
  # @param eventlog = EventLogger object
  # @return phy = ape::phylo object

  t_events <- eventlog$get.events('transmission', cumulative=FALSE)
  tips <- unlist(setdiff(t_events$compartment1, t_events$compartment2))
  internals <- unlist(intersect(t_events$compartment1, t_events$compartment2))
  root <- unlist(setdiff(t_events$compartment2, t_events$compartment1))
  
  # initializing attributes of an ape::phylo object
  tip.label <- vector()
  edge.length <- vector()
  Nnode <- length(unique(t_events$compartment2))
  edge <- matrix(nrow=nrow(t_events), ncol=2)
  node.label <- vector()
  
  tip.no <- 1
  root.no <- length(tips) + 1
  node.no <- root.no + 1
  
  for (x in 1:nrow(t_events)) {
    source <- t_events[x,]$compartment2
    recipient <- t_events[x,]$compartment1
    
    if (recipient %in% tips) {
      recipient.ind <- tip.no
      tip.label[recipient.ind] <- recipient
      tip.no <- tip.no + 1
    } else if (recipient %in% node.label) {
      recipient.ind <- which(sapply(node.label, function(y) {y == recipient}) == T)
    } else {
      recipient.ind <- node.no
      node.label[recipient.ind] <- recipient
      node.no <- node.no + 1
    }
    
    if (source %in% root) {
      source.ind <- root.no
      node.label[source.ind] <- source
    } else if (source %in% node.label) {
      source.ind <- which(sapply(node.label, function(y) {y == source}) == T)
    } else {
      source.ind <- node.no
      node.label[source.ind] <- source
      node.no <- node.no + 1
    }
    
    #edge[2*x-1,] <- as.numeric(c(source.ind, source.ind))      # source --> source
    edge[x,] <- as.numeric(c(source.ind, recipient.ind))     # source --> recipient
    #edge.length[2*x-1] = 
    edge.length[x] <- t_events[x,]$time
  }
  
  phy <- list(tip.label=tip.label, Nnode=Nnode, edge.length=as.numeric(edge.length), edge=edge, node.label=node.label[!is.na(node.label)])
  attr(phy, 'class') <- 'phylo'
  attr(phy, 'order') <- 'cladewise'
  phy
}

