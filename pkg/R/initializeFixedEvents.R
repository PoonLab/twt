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
    e$add.event('transmission', branching.time, obj1=NA, recipientLabel, sourceLabel)
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
    eventlog$add.event('transmission',  branching.time, lineage, x$get.name(), source)
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
  
  # for each CompartmentType:
  indiv.types <- sapply(unlist(comps), function(a){a$get.type()$get.name()})
  popn.totals <- lapply(model$get.types(), function(x) {
    # 1. enumerate active compartments, including unsampled infected hosts (U) at time t=0
    U <- x$get.unsampled()
    A <- length(which(indiv.types == x$get.name()))  # number of active sampled compartments at t=0
    
    # 2. enumerate number of susceptibles (S) at time t=0
    S <- x$get.susceptible()
  
    data.frame(U=U, A=A, S=S)
  })
  
  
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

  
  current.time <- 0.0
  
  while (length(comps) > 1) {
    # draw all possible recipients that has a max sampling time less than or equal to the current.time (NAs are excluded)
    qualified.sampled.recipients <- which(time.bands <= current.time)
    if (length(qualified.sampled.recipients) < 1) {
      stop ('There must be at least one Compartment containing a Lineage with a sampling time of t=0.')
    }
    
    # retrieve transmission types from dictionary of possibleSourceTypes for each qualified sampled recipient
    possibleTransmissions <- possibleSourceTypes[ qualified.sampled.recipients ]
    
    # calc transmission rates among all source-recipient pairings of CompartmentTypes
    indiv.rates.n.quantities <- sapply(model$get.types(), function(x) {
      sourceType <- x$get.name()
      recipientTypes <- names(x$get.branching.rates())
      sapply(recipientTypes, function(y) {
        # rate = B(Us+Is)(Sr+1)
        pairRate <- x$get.branching.rate(y) * (popn.totals[[y]]$S + 1) * (popn.totals[[sourceType]]$A + popn.totals[[sourceType]]$U)  
        
        qualified.r <- which(names(possibleTransmissions) == y)
        if (length(qualified.r) == 0) {
          nPairs <- 0
        } else {
          qualified.sr <- which(sapply(qualified.r, function(z) {
            sourceType %in% possibleTransmissions[z]
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
    
    # total rate of ANY transmission event occurring is the weighted sum of these rates in the dictionary
    total.rate <- sum(indiv.rates.n.quantities)
    
    # sample waiting time & update the current time to the new upper bound in time (waiting time)
    waiting.time <- current.time + rexp(n=1, rate=total.rate)
    current.time <- waiting.time
    
    if (length(which(time.bands <= current.time)) > length(qualified.sampled.recipients)) {
      # check if the waiting time exceeds any sampling time within the sampled compartments previously not qualifying as a recipient
      # re-start the filtering to include new qualified sampled infected recipients
      next
    } else {
      # determine source and recipient by relative transmission rate sums by type
      # sample individual recipient compartment within Types, uniformly distributed
      r_name <- sample(names(qualified.sampled.recipients), 1)
      r_ind_s_popn <- which(source.popn.names == r_name)
      r_ind_comps <- which(compnames == r_name)
      
      recipient <- comps[[ r_ind_comps ]]
      r_type <- recipient$get.type()$get.name()
      
      # remove recipient from relevant lists
      comps[[ r_ind_comps ]] <- NULL
      compnames <- compnames[-r_ind_comps]
      source.popn[[ r_ind_s_popn ]] <- NULL
      source.popn.names <- source.popn.names[-r_ind_s_popn]
      time.bands <- time.bands[ -which(names(time.bands) == r_name) ]
      
      # sample individual source compartment within Types, uniformly distributed
      source <- sample(source.popn, 1)[[1]]
      s_name <- source$get.name()
      
      # if source is a us_comp, now holds a sampled lineage we care about 
      s_ind_timebands <- which(names(time.bands) == s_name)
      if (is.na(time.bands[s_ind_timebands])) {
        # update "max sampling time" of particular us_comp from NA to current.time
        time.bands[s_ind_timebands] <- current.time
        # add us_comp source to list of comps (can now be a recipient)
        comps[[length(comps)+1]] <- source
        compnames[[length(compnames)+1]] <- s_name
      }
      
      # update recipient object `source` attr and `branching.time` attr
      recipient$set.source(source)
      recipient$set.branching.time(current.time)
      
      # add transmission event to EventLogger object
      eventlog$add.event('transmission', current.time, obj1='NA', r_name, s_name)
      
      # update all counts
      popn.totals[[r_type]]$S <- popn.totals[[r_type]]$S + 1
      popn.totals[[r_type]]$A <- popn.totals[[r_type]]$A - 1
    }
  }
  
  eventlog
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

