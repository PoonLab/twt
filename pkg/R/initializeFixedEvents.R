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
init.branching.events <- function(inputs, eventlog) {
  # @param inputs = MODEL object
  # @param eventlog = EventLogger object

  # if the user input includes a tree (host tree) then add transmission events
  comps <- inputs$get.compartments()
  lineages <- inputs$get.lineages()

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
}



# Case 3: no host tree provided, transmission events need to be generated
generate.transmission.events <- function(model, eventlog) {
  # treeswithintrees/Wiki/Simulation Pseudocode step 3 & 4; refer to issue # 12 discussion
  # simulate transmission events and fix them to the timeline of lineage sampled events
  # @param model = MODEL object
  # @param eventlog = EventLogger object
  
  sampled_comps <- model$get.extant_comps()                # init Compartments currently existing w/ Lineages at sampling.time = 0
  sampled_compnames <- model$get.names(sampled_comps)
  not_yet_sampled_comps <- model$get.non_extant_comps()    # init remaining Compartments w/ Lineages not yet sampled
  nys_compnames <- model$get.names(not_yet_sampled_comps)
  us_comps <- model$get.unsampled.hosts()
  us_compnames <- model$get.names(us_comps)
  
  # for each CompartmentType:
  indiv.types <- sapply(unlist(sampled_comps), function(a){a$get.type()$get.name()})
  active.lineages <- sapply(model$get.extant_lineages(), function(b){b$get.location()$get.type()$get.name()})
  popn.totals <- lapply(model$get.types(), function(x) {
    # 1. enumerate active compartments, including unsampled infected hosts (U) at time t=0
    U <- x$get.unsampled()
    A <- length(which(indiv.types == x$get.name()))
    
    # 2. enumerate active lineages of infected (I), pairs of active lineages within hosts at time t=0
    I <- length(which(active.lineages == x$get.name()))
    
    # 3. enumerate number of susceptibles (S) at time t=0
    S <- x$get.susceptible()
    
    data.frame(U=U, A=A, I=I, S=S)
  })
  
  
  
  # record relevant sampling times of lineages for each Compartment
  all.lineages <- model$get.lineages()
  lineage.locations <- sapply(all.lineages, function(x) {x$get.location()$get.name()})
  time.bands <- vector()         # vector of minimum sampling times for each Compartment
  max.s.times <- vector()        # vector of maximum sampling times for each Compartment
  for (x in 1:length(model$get.compartments())) {
    compName <- model$get.compartments()[[x]]$get.name()
    single.comp.lineages <- all.lineages[ which(lineage.locations == compName) ]
    single.comp.sampling.times <- sapply(single.comp.lineages, function(b) {
      b$get.sampling.time()
    })
    time.bands <- c(time.bands, min(single.comp.sampling.times))     # NOTE: only the unique ones matter
    max.s.times <- c(max.s.times, max(single.comp.sampling.times))
    names(time.bands) <- c(names(time.bands)[nzchar(x=names(time.bands))], compName)
    names(max.s.times) <- c(names(max.s.times)[nzchar(x=names(max.s.times))], compName)
  }

  
  
  # calc transmission rates among all source-recipient pairings of CompartmentTypes
  sr.pairings <- combn(sampled_compnames, 2)         # TODO: include unsampled infected hosts
  sr.pair.dict <- matrix(nrow=2, ncol=2*ncol(sr.pairings))
  sr.indiv.rates <- vector(mode='numeric', length=2*ncol(sr.pairings))
  
  for(x in 1:ncol(sr.pairings)) {
    comp1 <- sampled_comps[[which(sampled_compnames == sr.pairings[1,x])]]
    comp2 <- sampled_comps[[which(sampled_compnames == sr.pairings[2,x])]]
    comp1_type <- comp1$get.type()$get.name()
    comp2_type <- comp2$get.type()$get.name()
    
    # comp1 as source, comp2 as recipient
    sr.pair.dict[1, x*2-1] <- sr.pairings[1,x]          # row 1 = source
    sr.pair.dict[2, x*2-1] <- sr.pairings[2,x]          # row 2 = recipient
    sr.indiv.rates[x*2-1] <- model$get.types()[[comp1_type]]$get.branching.rate(comp2_type) * (1 / (popn.totals[[comp1_type]]$U + popn.totals[[comp1_type]]$A)) # weighted by number of source and recipient Types
    
    # comp2 as source, comp1 as recipient
    sr.pair.dict[1, x*2] <- sr.pairings[2,x]
    sr.pair.dict[2, x*2] <- sr.pairings[1,x]
    sr.indiv.rates[x*2] <- model$get.types()[[comp2_type]]$get.branching.rate(comp1_type) * (1 / (popn.totals[[comp2_type]]$U + popn.totals[[comp2_type]]$A))
  }
  
  # total rate of ANY transmission event occurring is the weighted sum of these rates in the dictionary
  total.rate <- sum(sr.indiv.rates)
  
  
  
  current.time <- 0.0
  
  # determine the "time band" (time interval between the last sampling event and the next sampling event across all Types)
  if (length(time.bands) == 0) {
    # the next time interval (bandwidth) is infinite
  } else {
    # remove the first element in the listed chunks of time.bands
    time.bands <- time.bands[-1]
    # set the minimum element in time.bands as the new upper bound in time
    bandwidth <- min(time.bands)
  }
  
  waiting.time <- current.time + rexp(n=1, rate=total.rate)
  # if waiting time exceeds the current time band, proceed to next time band, and update current.time
  if (waiting.time > bandwidth) {
    current.time <- bandwidth
    next
  } else {
    while (TRUE) {
      # determine source and recipient by relative transmission rate sums by type
      # sample individual source and recipient compartments within Types, uniformly distributed
      legitPairings <- sr.pair.dict[, which(sr.indiv.rates != 0)]
      sr.pair <- legitPairings[, sample(1:ncol(legitPairings), 1)]
      source <- sampled_comps[[which(sampled_compnames == sr.pair[1])]]
      recipient <- sampled_comps[[which(sampled_compnames == sr.pair[2])]]
      
      # check if the recipient has another sampling time that is earlier in time than the projected transmission time
      if (recipient$get.name() %in% names(max.s.times)) {
        if (max.s.times[recipient$get.name()] > current.time) {
          # this sr.pairing must be disqualified
          TRUE                                                  # FIXME: could go infinitely if all pairs have earlier sampling times than the current.time
        } else {
          FALSE
        }
      }
    }
    
    # update recipient object `source` attr and `branching.time` attr
    recipient$set.source(source)
    recipient$set.branching.time(current.time)
    
    # add transmission event to EventLogger object
    eventlog$add.event('transmission', current.time, obj1=NA, recipient$get.name(), source$get.name(), cumulative=T)
    
    # update all counts
    popn.totals[[recipient$get.type()$get.name()]]$S <- popn.totals[[recipient$get.type()$get.name()]]$S + 1
    popn.totals[[recipient$get.type()$get.name()]]$I <- popn.totals[[recipient$get.type()$get.name()]]$I - 1
    
    # remove all sr.pairs in dictionary that contain the recipient as a source
    
  }
  
  
}



.to.transmission.tree <- function(eventlog) {
  # require(ape, quietly=TRUE)
  # function converts the transmission events stored in an EventLogger object into a transmission tree

  t_events <- eventlog$get.events('transmission', cumulative=FALSE)
  tips <- unlist(setdiff(t_events$compartment1, t_events$compartment2))
  internals <- unlist(intersect(t_events$compartment1, t_events$compartment2))
  root <- unlist(setdiff(t_events$compartment2, t_events$compartment1))
  
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
  
  phy <- list(tip.label=unlist(tip.label), Nnode=Nnode, edge.length=as.numeric(unlist(edge.length)), edge=edge, node.label=unlist(node.label))
  attr(phy, 'class') <- 'phylo'
  attr(phy, 'order') <- 'cladewise'
  phy
}

