# treeswithintrees/Wiki/Simulation Pseudocode step 2
# after the objects are generated from user inputs, we need to initialize the list of fixed events


init.fixed.samplings <- function(inputs) {
  # retrieve sampling time and populate tip labels / times in ape::phylo object (building it tips up)
  # @param inputs = NestedCoalescent object

  # add lineage sampling events from Lineage objects
  lineages <- unlist(inputs$get.lineages())

  list(
    # store label w/ corresponding tip height in new ape::phylo object (not casted into `phylo` yet)
    tip.label = sapply(lineages, function(x) x$get.name()),

    # only used for calculating edge length
    tip.height = sapply(lineages, function(x) x$get.sampling.time())
  )
}


# Case 1 : User provides a host transmission tree
.to.eventlog <- function(newick) {
  # function converts an ape::phylo tree object into transmission events stored in a NEW EventLogger object
  e <- EventLogger$new()
  phy <- read.tree(text=newick)
  
  sapply(1:nrow(phy$edge), function(x) {
    s_ind <- phy$edge[x,1]
    r_ind <- phy$edge[x,2]
    
    if (is.null(phy$node.label)) {
      # generate unique name for internal node
      sourceLabel <- paste0('US_', s_ind-length(phy$tip.label))    ## FIXME: remove this
    } else {
      sourceLabel <- phy$node.label[s_ind]
    }
    
    if (r_ind <= length(phy$tip.label)) {
      recipientLabel <- phy$tip.label[r_ind]
    } else if (is.null(phy$node.label)) {
      recipientLabel <- paste0('US_', r_ind-length(phy$tip.label))   ## FIXME: remove this
    } else {
      recipientLabel <- phy$node.label[r_ind]
    }
    
    inf.time <- phy$edge.length[x]
    e$add.event('transmission', inf.time, obj1=NA, recipientLabel, sourceLabel)
  })
  
  e
}


# Case 2 : User manually input a host transmission tree into YAML format under Compartments header
init.fixed.transmissions <- function(inputs, eventlog) {
  # @param inputs = NestedCoalescent object
  # @param eventlog = EventLogger object

  # if the user input includes a tree (host tree) then add transmission events
  comps <- inputs$get.compartments()
  lineages <- inputs$get.lineages()

  transmissions <- sapply(comps, function(x) {
    inf.time <- x$get.inf.time()
      
    if (is.R6(x$get.source())) {
      source <- x$get.source()$get.name()
      xLin <- sapply(lineages, function(y){which(y$get.location()$get.name() == x$get.name())})
      lineage <- lineages[[ which(xLin == 1) ]]$get.name()
    } else {
      source <- x$get.source()
      lineage <- NA
    }

    # add transmission event to EventLogger object
    eventlog$add.event('transmission',  inf.time, lineage, x$get.name(), source)
  })
}




# Case 3: no host tree provided, transmission events need to be generated
generate.transmission.events <- function(inputs, eventlog) {
  # treeswithintrees/Wiki/Simulation Pseudocode step 3 & 4
  # simulate transmission events and fix them to the timeline of lineage sampled events
  # @param inputs = NestedCoalescent object
  # @param eventlog = EventLogger object
  
  # for each CompartmentType: 
  comps.types <- sapply(unlist(inputs$get.extant_comps()), function(a){a$get.type()$get.name()})  
  active.lineages <- sapply(inputs$get.extant_lineages(), function(b){b$get.location()$get.type()$get.name()})
  init.data <- lapply(inputs$get.types(), function(x) {
    
    # 1. enumerate active compartments, including unsampled hosts (U) at time t=0
    list.N_U <- unlist(x$get.unsampled.popns())
    list.N_A <- sapply(names(list.N_U), function(y) { length(which(comps.types == y)) })
    
    # 2. enumerate active lineages of infected (I), pairs of active lineages within hosts at time t=0
    list.N_I <- length(which(active.lineages == x$get.name()))
    
    # 3. enumerate number of susceptibles (S) at time t=0
    list.N_S <- unlist(x$get.susceptible.popns())
    
    data.frame(U=list.N_U, A=list.N_A, I=list.N_I, S=list.N_S)
  })
  
  extant_comps <- inputs$get.extant_comps()
  non_extant_comps <- inputs$get.non_extant_comps()
  us_comps <- inputs$get.unsampled.hosts()
  us_compnames <- sapply(us_comps, function(x){x$get.name()})
  # a source can transmit to multiple recipients, but cannot receive a transmission from one of its descendants
  # `private$transmission.history` attr for each Compartment ensures this does not happen  <- could be useless.. pending judgement
  
  # main loop
  while (length(c(extant_comps, non_extant_comps)) > 1) {
    # grab a recipient from list of extant
    r_ind <- sample(1:length(extant_comps), 1)
    recipient <- extant_comps[[r_ind]]
    extant_comps[[r_ind]] <- NULL            # remove chosen recipient from list of extant, as well as host_popn
    
    host_popn <- c(extant_comps, non_extant_comps, us_comps)
    s_ind <- sample(1:length(host_popn), 1)  # sample source from combined list of extant and unsampled-infected hosts
    source <- host_popn[[s_ind]]
    
    
    # calculate total event rate (lambda) = intrinsic base rate (ie. rate of TypeA --> TypeB transmission) x N_TypeA x N_TypeB
    # includes N_U and N_S for each type, and retrieves intrinsic base transmission rate of X --> Y transmission
    r_type <- recipient$get.type()$get.name()
    s_type <- source$get.type()$get.name()
    base.rate <- inputs$get.types()[[r_type]]$get.transmission.rate(s_type)
    popN <- init.data[[r_type]][s_type,]
    total.rate <- base.rate * (popN['U'] + popN['I']) * popN['S']  # total event rate   ... is popN['S'] for the recipient?
    
    #if (total.rate == 0) {}  #can't accept this rate
    
    # calculate change in time between source transmitting to recipient
    delta_t <- rexp(n = 1, rate = as.numeric(total.rate))
    
    
    # update recipient object `source` attr
    recipient$set.source(source)                

    # add transmission event to EventLogger object
    eventlog$add.event('transmission', delta_t, obj1=NA, recipient$get.name(), source$get.name(), cumulative=FALSE)  # argument `lineage` is determined later at coalescence
    
    # update all counts
    popN['S'] <- popN['S'] + 1
    popN['I'] <- popN['I'] - 1       # if source is US_host, US-- and I++ as it moves from host_popn into extant list (cancels out)
    
    if (source$get.name() %in% us_compnames) {
      extant_comps[[length(extant_comps)+1]] <- source
      us_ind <- which(us_compnames == source$get.name())
      us_compnames <- us_compnames[-us_ind]
      us_comps[[us_ind]] <- NULL          # if a US host enters the extant_comps, then it must be removed from the us_comps list to avoid duplication
    }
    
  }
  
  eventlog
}



.to.transmission.tree <- function(eventlog) {
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

