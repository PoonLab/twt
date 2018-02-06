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
.to.eventlog <- function() {
  # function converts an ape::phylo tree object into transmission events stored in a NEW EventLogger object
  
}


# Case 2 : User manually input a host transmission tree into YAML format under Compartment header
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




# Case 3: 
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
    total.rate <- base.rate * (popN['U'] + popN['I']) * popN['S']  # total event rate
    
    #if (total.rate == 0) {}  #can't accept this rate
    
    # claculate change in time between source transmitting to recipient
    delta_t <- rexp(n = 1, rate = as.numeric(total.rate))
    
    
    # update recipient object `source` attr
    recipient$set.source(source)                

    # add transmission event to EventLogger object
    eventlog$add.event('transmission', delta_t, obj1=NA, recipient$get.name(), source$get.name())  # argument `lineage` is determined later at coalescence
    
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
  require(igraph, quietly=TRUE)
  
  t_events <- eventlog$get.events('transmission')
  edges <- paste(t(cbind(t_events$compartment2, t_events$compartment1)))
  graph(edges=edges)
}

