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


# treeswithintrees/Wiki/Simulation Pseudocode step 3 & 4
# simulate transmission events and fix them to the timeline of lineage sampled events
generate.transmission.events <- function(inputs, eventlog) {
  # @param inpiuts = NestedCoalescent object
  # @param eventlog = EventLogger object
  
  # for each CompartmentType: 
  types <- inputs$get.types()
  comps.types <- sapply(unlist(inputs$get.compartments()), function(a){a$get.type()$get.name()})
  
  init.data <- lapply(types, function(x) {
    # 1. enumerate active compartments, including unsampled hosts (U) at time t=0
    list.N_U <- unlist(x$get.unsampled.popns())
    list.N_A <- sapply(names(list.N_U), function(y) {
      compY <- length(which(comps.types == y))
      y <- compY
    })
    
    # 2. enumerate active lineages of infected (I), pairs of active lineages within hosts at time t=0
    lineage.times <- sapply(inputs$get.lineages(), function(b){b$get.sampling.time()})
    list.N_I <- length(which(lineage.times == 0)) 
    
    # 3. enumerate number of susceptibles (S) at time t=0
    list.N_S <- unlist(x$get.susceptible.popns())
    
    data.frame(U=list.N_U, A=list.N_A, I=list.N_I, S=list.N_S)
  })
  
  
  # total event rate (lambda) = intrinsic base rate (ie. rate of TypeA --> TypeB transmission) x N_TypeA x N_TypeB
  # include N_U and N_S for each type 
  # retrieves intrinsic base transmission rate of X --> Y transmission
  .get.lambda <- function(source, recipient) {
    base.rate <- types[[recipient]]$get.transmission.rate(source)
    popN <- init.data[[recipient]][source,]
    base.rate * (popN['U'] + popN['I']) * popN['S']  # total event rate
  }
  
  extant_comps <- inputs$get.extant_comps()    # what about US?
 
  # are we only doing 1 transmission to a recipient per host, or are we allowing for multiple branching?
  
  while (length(extant_comps) > 1) {
    # first round, need to specify the final recipient where we can trace backwards from
    if (length(extant_comps) == length(inputs$get.extant_comps())) {
      recipient_ind <- sample(1:length(extant_comps), 1)
      recipient <- extant_comps[[recipient_ind]]
      recipient$set.inf.time(0)                 # FIXME: arbitrarily assigning final recipient attr `inf.time` of 0
      extant_comps[[recipient_ind]] <- NULL
    }
    
    source_ind <- sample(1:length(extant_comps), 1)
    source <- extant_comps[[source_ind]]
    
    # calculate total event rate
    total.rate <- .get.lambda(source$get.type()$get.name(), recipient$get.type()$get.name())
    
    # sample waiting time & calculate source `inf.time` relative to recipient
    wait <- rexp(n = 1, rate = as.numeric(total.rate))
    inf.time <- recipient$get.inf.time() + wait
    
    # sample event type (isn't this always going to be a transmission for now?)
    
    
    recipient$set.source(source)                # update recipient object `source` attr
    source$set.inf.time(inf.time)               # update source object `inf.time` attr
    
    # add transmission event to EventLogger object
    # eventlog$add.event('transmission', inf.time, lineage, recipient, source)  # what is arg `lineage` here?
    
    # update all counts
    popN <- init.data[[recipient$get.type()$get.name()]][source$get.type()$get.name(),]
    popN['S'] <- popN['S'] + 1
    popN['I'] <- popN['I'] - 1       # what if the source is an Unsampled host??
    
    # change source to be the new recipient
    recipient <- source
    extant_comps[[source_ind]] <- NULL
    
    extant_comps
  }
  
}


