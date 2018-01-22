# treeswithintrees/Wiki/Simulation Pseudocode step 2
# after the objects are generated from user inputs, we need to initialize the list of fixed events

# retrieve sampling time and populate tip labels / times in ape::phylo object (building it tips up)
# @param inputs = NestedCoalescent object
init.fixed.samplings <- function(inputs) {
  # add lineage sampling events from Lineage objects
  lineages <- unlist(inputs$get.lineages())
  tips.n.heights <- list(tip.label=character(), edge.length=numeric())
  for (x in 1:length(lineages)) {
    tip <- lineages[[x]]
    label <- tip$get.name()
    tip.height <- tip$get.sampling.time()
    # store label w/ corresponding tip height in new ape::phylo object (not casted into `phylo` yet)
    tips.n.heights$tip.label[x] <- label
    tips.n.heights$edge.length[x] <- tip.height
    tips.n.heights
  }
  tips.n.heights
}


# @param inputs = NestedCoalescent object
# @params eventlog = EventLogger object
init.fixed.transmissions <- function(inputs, eventlog) {
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
generate.transmission.events <- function(inputs) {
# for each CompartmentType
  types <- inputs$get.types()
  comps.types <- sapply(unlist(inputs$get.compartments()), function(a){a$get.type()$get.name()})
  
  init.data <- lapply(types, function(x) {
  # enumerate active compartments, including unsampled hosts (U) at time t=0
    list.U <- unlist(x$get.unsampled.popns())
    list.A <- sapply(names(list.U), function(y) {
      compY <- length(which(comps.types == y))
      y <- compY
    })
    
  # enumerate active lineages of infected (I), pairs of active lineages within hosts at time t=0
    lineage.times <- sapply(unlist(inputs$get.lineages()), function(b){b$get.sampling.time()})
    list.I <- length(which(lineage.times == 0)) 
    
  # enumerate number of susceptibles (S) at time t=0
    list.S <- unlist(x$get.susceptible.popns())
    
    enumerate <- data.frame(U=list.U, A=list.A, I=list.I, S=list.S)
    enumerate
  })
  init.data
}
