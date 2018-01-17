# treeswithintrees/Wiki/Simulation Pseudocode step 2
# after the objects are generated from user inputs, we need to initialize the list of fixed events

# retrieve sampling time and populate tip labels / times in ape::phylo object (building it tips up)
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


init.fixed.transmissions <- function(inputs, eventlog) {
  # if the user input includes a tree (host tree) then add transmission events
  comps <- unlist(inputs$get.compartments())
  lineages <- unlist(inputs$get.lineages())

  transmissions <- sapply(comps, function(x) {
    inf.time <- x$get.inf.time()
      
    if (is.R6(x$get.source())) {
      source <- x$get.source()$get.name()
      lineage <- lineages[[ which( sapply(lineages, function(y){which(y$get.location()$get.name() == x$get.name())}) == 1) ]]$get.name()
    } else {
      source <- x$get.source()
      lineage <- NA
    }

    # add transmission event to EventLogger object
    eventlog$add.event('transmission',  inf.time, lineage, x$get.name(), source)
  })
}
