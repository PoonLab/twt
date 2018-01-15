# treeswithintrees/Wiki/Simulation Pseudocode step 2
# after the objects are generated from user inputs, we need to initialize the list of fixed events

# retrieve sampling time and populate tip labels / times in ape::phylo object (building it tips up)
init.fixed.samplings <- function(inputs, new.tree) {
  # add lineage sampling events from Lineage objects
  lineages <- inputs$get.lineages()
  ind <- 1
  samplings <- sapply(lineages, function(x) {
    sapply(x, function(y) {
      label <- y$get.name()
      tip.height <- y$get.sampling.time()
      # create label with corresponding tip height in new ape::phylo object, and populate edge matrix w/ initial values
      new.tree$tip.label[ind] <- label
      new.tree$edge.length[ind] <- tip.height
      new.tree$edge[ind,2] <- ind
      new.tree
    })
  })
  samplings
}


init.fixed.transmissions <- function(inputs) {
  # if the user input includes a tree (host tree) then add transmission events
  comps <- inputs$get.compartments()
  fixed.events <- c()
  transmissions <- sapply(comps, function(x) {
    if (is.R6(x$get.source())) {
      source <- x$get.source()$get.name()
    } else {
      source <- x$get.source()
    }
    inf.time <- x$get.inf.time()
    # add row in format c(inf.time, source, x)
    fixed.events <- rbind(fixed.events, c(inf.time, source, x$get.name()))
  })
  t(transmissions)
  
}