# treeswithintrees/Wiki/Simulation Pseudocode step 2
# after the objects are generated from user inputs, we need to initialize the list of fixed events
init.fixed.events <- function(inputs) {
  
  # add lineage sampling events from Lineage objects
  lineages <- inputs$get.lineages()
  samplings <- lapply(lineages, function(x) {
    indiv.sampled <- c()
    sampling.times <- x$get.sampling.time()     # could either be a single sampling time or a vector of multiple sampling times
    res <- sapply(sampling.times, function(x) {
      indiv.sampled <- rbind(indiv.sampled, c(x))
    })
    res
  })
  samplings
  
  
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