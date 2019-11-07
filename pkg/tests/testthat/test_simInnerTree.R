require(twt)

# load test fixtures
path <- system.file('extdata', 'chain.yaml', package='twt')
settings <- yaml.load_file(path)
model <- Model$new(settings)

#  A-----V-----------3  lineages
#        |
#  B     +-------V---3
#                |
#  C             +---3
# ______________________
#        3   2   1   0  time



test_that("generate coalescent", {
  run <- init.branching.events(model)
  
  # cause the first two Lineages in Compartment 1 to coalesce
  comp <- run$get.compartments()[[1]]
  lineages <- comp$get.lineages()
  line1 <- lineages[[1]]
  line2 <- lineages[[2]]
  
  generate.coalescent(run, line1, line2, 0.5)
  
  result <- run$get.lineages()
  expect_equal(8, length(result))
  expect_equal(9, length(model$get.lineages()))  # confirm deep copy
  
  result <- run$get.eventlog()$get.events('coalescent')
  rownames(result) <- NULL  # ignore row names
  
  expected <- data.frame(
    event.type=rep('coalescent', 2),
    time=rep(0.5, 2),
    lineage1=c(line1$get.name(), line2$get.name()),
    lineage2=rep('Node1', 2),
    compartment1=rep(comp$get.name(), 2),
    compartment2=rep(NA, 2),
    stringsAsFactors = FALSE
  )
  expected$compartment2 <- as.character(expected$compartment2)
  
  expect_equal(expected, result)
  
  # is Run a deep copy of Model?
  comp <- model$get.compartments()[[1]]
  model.lineages <- comp$get.lineages()
  expect_equal(3, length(model.lineages))  # should be unchanged
    
  # if we re-initialize a Run from the same Model, these
  # lineages should be unchanged
  run2 <- init.branching.events(model)
  comp <- run2$get.compartments()[[1]]
  run2.lineages <- comp$get.lineages()
  expect_equal(3, length(run2.lineages))
})


test_that("generate bottleneck", {
  run <- init.branching.events(model)
  comp <- run$get.compartments()[[3]]
  expect_equal(1, comp$get.type()$get.bottleneck.size())
  
  result <- generate.bottleneck(run, comp)
  expect_equal(1, length(result))
})


test_that("update transmission", {
  #model <- Model$new(settings)
  
  run <- init.branching.events(model)
  
  inf <- c(run$get.compartments(), run$get.unsampled.hosts())
  eventlog <- run$get.eventlog()
  
  transm.events <- eventlog$get.events('transmission')
  transm.event <- transm.events[which.min(transm.events$time), ]
  
  # first transmission event is B->C
  update.transmission(run, inf, transm.event)
  result <- run$get.eventlog()$get.all.events()
  
  expected <- c(rep('transmission', 2), rep('bottleneck', 4))
  expect_equal(expected, result$event.type)  
  
  expected <- c(3, 1, rep(1, 4))
  expect_equal(expected, result$time)
})

