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

path <- system.file('extdata', 'structSI.yaml', package='twt')
settings <- yaml.load_file(path)
structSI <- Model$new(settings)



test_that("draw coalescent wait times", {
  skip("refactoring coalescent")
  run <- Run$new(model)
  
})


test_that("resolve transition", {
  set.seed(1)  # should result in 5 transitions
  run <- sim.outer.tree(structSI)
  # comps should be in their most ancestral Types
  
  eventlog <- run$get.eventlog()$get.all.events()
  e <- eventlog[eventlog$event.type=='transition', ][1,]
  expect_true(e$type1 != e$type2)
  
  types <- run$get.types()
  derived.type <- types[[e$type1]]
  ancestral.type <- types[[e$type2]]
  comps <- c(run$get.compartments(), run$get.unsampled.hosts())
  comp <- comps[[e$compartment1]]
  
  # Compartments should still be in ancestral types due to .assign.events
  expect_true(comp$get.type()$get.name() == ancestral.type$get.name())
  # manually switch it for this test
  comp$set.type(derived.type)
  expect_true(comp$get.type()$get.name() == derived.type$get.name())
  
  # Compartment should be carrying Lineages
  result <- length(comp$get.lineages())
  expect_true(result > 0)
  
  # Lineages should all be extant because we skipped all other events
  expect_equal(length(comp$get.lineages()), 
               length(run$get.extant.lineages(e$time, comp)))
  
  .resolve.transition(run, e)
  expect_true(comp$get.type()$get.name() == ancestral.type$get.name())
  
})


test_that("resolve_migration", {
  set.seed(1)  # most recent event is a migration
  run <- sim.outer.tree(structSI)
  
  eventlog <- run$get.eventlog()$get.all.events()
  e <- eventlog[1,]
  expect_equal('migration', e$event.type)
  
  comps <- c(run$get.compartments(), run$get.unsampled.hosts())
  recipient <- comps[[e$compartment1]]
  
  result <- length(recipient$get.lineages())
  expect_equal(5, result)
  
  # effective population size = 10
  expect_equal(0.1, recipient$get.type()$get.coalescent.rate())
  expect_null(recipient$get.type()$get.popn.growth.dynamics())
  expect_equal(1, recipient$get.type()$get.bottleneck.size())
  
  # if we set number of extant lineages to 10, migration should 
  # be guaranteed
  for (i in 1:5) {
    l <- Lineage$new(name=paste("temp", i, sep='_'),
                     sampling.time=0, location=recipient)
    run$add.lineage(l)
    # note if we use Compartment$add.lineage() then Run object
    # will NOT be aware of these Lineage objects and the next
    # test will fail
  }

  result <- length(run$get.extant.lineages(0, recipient))
  expect_equal(10, result)
  
  result <- length(recipient$get.lineages())
  expect_equal(10, result)
  
  # source Compartment is unsampled
  source <- comps[[e$compartment2]]
  expect_true(source$is.unsampled())
  result <- length(source$get.lineages())
  expect_equal(0, result)
  
  .resolve.migration(run, e)
  
  result <- length(source$get.lineages())
  expect_equal(1, result)  # bottleneck is 1
  
})


test_that("resolve coalescent", {
  run <- init.branching.events(model)
  lineages <- run$get.lineages()
  lineage.names <- sapply(lineages, function(x) x$get.name())
  
  # cause two Lineages in Compartment 1 to coalesce
  comp <- run$get.compartments()[[1]]
  
  .resolve.coalescent(run, comp, 0.5)
  
  result <- run$get.lineages()
  expect_equal(8, length(result))
  
  expect_equal(9, length(model$get.lineages()))  # confirm deep copy
  
  result <- run$get.eventlog()$get.events('coalescent')
  expect_equal(2, nrow(result))
  expect_true(all(result$lineage2=='Node1'))
  expect_true( all(is.element(result$lineage1, lineage.names)) )
  expect_true( all(result$time == 0.5) )
  expect_true( all(result$compartment1 == comp$get.name()) )
  expect_true( all(is.na(result$compartment2)) )

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


test_that("resolve bottleneck", {
  run <- init.branching.events(model)
  comp <- run$get.compartments()[[3]]
  expect_equal(1, comp$get.type()$get.bottleneck.size())
  
  result <- .resolve.bottleneck(run, comp)
  expect_equal(1, length(result))
  
  # surviving lineage should still be in original Compartment
  expect_equal(result[[1]]$get.location()$get.name(), comp$get.name())
})


test_that("resolve transmission", {
  run <- init.branching.events(model)
  
  # retrieve most recent transmission event
  eventlog <- run$get.eventlog()
  transm.events <- eventlog$get.events('transmission')
  transm.event <- transm.events[which.min(transm.events$time), ]
  
  # first transmission event is: [C] <- [B]
  .resolve.transmission(run, transm.event)
  result <- run$get.eventlog()$get.all.events()
  
  expected <- c(rep('transmission', 2), rep('bottleneck', 4))
  expect_equal(expected, result$event.type)  
  
  expected <- c(3, 1, rep(1, 4))
  expect_equal(expected, result$time)
  
  expected <- c('B_1', 'C_1', rep('C_1', 4))
  expect_equal(expected, result$compartment1)
  
  expected <- c(TRUE, FALSE, rep(FALSE, 4))
  expect_equal(expected, is.na(result$lineage1))
})


test_that("simulate inner tree", {
  run <- init.branching.events(model)
  
  set.seed(1)
  tree <- sim.inner.tree(run)
  result <- as.phylo(tree)
})

