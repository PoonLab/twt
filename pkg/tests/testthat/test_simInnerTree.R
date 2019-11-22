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


test_that("generate coalescent", {
  skip('refactor')
  run <- init.branching.events(model)
  
  # cause the first two Lineages in Compartment 1 to coalesce
  comp <- run$get.compartments()[[1]]
  lineages <- comp$get.lineages()
  line1 <- lineages[[1]]
  line2 <- lineages[[2]]
  
  resolve.coalescent(run, line1, line2, 0.5)
  
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
    type1=rep(NA, 2),
    type2=rep(NA, 2),
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
  skip('refactor')
  run <- init.branching.events(model)
  comp <- run$get.compartments()[[3]]
  expect_equal(1, comp$get.type()$get.bottleneck.size())
  
  result <- generate.bottleneck(run, comp)
  expect_equal(1, length(result))
})


test_that("update transmission", {
  #model <- Model$new(settings)
  skip('refactor')
  
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




test_that("generate migration between sampled Compartments", {
  skip('refactor')
  
  run <- init.branching.events(model)
  compartments <- run$get.compartments()
  source <- compartments[[1]]  # A_1
  recipient <- compartments[[2]]  # B_1
  lineage <- recipient$get.lineages()[[1]]  # B_1__2_1
  
  # note, source *gains* a Lineage because we're going back in time!
  generate.migration(run, source, recipient, lineage, 2)
  
  result <- run$get.eventlog()$get.events('migration')
  expect_equal(1, nrow(result))
  
  row.names(result) <- NULL
  expected <- data.frame(
    event.type='migration', 
    time=2.0, 
    lineage1='B_1__2_1', 
    lineage2=NA, 
    compartment1='B_1', 
    compartment2='A_1',
    stringsAsFactors = FALSE)
  expect_equal(expected, result)
})


test_that("incorrect migration", {
  skip('refactor')
  
  run <- init.branching.events(model)
  compartments <- run$get.compartments()
  source <- compartments[[1]]  # A_1
  recipient <- compartments[[2]]  # B_1
  lineage <- source$get.lineages()[[1]]  # A_1__1_1
  
  
  expect_error(generate.migration(run, source, recipient, lineage, 2))
})


test_that("generate migration from unsampled Compartment", {
  skip('refactor')
  
  set.seed(1)
  run <- sim.outer.tree(model)
  
  eventlog <- run$get.eventlog()
  print(eventlog)
  
  source <- run$get.unsampled.hosts()[[1]]
  
  recipient <- run$get.compartments()[[2]]  # B_1
  migrating.lineage <- recipient$get.lineages()[[1]]  # B_1__2_1
  
  generate.migration(run, source, recipient, migrating.lineage, 2)
  
  result <- run$get.eventlog()$get.events('migration')
  expect_equal(1, nrow(result))
  
  row.names(result) <- NULL
  expected <- data.frame(
    event.type='migration', 
    time=2.0, 
    lineage1='B_1__2_1', 
    lineage2=NA, 
    compartment1='B_1', 
    compartment2='US_host_1',
    stringsAsFactors = FALSE)
  expect_equal(expected, result)
})



test_that("resolve migration", {
  skip('refactor')
})
