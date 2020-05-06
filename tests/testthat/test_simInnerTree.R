require(twt)
#setwd('~/git/twt/pkg/tests/testthat')

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
  # FIXME: actually an unsampled compartment may be active
  #        but it will not start with Lineages
  #expect_true(result > 0)
  
  # Lineages should all be extant because we skipped all other events
  #expect_equal(length(comp$get.lineages()), 
  #             length(run$get.extant.lineages(e$time, comp)))
  
  .resolve.transition(run, e)
  expect_true(comp$get.type()$get.name() == ancestral.type$get.name())
  
})


test_that("resolve_migration", {
  #set.seed(1)  # most recent event is a migration
  run <- sim.outer.tree(structSI)
  
  eventlog <- run$get.eventlog()$get.all.events()
  e <- eventlog[1,]
  expect_equal('migration', e$event.type)
  
  comps <- c(run$get.compartments(), run$get.unsampled.hosts())
  recipient <- comps[[e$compartment1]]
  
  # FIXME: this is not the case for unsampled hosts
  #        note the migration events have not been resolved
  #result <- length(recipient$get.lineages())
  #expect_equal(5, result)
  
  # effective population size = 10
  expect_equal(10, recipient$get.type()$get.effective.size())
  expect_null(recipient$get.type()$get.popn.growth.dynamics())
  expect_equal(1, recipient$get.type()$get.bottleneck.size())
  
  # if we set number of extant lineages to 10, migration should 
  # be guaranteed
  for (i in 1:(10-length(recipient$get.lineages()))) {
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
  # FIXME: this is not always true
  source <- comps[[e$compartment2]]
  #expect_true(source$is.unsampled())
  result <- length(source$get.lineages())
  #expect_equal(0, result)
  
  before <- length(source$get.lineages())
  .resolve.migration(run, e)
  
  after <- length(source$get.lineages())
  expect_equal(1, after-before)  # bottleneck is 1
  
})


test_that("resolve coalescent", {
  run <- load.outer.tree(model)
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
  run2 <- load.outer.tree(model)
  comp <- run2$get.compartments()[[1]]
  run2.lineages <- comp$get.lineages()
  expect_equal(3, length(run2.lineages))
})


test_that("resolve bottleneck", {
  run <- load.outer.tree(model)
  comp <- run$get.compartments()[[3]]
  expect_equal(1, comp$get.type()$get.bottleneck.size())
  
  result <- .resolve.bottleneck(run, comp)
  expect_equal(1, length(result))
  
  # surviving lineage should still be in original Compartment
  expect_equal(result[[1]]$get.location()$get.name(), comp$get.name())
})


test_that("resolve transmission", {
  run <- load.outer.tree(model)
  
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


test_that("random exponential deviate under linear decay", {
  # SI model, 10 compartments with 3 Lineages each sampled at time 0
  settings <- yaml.load_file("example2.yaml")
  model <- Model$new(settings)
  
  set.seed(1)
  run <- sim.outer.tree(model)
  comps <- c(run$get.compartments(), run$get.unsampled.hosts())
  
  t.events <- run$get.eventlog()$get.events('transmission')
  e <- t.events[which.min(t.events$time), ]
  
  comp <- comps[[e$compartment1]]
  expect_false(comp$is.unsampled())
  expect_equal(3, length(run$get.extant.lineages(0, comp)))
  
  ctype <- comp$get.type()
  pieces <- as.data.frame(ctype$get.popn.growth.dynamics())
  piece <- pieces[1,]  # beta = -3.8, N_0 = 20, k=3
  
  k <- run$get.num.extant(0, comp$get.name())
  result <- replicate(1e3, .rexp.coal(k, comp, 0))
  expected <- -(piece$endPopn + piece$slope * (piece$endTime-comp$get.branching.time())) /
    piece$slope
  expect_gt(0.1, expected - max(result))
  
  # expected mean from integration of pdf*t for t from 0 to Ne/\beta, given \beta<0
  func <- function(k, n0, b) { 
    k2 <- choose(k, 2)
    -(n0^(k2/b+1) * k2 * exp(-k2*log(n0)/b) / (b*k2 - k2^2))
  }
  delta.t <- piece$endTime - e$time
  expected <- func(k=3, n0=delta.t*piece$slope + piece$endPopn, b = piece$slope)
  expect_gt(0.1, expected-mean(result))
})


test_that("simulate inner tree", {
  # user-specified transmission chain (A->B->C)
  run <- load.outer.tree(model)
  tree <- sim.inner.tree(run)
  result <- as.phylo(tree)
  expect_true(is.rooted(result))
})

