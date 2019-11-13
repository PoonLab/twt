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



test_that("remove lineage pairs", {
  run <- init.branching.events(model)
  lineage <- run$get.lineages()[[1]]  # A_1__1_1

  remove.lineage.pairs(run, lineage)
  result <- run$get.pairs()  # list keyed by Lineage name tuples
  
  expected <- list(
    'A_1__1_2,A_1__1_3' = "A_1",
    'B_1__2_1,B_1__2_2' = "B_1",
    'B_1__2_1,B_1__2_3' = "B_1",
    'B_1__2_2,B_1__2_3' = "B_1",
    'C_1__3_1,C_1__3_2' = "C_1",
    'C_1__3_1,C_1__3_3' = "C_1",
    'C_1__3_2,C_1__3_3' = "C_1"
    )
  expect_equal(expected, result)
})


test_that("add lineage pairs", {
  run <- init.branching.events(model)
  comp <- run$get.compartments()[[1]]
  lineage <- Lineage$new(name="test", sampling.time=0, location=comp)
  comp$add.lineage(lineage)
  
  add.lineage.pairs(run, lineage)
  result <- run$get.pairs()
  
  expect_equal(12, length(result))
  expect_true(is.element(list('A_1__1_1,test'='A_1'), result))
  expect_true(is.element(list('A_1__1_2,test'='A_1'), result))
  expect_true(is.element(list('A_1__1_3,test'='A_1'), result))
})


test_that("generate migration between sampled Compartments", {
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
  run <- init.branching.events(model)
  compartments <- run$get.compartments()
  source <- compartments[[1]]  # A_1
  recipient <- compartments[[2]]  # B_1
  lineage <- source$get.lineages()[[1]]  # A_1__1_1
  
  
  expect_error(generate.migration(run, source, recipient, lineage, 2))
})


test_that("generate migration from unsampled Compartment", {
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
  
})
