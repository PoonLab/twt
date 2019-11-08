require(twt)
#setwd('~/git/twt/pkg/tests/testthat/')
#source('pkg/R/simOuterTree.R')

# chain
model.chain <- Model$new(yaml.load_file('example1.yaml'))

# SI model with 10 sampled hosts
model.SI <- Model$new(yaml.load_file('example2.yaml'))

# load test fixtures
path <- system.file('extdata', 'structSI.yaml', package='twt')
settings <- yaml.load_file(path)
model.structSI <- Model$new(settings)


test_that("build eventlog from tree string", {
  tree <- "((A:0.1,B:0.2)A:0.3,C:0.4)A:0.5;"
  e <- eventlog.from.tree(tree)
  result <- e$get.all.events()
  
  expected <- data.frame(
    event.type=rep('transmission', 2),
    time=c(0.2, 0.5),
    lineage1=rep(NA, 2),
    lineage2=rep(NA, 2),
    compartment1=c('B', 'C'),  # recipient (reverse time)
    compartment2=c('A', 'A'),  # source
    stringsAsFactors = FALSE
  )
  
  expect_equal(result, expected)
})


test_that("build eventlog from YAML", {
  run <- init.branching.events(model.chain)
  result <- run$get.eventlog()$get.all.events()
  
  expected <- data.frame(
    event.type=rep('transmission', 2),
    time=c(3.0, 1.0),
    lineage1=rep(NA, 2),
    lineage2=rep(NA, 2),
    compartment1=c('B_1', 'C_1'),  # recipient (reverse time)
    compartment2=c('A_1', 'B_1'),  # source
    stringsAsFactors = FALSE
  )
  
  expect_equal(result, expected)
})


test_that("simulate outer tree", {
  r <- sim.outer.tree(model.SI)
  e <- r$get.eventlog()
  result <- e$get.all.events()
  
  # TODO: transmission events should define a DAG
  
  # every recipient (compartment1) should appear once only
  expect_true(all(result$event.type=='transmission'))
  expect_equal(length(unique(result$compartment1)), nrow(result))
  
})


test_that("calculate population rates", {
  run <- Run$new(model.structSI)
  expect_equal(2, length(run$get.types()))
  result <- .calc.popn.rates(run$get.types())
  
  expected <- matrix(c(0.01, 0.01, 0.02, 0.02), byrow=T, nrow=2, 
                     dimnames=list(c('host1', 'host2'), 
                                   c('host1', 'host2')))
  expect_equal(expected, result)
})


test_that("store initial samplings", {
  # Compartment "A" has 3 Lineages sampled at t=0, 2 and 3
  run <- Run$new(model.chain)
  infected <- run$get.compartments()
  lineages <- run$get.lineages()
  popn.rates <- .calc.popn.rates(run$get.types())
  
  result <- .store.initial.samplings(infected, lineages, popn.rates)
  expected <- c("A_1"=3.0, "B_1"=0, "C_1"=0)
  expect_equal(expected, result)
})


test_that("calculate transmission events", {
  
  model <- Model$new(settings)
  run <- Run$new(model)
  
  #t.events <- .calc.transmission.events(
  #  run$get.initial.conds(), 
  #  )
  
})
