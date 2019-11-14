require(twt)
setwd('~/git/twt/pkg/tests/testthat/')
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
  skip("not yet implemented")
  r <- sim.outer.tree(model.SI)
  e <- r$get.eventlog()
  result <- e$get.all.events()
  
  # TODO: transmission events should define a DAG
  
  # every recipient (compartment1) should appear once only
  expect_true(all(result$event.type=='transmission'))
  expect_equal(length(unique(result$compartment1)), nrow(result))
  
})


test_that("get rate matrices", {
  run <- Run$new(model.structSI)
  expect_equal(2, length(run$get.types()))
  result <- .get.rate.matrices(run$get.types())
  
  expected <- matrix(c(0.011, 0.012, 0.021, 0.022), byrow=T, nrow=2, 
                     dimnames=list(c('host1', 'host2'), 
                                   c('host1', 'host2')))
  expect_equal(expected, result[['transmission']])
  
  # function should zero out diagonal
  expected <- matrix(c(0.0, 0.014, 0.023, 0.0), byrow=T, nrow=2, 
                     dimnames=list(c('host1', 'host2'), 
                                   c('host1', 'host2')))
  expect_equal(expected, result[['transition']])
  
  expected <- matrix(c(0.001, 0.002, 0.003, 0.004), byrow=T, nrow=2, 
                     dimnames=list(c('host1', 'host2'), 
                                   c('host1', 'host2')))
  expect_equal(expected, result[['migration']])
  
})


test_that("get initial samplings", {
  # Compartment "A" has 3 Lineages sampled at t=0, 2 and 3
  run <- Run$new(model.chain)
  infected <- run$get.compartments()
  
  result <- .get.initial.samplings(infected)
  expected <- data.frame(
    type=rep('host', 3),
    init.sample=c(3.0, 0, 0)
    )
  row.names(expected) <- c('A_1', 'B_1', 'C_1')
  expect_equal(expected, result)
})


test_that("sample outer events", {
  run <- Run$new(model.structSI)
  types <- run$get.types()
  init.conds <- run$get.initial.conds()
  init.conds$originTime <- 20  # make simulation longer
  popn.rates <- .get.rate.matrices(types)
  init.samplings <- .get.initial.samplings(run$get.compartments())
  
  result <- .sample.outer.events(types, init.conds, popn.rates, init.samplings)
  
  # kind of a dumb test
  expect_true(all(is.element(result$event.type, 
                             c('transmission', 'transition', 'migration'))))
  
  # epidemic starts with 1 infected, 99 susceptible
  result <- sum(result$event.type == 'transmission')
  expect_equal(99, result)  # IF enough time in simulation
})


test_that("check simple SI model", {
  run <- Run$new(model.SI)
  types <- run$get.types()
  expect_equal(1, length(types))
  
  init.conds <- run$get.initial.conds()
  popn.rates <- .get.rate.matrices(types)
  init.samplings <- .get.initial.samplings(run$get.compartments())
  expect_true(all(init.samplings$init.sample==0))
  
  # no migration, no transition)
  expect_equal(0, as.vector(popn.rates[['transition']]))
  expect_equal(0, as.vector(popn.rates[['migration']]))
  
  # look at distribution of waiting times to first transmission
  ot <- init.conds$originTime
  result <- sapply(1:50, function(i) {
    sim <- .sample.outer.events(types, init.conds, popn.rates, init.samplings)
    c(ot - sim$time[sim$event.type=='transmission'][1],
      ot - sim$time[49])
  })
  # some tolerance due to low number of replicates
  expect_gt( 0.2, abs((1/(0.01*99)) - mean(result[1,])) )
  
  
  # exact solution of deterministic SI model from 
  # https://davidearn.github.io/math4mb/2018/lectures/4mbl05_2018.pdf
  solve.SI <- function(t, i0, n, beta) {
    (i0*exp(n*beta*t)) / (1+(i0/n)*(exp(n*beta*t)-1))
  }
  det.sol <- solve.SI(t=seq(0, 10, 0.2), i0=1, n=sum(unlist(init.conds$size)), 
                      beta=as.vector(popn.rates[['transmission']]))
  # time that I reaches 50
  expected <- approx(x=det.sol, y=seq(0, 10, 0.2), xout=50)$y
  
  # re-use simulations from before to save time
  expect_gt( 1.0, abs(expected - mean(result[2,], na.rm=T)) )
  
})
  






