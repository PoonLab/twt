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
    type1=rep(NA, 2),
    type2=rep(NA, 2),
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
    type1=rep(NA, 2),
    type2=rep(NA, 2),
    stringsAsFactors = FALSE
  )
  
  expect_equal(result, expected)
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




test_that("check simple SI model", {
  run <- Run$new(model.SI)
  types <- run$get.types()
  expect_equal(1, length(types))
  
  init.conds <- run$get.initial.conds()
  popn.rates <- .get.rate.matrices(types)
  
  init.samplings <- sapply(run$get.compartments(), function(x) x$get.sampling.time())
  expect_true(all(init.samplings==0))
  
  # no migration, no transition)
  expect_equal(0, as.vector(popn.rates[['transition']]))
  expect_equal(0, as.vector(popn.rates[['migration']]))
  
  # look at distribution of waiting times to first transmission
  ot <- init.conds$originTime
  result <- sapply(1:50, function(i) {
    sim <- .sample.outer.events(run)
    c(ot - sim$time[sim$event.type=='transmission'][1],
      ot - sim$time[49])
  })
  # some tolerance due to low number of replicates
  expected <- 1 / (0.01*99)  # 1.01
  expect_gt( 0.5, abs(expected - mean(result[1,])) )
  
  
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
  

test_that("assignment of transmission events", {
  run <- Run$new(model.SI)
  events <- .sample.outer.events(run)
  
  .assign.events(run, events)
  result <- run$get.eventlog()$get.all.events()
  
  # for this model, only outer events are transmission
  expect_true(all(result$event.type=='transmission'))
  
  # output should represent a complete tree relating the sampled Compartments
  sampled <- names(run$get.compartments())
  all.cases.are.tips <- all(is.element(sampled, result$compartment1))
  
  root <- unique(result$compartment2[!is.element(result$compartment2, result$compartment1)])
  expect_equal(1, length(root))
  sampled.index.case <- is.element(root, sampled)

  expect_true( all.cases.are.tips || sampled.index.case )

  
  # transmissions should define a DAG:
  
  # 1. only one parent per node - recipients must be unique
  expect_equal(nrow(result), length(unique(result$compartment1)))
  
  # 2. transmission events of parents should always precede children
  expect_true(
    all(
      sapply(1:nrow(result), function(i) {
        recipient <- result$compartment1[i]
        source <- result$compartment2[i]
        r.time <- result$time[i]
      
        if ( is.element(source, result$compartment1) ) {
          s.time <- result$time[ which(result$compartment1 == source) ]
          (s.time > r.time)  # further back in time
        } 
        else {
          TRUE  # skip comparison
        }
      })
    ))

  # 3. no cycles
  expect_true(all(
    sapply(1:nrow(result), function(i) {
      recipient <- result$compartment1[i]
      source <- result$compartment2[i]
      while(is.element(source, result$compartment1)) {
        next.row <- which(result$compartment1==source)
        source <- result$compartment2[next.row]  # next source
        if (source == recipient) {
          # cycle!
          return(FALSE)
        }
      }
      TRUE
    })
  ))
  
  # end of tests
})



test_that("assignment of outer events", {
  run <- Run$new(model.structSI)
  events <- .sample.outer.events(run)
  
  .assign.events(run, events)
  log <- run$get.eventlog()$get.all.events()
  
  # recipient should never be its own source for contact events
  contacts <- log[log$event.type != 'transition', ]
  expect_true(all(contacts$compartment1 != contacts$compartment2))
  
  # migrations are only between infected compartments (preceded by transmission)
  migrations <- log[log$event.type == 'migration', ]
  transmissions <- log[log$event.type == 'transmission', ]
  root <- unique(transmissions$compartment2[
    !is.element(transmissions$compartment2, transmissions$compartment1)
    ])
  expect_equal(1, length(root))
  expect_equal(nrow(transmissions), length(unique(transmissions$compartment1)))
  
  expect_true(all(
    sapply(1:nrow(migrations), function(i) {
      mig.time <- migrations$time[i]
      recipient <- migrations$compartment1[i]
      source <- migrations$compartment2[i]
      
      inf.time1 <- transmissions$time[which(transmissions$compartment1 == recipient)]
      if (length(inf.time1) == 0) inf.time1 <- Inf
      
      inf.time2 <- transmissions$time[which(transmissions$compartment1 == source)]
      if (length(inf.time2) == 0) inf.time2 <- Inf
      
      return( min(inf.time1, inf.time2) > mig.time )
      })
    ))
  
})



test_that("simulate outer tree", {
  # check if any of the example models raise warnings or errors
  for (f in Sys.glob('example*.yaml')) {
    model <- Model$new(yaml.load_file('example1.yaml'))  
    expect_silent(sim.outer.tree(model))
  }
  
  expect_silent(sim.outer.tree(model.chain))
  expect_silent(sim.outer.tree(model.SI))
  expect_silent(sim.outer.tree(model.structSI))
})

