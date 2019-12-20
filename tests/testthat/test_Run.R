require(twt)
#setwd('~/git/twt/tests/testthat')
#source('pkg/R/loadInputs.R')


# a system of 2 compartments, 3 sampled lineages each
settings <- yaml.load_file('test.yaml')
m <- Model$new(settings)
r <- Run$new(m)


test_that("Run is deep copy of Model compartments", {
  m.comps <- m$get.compartments()
  r.comps <- r$get.compartments()
  
  m.c1 <- m.comps[[1]]
  m.c1$set.branching.time(1.0)
  r.c1 <- r.comps[[1]]
  expect_null(r.c1$get.branching.time())
  
  r.c1$set.branching.time(2.0)
  expect_equal(r.c1$get.branching.time(), 2.0)
  expect_equal(m.c1$get.branching.time(), 1.0)
})


test_that("Run is deep copy of Model lineages", {
  m.lines <- m$get.lineages()
  r.lines <- r$get.lineages()
  expect_equal(length(m.lines), length(r.lines))
  expect_equal(sapply(m.lines, function(x) x$get.sampling.time()),
               sapply(r.lines, function(x) x$get.sampling.time()),
               check.names=FALSE)
  
  # move one Lineage from I_1 to I_2
  r.l1 <- r.lines[[1]]
  r.l1$set.location.by.name(r$get.compartments(), "I_2")
  expect_equal(r.l1$get.location()$get.name(), "I_2")
  
  # this should be unchanged
  m.l1 <- m.lines[[1]]
  expect_equal(m.l1$get.location()$get.name(), "I_1")
})



test_that("Run extracts pairs of lineages", {
  skip("Deprecated")
  result <- r$get.pairs()
  
  # there are only one pair per host because the third lineage is not yet extant
  # at time 0 (sampling times = 0.2)
  expected <- list(
    'I_1__I_2,I_1__I_3'='I_1',
    'I_2__I_2,I_2__I_3'='I_2'
  )
  
  expect_equal(result, expected)
})

