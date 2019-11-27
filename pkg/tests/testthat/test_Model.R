context("load inputs")

library(twt)
#setwd('~/git/treeswithintrees')
#source('pkg/R/loadInputs.R')


# a system of 2 compartments, 3 sampled lineages each
settings <- yaml.load_file('test.yaml')
m <- Model$new(settings)


test_that("Model loads initial conditions", {
  result <- m$get.initial.conds()
  
  expected <- c('originTime', 'size', 'indexType')
  expect_equal(names(result), expected)
  
  expect_equal(result$originTime, 10)
  expect_equal(result$indexType, 'host')
  
  expected <- list(host=1000)
  expect_equal(result$size, expected)
})


test_that("Model loads types", {
  result <- m$get.types()
  
  expect_equal(class(result), 'list')
  expect_equal(length(result), 1)
  expect_equal(names(result), c('host'))
  
  # TODO: run tests for each member of <result>
  expect_is(result[[1]], "CompartmentType")
  expect_true(is.expression(parse(text=result[[1]]$get.wait.time.distr())))
})


test_that("compartments have names", {
  result <- m$get.names(m$get.compartments())
  expected <- c("I_1", "I_2")
  expect_equal(result, expected)
})


test_that("lineages have names", {
  result <- m$get.names(m$get.lineages())
  # lineage name is generated from <location>__<label>_<obj>
  expected <- c("I_1__I_1", "I_1__I_2", "I_1__I_3", "I_2__I_1", "I_2__I_2", "I_2__I_3")
  expect_equal(result, expected)
})

test_that("Model assigns sampling times", {
  result <- m$get.fixed.samplings()
  expected <- rep(c(0.2, 0, 0), times=2)
  expect_equal(result, expected, check.names=FALSE)
})

test_that("Model parse population growth dynamics", {
  result <- m$get.types()[[1]]$get.popn.growth.dynamics()
  expect_equal(class(result), "matrix")
  expect_equal(nrow(result), 4)
  expect_equal(as.vector(result[,1]), c(0, 0.75, 1.5, 2.25))  # startTime
  expect_equal(as.vector(result[,2]), c(1, 50, 100, 150))  # startPopn
  expect_equal(as.vector(result[,3]), c(0.75, 1.5, 2.25, NA))  # endTime
  expect_equal(as.vector(result[,4]), c(45, 85, 125, 150))  # endPopn
  expect_equal(as.vector(result[,5]), c(-44/0.75, -35/0.75, -25/0.75, 0))  # slope
})



