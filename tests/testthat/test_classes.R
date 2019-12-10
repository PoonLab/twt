library(twt)

test_that("Compartment cloning", {
  comp <- Compartment$new()
  line <- Lineage$new(location=comp)
  comp$add.lineage(line)
  
  result <- comp$get.lineages()
  expect_true(identical(result[[1]], line))
  
  comp2 <- comp$copy(deep=TRUE)
  result <- comp2$get.lineages()
  
  expect_false(identical(result, line))
  
})

test_that("set sampling time", {
  comp <- Compartment$new()
  expect_equal(NA, comp$get.sampling.time())
  
  line <- Lineage$new(location=comp, sampling.time=1)
  comp$add.lineage(line)
  expect_equal(1, comp$get.sampling.time())
  
  line <- Lineage$new(location=comp, sampling.time=3)
  comp$add.lineage(line)
  expect_equal(3, comp$get.sampling.time())
  
  line <- Lineage$new(location=comp, sampling.time=2)
  comp$add.lineage(line)
  expect_equal(3, comp$get.sampling.time())
})

