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
