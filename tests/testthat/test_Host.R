library(twt)

test_that("Count host types in set", {
  h1 <- Host$new(compartment="I1")
  h2 <- Host$new(compartment="I1")
  h3 <- Host$new(compartment="I2")
  
  hset <- HostSet$new(hosts=list(h1, h2, h3))
  
  result <- hset$get.types()
  expected <- c("I1", "I1", "I2")
  expect_equal(expected, result)
  
  # Hosts added on initialization should have same names (default NA)
  result <- names(result)
  expected <- c(NA, NA, NA)
  expect_equal(expected_result)
  
  result <- hset$count.type()
  expected <- 3
  expect_equal(expected, result)
  
  result <- hset$count.type("I1")
  expected <- 2
  expect_equal(expected, result)
})


test_that("Add/remove hosts from set", {
  h1 <- Host$new(compartment="I")
  hset <- HostSet()
  hset$add.host()
})