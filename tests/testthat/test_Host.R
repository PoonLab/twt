library(twt)

test_that("Count host types in set", {
  h1 <- Host$new(compartment="I1")
  h2 <- Host$new(compartment="I1")
  h3 <- Host$new(compartment="I2")
  
  hset <- HostSet$new(hosts=list(h1, h2, h3))
  
  result <- hset$get.types()
  expected <- c("I1", "I1", "I2")
  expect_equal(expected, result)
  
  result <- hset$count.type()
  expected <- 3
  expect_equal(expected, result)
  
  result <- hset$count.type("I1")
  expected <- 2
  expect_equal(expected, result)
})

test_that("Add/remove hosts from set", {
  h1 <- Host$new(compartment="I")
  hset <- HostSet$new()
  hset$add.host(h1)
  
  result <- hset$get.names()
  expected <- c("I_1")
  expect_equal(expected, result)
  
  h2 <- Host$new(compartment="I", unsampled=TRUE)
  hset$add.host(h2)
  
  result <- hset$get.names()
  expected <- c("I_1", "US_I_2")
  expect_equal(result, expected)
  
  result <- hset$remove.host.by.idx(1)
  expect_equal(result, h1)
  
  result <- hset$get.names()
  expected <- c("US_I_2")
  expect_equal(result, expected)
  
  hset$add.host(h1)
  result <- hset$get.names()
  # previous name should NOT get overwritten with new index
  expected <- c("US_I_2", "I_1")
  expect_equal(result, expected)
})

test_that("Sample hosts from set", {
  h1 <- Host$new(compartment="I1")
  h2 <- Host$new(compartment="I1")
  h3 <- Host$new(compartment="I2")
  
  hset <- HostSet$new(hosts=list(h1, h2, h3))
  
  # sampling with replacement
  sample <- sapply(1:1000, function(i) { 
    hset$sample.host()$get.compartment() })
  result <- sum(sample=="I1") / length(sample)
  expected <- 0.667
  expect_equal(result, expected, tolerance=0.065)
})


test_that("Superinfection of host", {
  host <- Host$new(name="recipient", compartment="I")
  source <- Host$new(name="source", compartment="I")
  host$set.source(source)
  host$set.transmission.time(1.0)
  
  result <- host$get.source()
  expected <- source
  expect_equal(result, expected)
  
  result <- host$get.transmission.time()
  expected <- 1.0
  expect_equal(result, expected)
  
  source2 <- Host$new(name="source2", compartment="I")
  host$set.source(source2)
  host$set.transmission.time(2.0)
  
  result <- host$get.source()
  expected <- c(source, source2)
  expect_equal(result, expected)
  
  result <- host$get.transmission.time()
  expected <- c(1.0, 2.0)
  expect_equal(result, expected)
})

