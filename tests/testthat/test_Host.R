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


test_that("Pool: activate.new.slot enforces pool.size", {
  h <- Host$new(compartment="I")
  h$init.pool(2)
  h$activate.new.slot(Pathogen$new())
  h$activate.new.slot(Pathogen$new())

  expect_error(h$activate.new.slot(Pathogen$new()), "capacity")
})

test_that("Pool: activate.new.slot requires init.pool first", {
  h <- Host$new(compartment="I")
  expect_error(h$activate.new.slot(Pathogen$new()), "initializing")
})

test_that("Pool: activate.slot updates the pathogen's slot id", {
  h <- Host$new(compartment="I")
  h$init.pool(2)
  p <- Pathogen$new()
  h$activate.slot(1, p)

  result <- p$get.slot.id()
  expected <- 1
  expect_equal(result, expected)
})

test_that("Pool: sample.other.slot(include.active=FALSE) never reuses an active lineage", {
  h <- Host$new(compartment="I")
  h$init.pool(3)
  h$activate.new.slot(Pathogen$new())
  h$activate.new.slot(Pathogen$new())

  for (i in 1:20) {
    draw <- h$sample.other.slot(include.active = FALSE)
    expect_false(draw$active)
    expect_null(draw$pathogen)
  }
})

test_that("Pool: sample.other.slot(include.active=TRUE) can reuse an active lineage", {
  h <- Host$new(compartment="I")
  h$init.pool(2)
  h$activate.new.slot(Pathogen$new())
  h$activate.new.slot(Pathogen$new())
  # both slots active, excluding one leaves exactly one candidate, which
  # is active -- deterministic, no reliance on RNG luck
  draw <- h$sample.other.slot(exclude.slot.id = 1, include.active = TRUE)
  expect_true(draw$active)
})

test_that("Pool: init.pool rejects invalid n", {
  h <- Host$new(compartment="I")
  expect_error(h$init.pool(-1), "positive scalar integer")
  expect_error(h$init.pool(0), "positive scalar integer")
  expect_error(h$init.pool(1.5), "positive scalar integer")
  expect_error(h$init.pool(c(1, 2)), "positive scalar integer")
  expect_error(h$init.pool(NA_real_), "positive scalar integer")
  expect_error(h$init.pool("2"), "positive scalar integer")
})

test_that("Pool: methods error clearly before init.pool is called", {
  h <- Host$new(compartment="I")
  expect_error(h$activate.new.slot(Pathogen$new()), "initializing the host pool")
  expect_error(h$activate.slot(1, Pathogen$new()), "initializing the host pool")
  expect_error(h$deactivate.slot(1), "initializing the host pool")
  expect_error(h$sample.other.slot(), "before initializing it")
  expect_error(h$is.slot.active(1), "initializing the host pool")
  expect_error(h$get.slot.occupant(1), "initializing the host pool")
  expect_error(h$get.active.slot.ids(), "initializing the host pool")
  expect_error(h$count.active.slots(), "initializing the host pool")
  expect_error(h$count.inactive.slots(), "initializing the host pool")
})

test_that("Pool: sample.other.slot rejects a stale/invalid exclude.slot.id", {
  h <- Host$new(compartment="I")
  h$init.pool(3)
  h$activate.new.slot(Pathogen$new())
  h$activate.new.slot(Pathogen$new())

  # slot id 99 was never activated
  expect_error(h$sample.other.slot(exclude.slot.id = 99),
               "not a currently active slot")
})
