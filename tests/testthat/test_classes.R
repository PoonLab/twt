library(twt)

test_that("Host cloning", {
  host <- Host$new()
  patho <- Pathogen$new(location=host)
  host$add.pathogen(patho)
  
  result <- host$get.pathogens()
  expect_true(identical(result[[1]], patho))
  
  host2 <- host$copy(deep=TRUE)
  result <- host2$get.pathogens()
  
  expect_false(identical(result, patho))
  
})

test_that("set sampling time", {
  host <- Host$new()
  expect_equal(NA, host$get.sampling.time())
  
  patho <- Pathogen$new(location=host, sampling.time=1)
  host$add.pathogen(patho)
  expect_equal(1, host$get.sampling.time())
  
  patho <- Pathogen$new(location=host, sampling.time=3)
  host$add.pathogen(patho)
  expect_equal(3, host$get.sampling.time())
  
  patho <- Pathogen$new(location=host, sampling.time=2)
  host$add.pathogen(patho)
  expect_equal(3, host$get.sampling.time())
})

