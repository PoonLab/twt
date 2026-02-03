library(twt)

test_that("Load SIR model", {
  settings <- yaml.load_file("test_SIR.yaml")
  mod <- Model$new(settings)
  
  result <- mod$get.parameters()
  expected <- list(simTime=10, beta=0.002, gamma=0.01, psi=0.005)
  expect_equal(result, expected)
  
  result <- mod$get.compartments()
  expected <- c("S", "I", "R", "I_samp")
  expect_setequal(result, expected)
  
  result <- mod$get.sampling()
  expected <- list(mode="compartment", targets=list("I_samp"=10))
  expect_equal(result, expected)
  
  result <- mod$get.init.sizes()
  expected <- c("S"=1000, "I"=1, "I_samp"=0, "R"=0)
  expect_equal(result, expected)
  
  result <- mod$get.infected()
  expected <- c(S=FALSE, I=TRUE, I_samp=TRUE, R=FALSE)
  expect_equal(result, expected)
  
  expect_true(mod$get.infected('I_samp'))
  expect_false(mod$get.infected('S'))
  
  # unspecified rates
  result <- mod$get.birth.rates()
  expected <- c(S="0", I="0", I_samp="0", R="0")
  expect_equal(result, expected)
  
  result <- mod$get.death.rates()
  expect_equal(result, expected)
  
  result <- mod$get.bottleneck.sizes()
  expect_equal(result, expected)
  
  result <- mod$get.coalescent.rates()
  expect_equal(result, expected)
  
  # rate matrices
  result <- mod$get.migration.rates()
  cnames <- c("S", "I", "I_samp", "R")
  expected <- matrix(
    c("0", "0", "0", "0",  # from S
      "0", "0", "psi*I", "gamma*I",  # from I
      "0", "0", "0", "0",  # from I_samp
      "0", "0", "0", "0"   # from R
    ), nrow=4, ncol=4, byrow=TRUE, dimnames=list(cnames, cnames))
  expect_equal(result, expected)
  
  result <- mod$get.transmission.rates()
  expected <- array("0", dim=c(4,4,4), dimnames=list(cnames, cnames, cnames))
  expected["S", "I", "I"] <- "beta*S*I"
  expect_equal(result, expected)
  
  result <- mod$get.graph()
  expected <- igraph::graph_from_literal(S-+I-+I_samp, I-+R)
  expect_true(igraph::isomorphic(result, expected))
})

test_that("Reject misspecified model", {
  settings <- yaml.load_file("test_SIR.yaml")
  
  # settings must contain Parameters field
  set1 <- settings
  set1$Parameters <- NULL  # copy on modify
  f <- function(s) { mod <- Model$new(s) }
  expect_error(f(set1), regexp="Parameters key missing")
  
  # simTime must be declared as positive numeric value
  set1 <- settings
  set1$Parameters$simTime <- -1
  expect_error(f(set1), regexp="must be positive")
  set1$Parameters$simTime <- 0
  expect_error(f(set1), regexp="must be positive")
  set1$Parameters$simTime <- "1"
  expect_error(f(set1), regexp="must be numeric")
  set1$Parameters$simTime <- NULL
  expect_error(f(set1), regexp="missing")
  
  # settings must contain Compartments field
  set1 <- settings
  set1$Compartments <- NULL
  expect_error(f(set1), regexp="Missing")
  set1$Compartments <- settings$Compartments  # restore
  set1$Compartments$R <- NULL
  expect_error(f(set1), regexp="undeclared")
  
  # settings must contain valid R expressions for rates
  set1 <- settings
  set1$Compartments[['I']][['migration']][['R']] <- "$%!@)"
  expect_error(f(set1), regexp="Invalid")
  set1$Compartments[['I']][['migration']][['R']] <- "banana*2"
  expect_error(f(set1), regexp="undeclared")
  
  # settings must contain Sampling field
  set1 <- settings
  set1$Sampling <- NULL
  expect_error(f(set1), regexp="required field")
  set1$Sampling <- settings$Sampling
  set1$Sampling$mode <- NULL
  expect_error(f(set1), regexp="required field")
  set1$Sampling$mode <- "willynilly"
  expect_error(f(set1), regexp="not recognized")
  set1$Sampling$mode <- "compartment"
  set1$Sampling$targets$I_samp <- 0
  expect_error(f(set1), regexp="must be a positive integer")
  
})
