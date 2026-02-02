library(twt)

test_that("Increment compartment size", {
  comp <- "A"
  init.size <- 1
  e <- new.env()
  eval(parse(text=paste(comp, "<-", 1)), envir=e)
  .plus.one(comp, e)
  result <- eval(parse(text=comp), envir=e)
  expected <- 2
  expect_equal(result, expected)
  
  .minus.one(comp, e)
  result <- eval(parse(text=comp), envir=e)
  expected <- 1
  expect_equal(result, expected)
  rm(e)
})

test_that("Instantiate and update model rates", {
  settings <- yaml.load_file("test_SIR.yaml")
  mod <- Model$new(settings)
  cnames <- mod$get.compartments()
  e <- .init.model(mod)
  expect_true(is.environment(e))
  
  result <- eval(parse(text="beta"), envir=e)
  expected <- 0.02
  expect_equal(result, expected)
  
  # sets initial sizes and rates
  rates <- .update.rates(mod, e, reset=TRUE)
  
  result <- eval(parse(text="S"), envir=e)
  expected <- 1000
  expect_equal(result, expected)
  result <- eval(parse(text="I"), envir=e)
  expected <- 1
  expect_equal(result, expected)
  result <- eval(parse(text="R"), envir=e)
  expected <- 0
  expect_equal(result, expected)
  
  expect_true(is.list(rates))
  result <- names(rates)
  expected <- c("birth", "death", "migration", "transmission")
  expect_equal(result, expected)
  
  result <- rates[["birth"]]
  expected <- c(S=0, I=0, I_samp=0, R=0)
  expect_equal(result, expected)
  result <- rates[['death']]
  expect_equal(result, expected)
  
  result <- rates[['migration']]
  expected <- matrix(c(
    0, 0, 0, 0,
    0, 0, 0.005, 0.01,
    0, 0, 0, 0,
    0, 0, 0, 0), nrow=4, byrow=TRUE, dimnames=list(cnames, cnames))
  expect_equal(result, expected)
  
  .minus.one('S', e)  # force transmission
  .plus.one('I', e)
  new.rates <- .update.rates(mod, e)
  result <- new.rates[['migration']]
  expected['I', 'I_samp'] <- 0.01
  expected['I', 'R'] <- 0.02
  expect_equal(result, expected)
  
  result <- new.rates[['transmission']]
  # array populates columns before rows
  expected <- array(data=c(
    # to S       to I
    0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  # by S
    0, 0, 0, 0,  0.02*999*2, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  # by I
    0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  # by I_samp
    0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0  # by R
  ), dim=c(4, 4, 4), dimnames=list(cnames, cnames, cnames))
  expect_equal(result, expected)
})


test_that("Calculate counts from event log", {
  settings <- yaml.load_file("test_SIR.yaml")
  mod <- Model$new(settings)
  eventlog <- data.frame(
    time=c(0.1, 0.2, 0.3, 0.4),
    event=c('transmission', 'transmission', 'migration', 'migration'),
    from.comp=c('S', 'S', 'I', 'I'),
    to.comp=c('I', 'I', 'I_samp', 'R'),
    source=c('I', 'I', NA, NA)
  )
  result <- get.counts(eventlog, mod, chunk.size=4)
  expected <- data.frame(
    time=c(0, 0.1, 0.2, 0.3, 0.4),
    S=c(1000, 999, 998, 998, 998),
    I=c(1, 2, 3, 2, 1),
    I_samp=c(0, 0, 0, 1, 1),
    R=c(0, 0, 0, 0, 1)
  )
  class(expected) <- c('twt.counts', 'data.frame')
  expect_equal(result, expected)
})


test_that("Expected SI dynamics", {
  settings <- yaml.load_file("test_SIR.yaml")
  
  # modify settings to SI model
  set.SI <- settings
  #set.SI$Sampling$targets$I_samp <- 100  # run full sim time
  set.SI$Compartments$I$migration$R <- NULL
  set.SI$Compartments$R <- NULL
  
  SI.mod <- Model$new(settings)
  event.log <- sim.dynamics(SI.mod)
})