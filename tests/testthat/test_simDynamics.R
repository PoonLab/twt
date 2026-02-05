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
  expected <- 0.002
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
    0, 0, 0, 0,  0.002*999*2, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  # by I
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
  set.SI$Parameters$simTime <- 10.0
  set.SI$Compartments$I$migration$R <- NULL
  set.SI$Compartments$R <- NULL
  set.SI$Sampling$targets[['I_samp']] <- 3
  
  SI.mod <- Model$new(set.SI)
  event.log <- sim.dynamics(SI.mod)
  
  expect_true(is.data.frame(event.log))
  result <- names(event.log)
  expected <- c("time", "event", "from.comp", "to.comp", "source")
  expect_equal(result, expected)
  
  # simulation should have full sample
  result <- sum(event.log$event=='migration' & event.log$to.comp=='I_samp')
  expected <- 3
  expect_equal(result, expected)
  
  # accumulation of infected hosts should be roughly exponential
  counts <- get.counts(event.log, SI.mod)
  y <- log(counts$I + counts$I_samp)
  x <- counts$time
  fit <- lm(y~x)
  result <- summary(fit)$adj.r.squared
  expect_gte(result, 0.88)  # this is not very sensitive
})


test_that("More thorough test of SI dynamics", {
  skip("This test requires a minute to run, skipping for routine testing")
  settings <- yaml.load_file("test_SIR.yaml")
  
  # modify settings to SI model
  set.SI <- settings
  set.SI$Parameters$simTime <- 3.5
  set.SI$Compartments$I$migration$R <- NULL
  set.SI$Compartments$R <- NULL
  set.SI$Sampling$targets$I_samp <- 100  # not feasible, run full time
  
  SI.mod <- Model$new(set.SI)
  
  # exact solution of deterministic SI model from 
  # https://davidearn.github.io/math4mb/2018/lectures/4mbl05_2018.pdf
  solve.SI <- function(t, i0, n, beta) {
    (i0*exp(n*beta*t)) / (1+(i0/n)*(exp(n*beta*t)-1))
  }
  time.pts <- seq(0.2, 3, 0.2)
  expected <- solve.SI(t=time.pts, i0=1, n=1001, 
                       beta=set.SI$Parameters$beta)
  # TODO: it would be nice to compute the stochastic mean
  
  # average counts at fixed time points
  x <- sapply(1:100, function(i) {
    suppressMessages(suppressWarnings(
      event.log <- sim.dynamics(SI.mod, max.attempts=1)
    ))
    counts <- get.counts(event.log[event.log$time <= 3,], SI.mod)
    sapply(time.pts, function(tp) {
      (counts$I + counts$I_samp)[sum(counts$time<=tp)]
    })
  })
  result <- apply(x, 1, mean)
  
  # check RMSE - note stochastic mean tends to lag
  expect_lt(sqrt(mean((expected-result)^2)), 20)
  
  # check median absolute error
  expect_lt(median(abs(result-expected)), 5)
})
