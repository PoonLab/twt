require(twt)

test_that("do migration", {
  settings <- yaml.load_file("test_SIR.yaml")
  mod <- Model$new(settings)
  outer <- OuterTree$new(mod)
  expect_equal(outer$get.nrow(), 0)
  active <- outer$get.active()
  
  e <- list(
    time=1.0, event='migration', from.comp='I', to.comp='R', source=NA,
    S=98, I=1, I_samp=0, R=1)
  
  # if no Hosts are active, then the migration is not recorded
  .do.migration(e, outer)
  
  outer.log <- outer$get.log()
  expect_equal(nrow(outer.log), 0)
  
  settings$Compartments$R$size <- 1  # change initial size
  mod <- Model$new(settings)
  outer <- OuterTree$new(mod)
  active <- outer$get.active()
  
  migrant <- Host$new(name="R_1", compartment='R')
  active$add.host(migrant)
  
  # only one member in R, so this should be recorded
  .do.migration(e, outer)
  
  outer.log <- outer$get.log()
  expect_equal(nrow(outer.log), 1)
  
  result <- outer.log[1,]
  expected <- data.frame(
    time=1.0, event='migration', from.comp='I', to.comp='R',
    from.host="R_1", to.host=as.character(NA)
    )
  expect_equal(result, expected)
  
})


test_that("do transmission", {
  # set up the model
  settings <- yaml.load_file("test_SIR.yaml")
  mod <- Model$new(settings)
  outer <- OuterTree$new(mod)
  active <- outer$get.active()
  
  # one transmission that is immediately sampled (note this includes counts)
  event.log <- data.frame(
    time = c(0, 1.0, 1.1), 
    event = c(NA, 'transmission', 'migration'), 
    from.comp = c(NA, 'S', 'I'), 
    to.comp = c(NA, 'I', 'I_samp'), 
    source = c(NA, 'I', NA), 
    S = c(100, 99, 98), 
    I = c(1, 2, 1), 
    I_samp = c(0, 0, 1), 
    R = c(0, 0, 0)
    )
  
  expect_equal(outer$nsamples(), 0)
  result <- sum(unlist(outer$get.targets()))
  expect_equal(result, 10)
  expect_equal(active$count.type(), 0)
  
  # this should activate a Host
  .do.migration(event.log[3,], outer)
  
  expect_equal(outer$nsamples(), 1)
  result <- outer$get.sampled()$count.type("I")
  expect_equal(result, 1)
  
  # sampling activates a host
  expect_equal(active$count.type(), 1)
  expect_equal(active$get.types(), c('I'))
  expect_equal(active$get.names(), c('I_1'))
  host <- active$sample.host()
  expect_equal(host$get.name(), 'I_1')
  expect_equal(host$get.sampling.time(), 1.1)
  expect_equal(host$get.transmission.time(), numeric())
  
  # at this point, there are two individuals in I.  One of them is the 
  # source and one is the recipient (I_1).
  .do.transmission(event.log[2,], event.log[1,], outer)
  
  expect_equal(outer$nsamples(), 1)  # should be unchanged
  expect_equal(active$count.type(), 1)
  
  retired <- outer$get.retired()
  
  result <- active$get.names()
  if (result == c("US_I_2")) {
    # active host was recipient and got retired
    expect_equal(retired$get.names(), "I_1")
    result <- outer$get.log()
    expected <- data.frame(
      time=c(1.1, 1.0),
      event=c('migration', 'transmission'),
      from.comp=c('I', 'S'),
      to.comp=c('I_samp', 'I'),
      from.host=c('I_1', 'US_I_2'),
      to.host=c(as.character(NA), 'I_1')
    )
    expect_equal(result, expected)
    
  } else if (result == c("I_1")) {
    # active host was not the recipient
    expected <- data.frame(
      time=c(1.1),
      event=c('migration'),
      from.comp=c('I'),
      to.comp=c('I_samp'),
      from.host=c('I_1'),
      to.host=c(as.character(NA))
    )
    result <- outer$get.log()
    expect_equal(result, expected)
    expect_equal(retired$count.type(), 0)
    
  } else {
    fail("Unexpected membership in active HostSet: ", result)
  }
  
})


test_that("full outer tree simulation", {
  settings <- yaml.load_file("test_SIR.yaml")
  settings$Sampling$targets$I_samp <- 1
  mod <- Model$new(settings)
  
  # one transmission that is immediately sampled
  event.log <- data.frame(
    time = c(1.0, 1.1), 
    event = c('transmission', 'migration'), 
    from.comp = c('S', 'I'), 
    to.comp = c('I', 'I_samp'), 
    source = c('I', NA)
  )
  
  result <- get.counts(event.log, mod)
  expected <- data.frame(
    time=c(0, 1.0, 1.1),
    S=c(1000, 999, 999),
    I=c(1, 2, 1),
    I_samp=c(0, 0, 1),
    R=c(0, 0, 0)
  )
  class(expected) <- c('twt.counts', 'data.frame')
  expect_equal(result, expected)
  
  outer <- sim.outer.tree(mod, event.log)
  expect_equal(class(outer), c("OuterTree", "R6"))
  result <- outer$get.log()
  
  # half the time transmission is from US_I_2 to I_1
  if (nrow(result) == 2) {
    expected <- data.frame(
      time=c(1.1, 1.0), 
      event=c('migration', 'transmission'), 
      from.comp=c('I', 'S'), 
      to.comp=c('I_samp', 'I'), 
      from.host=c('I_1', 'US_I_2'), 
      to.host=c(as.character(NA), 'I_1')
    )
    expect_equal(result, expected)
    
  } else {
    # other times transmission is from I_1 to unrecorded host
    expected <- data.frame(
      time=1.1, event='migration', from.comp='I', to.comp='I_samp', 
      from.host='I_1', to.host=as.character(NA)
    )
    expect_equal(result, expected)    
  }
})
