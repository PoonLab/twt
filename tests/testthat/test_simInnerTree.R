require(twt)

# generate test fixtures
settings <- yaml.load_file("test_SIR.yaml")
settings$Compartments$I$coalescent.rate <- 1.0
mod <- Model$new(settings)
set.seed(276)
dynamics <- sim.dynamics(mod)
outer <- sim.outer.tree(dynamics)


test_that("Add sampling event to inner tree", {
  inner <- InnerTree$new(outer)
  
  result <- inner$get.active()$count.type()
  expect_equal(result, 0)
  
  # activate a Host in InnerTree by sampling event
  events <- outer$get.log()
  e <- events[1,]
  .do.sampling(e, inner)
  result <- inner$get.active()$count.type()
  expect_equal(result, 1)
  
  log <- inner$get.log()
  expect_equal(nrow(log), 1)
  result <- log[1,]
  expect_equal(result$event, "sampling")
  expect_equal(result$pathogen1, "P_1")
  expect_equal(result$time, e$time)
  expect_equal(result$from.host, e$from.host)
})


test_that("Sample coalescent event", {
  inner <- InnerTree$new(outer)
  
  # manually set outer event log - note this may not agree with Host variables
  events <- data.frame(
    time=c(5, 4, 3),
    event=c('migration', 'migration', 'transmission'),
    from.comp=c("I", "I", "S"),
    src.comp=c(NA, NA, "I"),
    to.comp=c("I_samp", "I_samp", "I"),
    from.host=c("I_2", "I_1", "I_1"),
    to.host=c(NA, NA, "I_2")
  )
  
  .do.sampling(events[1,], inner)
  .do.sampling(events[2,], inner)
  
  # neither Host has more than one Pathogen
  active <- inner$get.active()
  expect_equal(active$count.type(), 2)
  result <- .rcoal(active, mod)
  expect_true(is.na(result))
  
  .do.infection(events[3,], inner)
  expect_equal(active$count.type(), 1)
  result <- active$get.names()
  expect_equal(result, "I_1")
  
  host <- active$get.host.by.name("I_1")
  expect_equal(host$count.pathogens(), 2)
  
  result <- sapply(host$get.pathogens(), function(p) p$get.name())
  expect_setequal(result, c("P_1", "P_2"))
  
  result <- .rcoal(active, mod)
  expect_equal(result$host, "I_1")
  expect_true(is.numeric(result$dt))
  expect_true(result$dt>0)
})


test_that("Coalesce pathogen lineages", {
  inner <- InnerTree$new(outer)
  
  # override outer events to get Pathogens into the same Host
  events <- data.frame(
    time=c(15, 14, 13, 12, 11),
    event=c('migration', 'migration', 'migration', 
            'transmission', 'transmission'),
    from.comp=c("I", "I", "I", "S", "S"),
    src.comp=c(NA, NA, NA, "I", "I"),
    to.comp=c("I_samp", "I_samp", "I_samp", "I", "I"),
    from.host=c("I_1", "I_2", "I_3", "I_1", "I_1"),
    to.host=c(NA, NA, NA, "I_2", "I_3")
  )
  .do.sampling(events[1,], inner)
  .do.sampling(events[2,], inner)
  .do.sampling(events[3,], inner)
  .do.infection(events[4,], inner)
  .do.infection(events[5,], inner)
  
  active <- inner$get.active()
  expect_equal(active$count.type(), 1)
  expect_equal(active$get.names(), "I_1")
  host <- active$get.host.by.name("I_1")
  expect_equal(host$count.pathogens(), 3)
  
  .do.coalescent("I_1", inner, 10.)
  expect_equal(host$count.pathogens(), 2)
  .do.coalescent("I_1", inner, 9.)
  expect_equal(host$count.pathogens(), 1)
  expect_warning(.do.coalescent("I_1", inner, 8.))
  expect_equal(host$count.pathogens(), 1)  # no coalescence
})


test_that("Transfer lineage by superinfection", {
  # make a very simple outer tree with three hosts
  mod <- Model$new(yaml.load_file("test_super.yaml"))
  set.seed(19)
  dyn <- sim.dynamics(mod)
  outer <- sim.outer.tree(dyn)
  inner <- InnerTree$new(outer)
  
  # override outer event log
  events <- data.frame(
    time=c(0.5, 0.4, 0.3),
    event=c('migration', 'transmission', 'transmission'),
    from.comp=c('I', 'I', 'S'),
    src.comp=c(NA, 'I', 'I'),
    to.comp=c('I_samp', 'I', 'I'),
    from.host=c('I_1', 'US_I_3', 'US_I_8'),
    to.host=c(NA, 'I_1', 'I_1')
  )
  
  active <- inner$get.active()
  expect_equal(active$count.type(), 0)
  .do.sampling(events[1,], inner)
  expect_equal(active$count.type(), 1)
  
  host <- active$get.host.by.name("I_1")
  expect_false(is.null(host))
  expect_equal(host$count.pathogens(), 1)
  
  .do.infection(events[2,], inner)
  expect_equal(host$count.pathogens(), 0)
  expect_false(host$get.name() %in% active$get.names())
  
  result <- .do.infection(events[3,], inner)
  expect_false(result)
  expect_equal(active$get.names(), c("US_I_3"))
})
