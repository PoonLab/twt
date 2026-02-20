require(twt)

# generate test fixtures
settings <- yaml.load_file("test_SIR.yaml")
settings$Compartments$I$coalescent.rate <- 1.0
mod <- Model$new(settings)
set.seed(276)
eventlog <- sim.dynamics(mod)
outer <- sim.outer.tree(mod, eventlog)


test_that("Add sampling event to inner tree", {
  inner <- InnerTree$new(outer, mod)
  
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
  inner <- InnerTree$new(outer, mod)
  
  # manually set outer event log - note this may not agree with Host variables
  events <- data.frame(
    time=c(5, 4, 3),
    event=c('migration', 'migration', 'transmission'),
    from.comp=c("I", "I", "S"),
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