require(twt)


test_that("InnerTree constructor", {
  settings <- yaml.load_file("test_SIR.yaml")
  # unspecified in YAML to test default values
  settings$Compartments$I$coalescent.rate <- 1.0
  
  mod <- Model$new(settings)
  set.seed(127)
  dynamics <- sim.dynamics(mod)
  outer <- sim.outer.tree(dynamics)
  inner <- InnerTree$new(outer)
  
  result <- inner$get.active()$count.type()
  expect_equal(result, 0)
})


test_that("Relabel inner events", {
  events <- data.frame(
    time=c(3, 2, 1, 0, 0),
    event=c('sampling', 'sampling', 'transmission', 'coalescent', 'coalescent'),
    from.comp=c("I", "I", "S", "I", "I"),
    to.comp=c("I_samp", "I_samp", "I", NA, NA),
    from.host=c("I_1", "I_2", "I_1", "I_1", "I_1"),
    to.host=c(NA, NA, "I_2", NA, NA),
    pathogen1=c("P_1", "P_2", "P_2", "P_3", "P_3"),
    pathogen2=c(NA, NA, NA, NA, NA)
  )
  result <- .relabel.inner.events(events)
  expected <- data.frame(
    # FIXME: work in progress!
  )
})