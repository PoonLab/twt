require(twt)


test_that("InnerTree constructor", {
  settings <- yaml.load_file("test_SIR.yaml")
  # unspecified in YAML to test default values
  settings$Compartments$I$coalescent.rate <- 1.0
  
  mod <- Model$new(settings)
  set.seed(127)
  eventlog <- sim.dynamics(mod)
  outer <- sim.outer.tree(mod, eventlog)
  inner <- InnerTree$new(outer, mod)
  
  result <- inner$get.active()$count.type()
  expect_equal(result, 0)
})
