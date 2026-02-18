require(twt)

test_that("Sample coalescent waiting times", {
  settings <- yaml.load_file("test_SIR.yaml")
  mod <- Model$new(settings)
  set.seed(1276)
  eventlog <- sim.dynamics(mod)
  outer <- sim.outer.tree(mod, eventlog)
})
