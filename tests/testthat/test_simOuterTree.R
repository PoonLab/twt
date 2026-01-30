require(twt)

test_that("do migration", {
  settings <- yaml.load_file("test_SIR.yaml")
  mod <- Model$new(settings)
  outer <- OuterTree$new(mod)
  expect_equal(outer$get.nrow(), 0)
  
  e <- c(time=1.0, event='migration', src='I', dest='R')
  
})