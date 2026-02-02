library(twt)

test_that("Initialize OuterTree object", {
  settings <- yaml.load_file("test_SIR.yaml")
  mod <- Model$new(settings)
  outer <- OuterTree$new(mod)
  
  result <- outer$get.targets()
  expected <- list('I_samp'=10)
  expect_equal(result, expected)
  
  outer.log <- outer$get.log()
  expect_true(is.data.frame(outer.log))
  expect_true(nrow(outer.log)==0)
  expect_equal(names(outer.log), c('time', 'event', 'host.anc', 'host.des', 
                                   'comp.anc', 'comp.des'))
  
  event <- c(time=1.0, event='transmission', host.anc='S',
             host.des='I', comp.anc='S', comp.des='I')
})
