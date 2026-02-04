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
  expect_equal(names(outer.log), c('time', 'event', 'from.comp', 'to.comp', 
                                   'from.host', 'to.host'))
  
  event <- list(time=1.0, event='transmission', from.comp='S', to.comp='I', 
                from.host='S', to.host='I')
  outer$add.event(event)
  expect_true(outer$get.nrow()==1)
})
