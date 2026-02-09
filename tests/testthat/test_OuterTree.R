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

test_that("Convert OuterTree to phylo", {
  # modify settings for smaller event log
  settings <- yaml.load_file("test_SIR.yaml")
  settings$Parameters$beta <- 0.02
  settings$Parameters$psi <- 2.0
  settings$Sampling$targets$I_samp <- 3
  mod <- Model$new(settings)
  
  set.seed(13)
  event.log <- sim.dynamics(mod)
  
  outer <- sim.outer.tree(mod, event.log)
  phy <- as.phylo(outer)
  expect_equal(Ntip(phy), 3)
  expect_true(is.rooted(phy))
  expect_false(is.binary(phy))
  expect_true(has.singles(phy))
  
  phy2 <- collapse.singles(phy)
  expected <- read.tree(
    text="((I_2:0.04937079,I_1:0.07238023):0.1641129,I_4:0.2075166);")
  result <- comparePhylo(phy2, expected)
  expect_true(is.element("Both trees have the same tip labels.", 
              result$messages))
  expect_false(is.element("No split in common.", result$messages))
})