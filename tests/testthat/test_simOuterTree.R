require(twt)

test_that("do migration", {
  settings <- yaml.load_file("test_SIR.yaml")
  mod <- Model$new(settings)
  outer <- OuterTree$new(mod)
  expect_equal(outer$get.nrow(), 0)
  
  active <- HostSet$new(settings)
  e <- list(time=1.0, event='migration', from.comp='I', to.comp='R', source=NA,
         S=98, I=1, I_samp=0, R=1)
  
  # if no Hosts are active, then the migration is not recorded
  .do.migration(e, active, outer)
  outer.log <- outer$get.log()
  expect_equal(nrow(outer.log), 0)
  
  settings$Compartments$R$size <- 1  # change initial size
  mod <- Model$new(settings)
  outer <- OuterTree$new(mod)
  
  migrant <- Host$new(name="R_1", compartment='R')
  active$add.host(migrant)
  
  # only one member in R, so this should be recorded
  .do.migration(e, active, outer)
  outer.log <- outer$get.log()
  expect_equal(nrow(outer.log), 1)
})