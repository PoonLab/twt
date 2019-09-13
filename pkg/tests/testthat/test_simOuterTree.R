require(twt)
#setwd('~/git/treeswithintrees')
#source('pkg/R/simOuterTree.R')

file1 <- yaml.load_file('example2.yaml')
file2 <- yaml.load_file('example3.yaml')

test_that("build eventlog from tree string", {
  tree <- "((A:0.1,B:0.2)A:0.3,C:0.4)A:0.5;"
  e <- eventlog.from.tree(tree)
  result <- e$get.all.events()
  
  expected <- data.frame(
    event.type=rep('transmission', 2),
    time=c(0.2, 0.5),
    lineage1=rep(NA, 2),
    lineage2=rep(NA, 2),
    compartment1=c('B', 'C'),  # recipient (reverse time)
    compartment2=c('A', 'A'),  # source
    stringsAsFactors = FALSE
  )
  
  expect_equal(result, expected)
})


test_that("build eventlog from YAML", {
  settings <- yaml.load_file('example1.yaml')
  m <- MODEL$new(settings)
  e <- EventLogger$new()
  init.branching.events(m, e)
  result <- e$get.all.events()
  
  expected <- data.frame(
    event.type=rep('transmission', 2),
    time=c(3.0, 1.0),
    lineage1=c('B_1__2_1', 'C_1__3_1'),
    lineage2=rep(NA, 2),
    compartment1=c('B_1', 'C_1'),  # recipient (reverse time)
    compartment2=c('A_1', 'B_1'),  # source
    stringsAsFactors = FALSE
  )
  
  expect_equal(result, expected)
})


test_that("simulate outer tree", {
  # SI model with 10 sampled hosts
  settings <- yaml.load_file('example2.yaml')
  m <- MODEL$new(settings)
  e <- EventLogger$new()
  sim.outer.tree(m, e)
  result <- e$get.all.events()
  
  # event log should contain 9 transmission events
  expect_equal(nrow(result), 9)
})

