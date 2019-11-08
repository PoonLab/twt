require(twt)
#setwd('~/git/treeswithintrees')
#source('pkg/R/drawWaitingTime.R')

settings <- yaml.load_file('test.yaml')
model <- Model$new(settings)
run <- Run$new(model)


test_that('extant compartments', {
  
})