require(twt)
setwd('~/git/treeswithintrees')
source('pkg/R/simOuterTree.R')

file1 <- yaml.load_file('tests/fixtures/example2.yaml')
file2 <- yaml.load_file('tests/fixtures/example3.yaml')

# Case 1: single Type, multiple sampled infected, no unsampled infected
settings <- file1
settings$CompartmentTypes[[1]]$unsampled <- 0
model <- MODEL$new(settings)
eventlog <- EventLogger$new()
case1 <- sim.outer.tree(model, eventlog)

# Case 2: single Type, multiple sampled infected, multiple unsampled infected
settings <- file1
model <- MODEL$new(settings)
eventlog <- EventLogger$new()
case2 <- sim.outer.tree(model, eventlog)

# Case 3: multiple Types, multiple sampled infected, no unsampled infected, selective zero branching rates
settings <- file2
settings$CompartmentTypes[[1]]$unsampled <- 0
settings$CompartmentTypes[[2]]$unsampled <- 0
model <- MODEL$new(settings)
eventlog <- EventLogger$new()
case3 <- sim.outer.tree(model, eventlog)

# Case 4: multiple Types, multiple sampled infected, multiple unsampled infected, selective zero branching rates
settings <- file2
model <- MODEL$new(settings)
eventlog <- EventLogger$new()
case4 <- sim.outer.tree(model, eventlog)

