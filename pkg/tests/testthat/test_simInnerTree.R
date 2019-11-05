require(twt)

# load test fixtures
path <- system.file('extdata', 'chain.yaml', package='twt')
settings <- yaml.load_file(path)
model <- Model$new(settings)

#  A-----V-----------3  lineages
#        |
#  B     +-------V---3
#                |
#  C             +---3
# ______________________
#        3   2   1   0  time

test_that("generate coalescent") {
  run <- init.branching.events(model)
  
}

test_that("generate bottleneck") {
  run <- init.branching.events(model)
  comp <- run$get.compartments()[[3]]
  
  result <- generate.bottleneck(run, comp)
}

test_that("update transmission", {
  run <- init.branching.events(model)
  inf <- c(run$get.compartments(), run$get.unsampled.hosts())
  eventlog <- run$get.eventlog()
  
  transm.events <- eventlog$get.events('transmission')
  event <- transm.events[which.min(transm.events$time),]
  
  update.transmission(run, eventlog, inf, event)
  
  expected <- data.frame(
    event.type=c('transmission', 'transmission', 'bottleneck'),
    time=c(),
    lineage1=c(),
    lineage2=c(),
    compartment1=c(),
    compartment2=c()
    )
})

