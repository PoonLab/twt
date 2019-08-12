require(twt)

#setwd('~/git/treeswithintrees')
#source('pkg/R/classes.R')

# FIXME: why is this here?
#settings <- yaml.load_file('test.yaml')
#test <- MODEL$new(settings)


# create an EventLogger object for testing
e <- EventLogger$new()
e$add.event("transmission", 4, "NA", "I_95", "I_63")
e$add.event("transmission", 6, "NA", "I_73", "I_95")
e$add.event("transmission", 3, "NA", "I_20", "I_73")
e$add.event("transmission", 2, "NA", "I_94", "I_20")
e$add.event("transmission", 1, "NA", "I_97", "I_20")


test_that('get.all.events(cumul=F) returns correct df', {
  result.noncumul <- e$get.all.events(cumulative = F)
  expected.noncumul <- data.frame(
    event.type=rep('transmission', 5),
    time=c(4, 6, 3, 2, 1),
    lineage1=rep('NA', 5),
    lineage2=rep(NA, 5),
    compartment1=c('I_95', 'I_73', 'I_20', 'I_94', 'I_97'),
    compartment2=c('I_63', 'I_95', 'I_73', 'I_20', 'I_20'),
    stringsAsFactors = FALSE
  )
  row.names(expected.noncumul) <- 1:5
  expect_equal(result.noncumul, expected.noncumul)
})


test_that('get.all.events(cumul=T) returns correct df', {
  result.noncumul <- e$get.all.events(cumulative = T)
  expected.noncumul <- data.frame(
    event.type=rep('transmission', 5),
    time=c(11, 5, 2, 0, 1),
    lineage1=rep('NA', 5),
    lineage2=rep(NA, 5),
    compartment1=c('I_95', 'I_73', 'I_20', 'I_94', 'I_97'),
    compartment2=c('I_63', 'I_95', 'I_73', 'I_20', 'I_20'),
    stringsAsFactors = FALSE
  )
  row.names(expected.noncumul) <- 1:5
  expect_equal(result.noncumul, expected.noncumul)
})


# TODO: test get_events() by adding other types of events in fixture
