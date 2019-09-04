require(twt)


# create an EventLogger object for testing
# type, time, line1, line2, comp1, comp2
e <- EventLogger$new()
e$add.event("transmission", time=1, line1="NA", comp1="I_97", comp2="I_20")
e$add.event("transmission", time=2, line1="NA", comp1="I_94", comp2="I_20")
e$add.event("transmission", time=3, line1="NA", comp1="I_20", comp2="I_73")
e$add.event("transmission", time=4, line1="NA", comp1="I_95", comp2="I_63")
e$add.event("transmission", time=6, line1="NA", comp1="I_73", comp2="I_95")



test_that('get.all.events returns correct df', {
  result <- e$get.all.events()
  
  expected <- data.frame(
    event.type=rep('transmission', 5),
    time=c(1, 2, 3, 4, 6),
    lineage1=rep('NA', 5),
    lineage2=rep(NA, 5),
    compartment1=c('I_97', 'I_94', 'I_20', 'I_95', 'I_73'),
    compartment2=c('I_20', 'I_20', 'I_73', 'I_63', 'I_95'),
    stringsAsFactors = FALSE
  )
  
  row.names(expected) <- 1:5
  expect_equal(result, expected)
})


# TODO: test get_events() by adding other types of events in fixture
