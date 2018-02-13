require(R6)
require(yaml)
source('loadInputs.R')



test.get.leaves.names <- function() {
  e <- EventLogger$new()
  e$add.event("transmission", 4.731931e-05, "NA", "I_87", "I_61")
  e$add.event("transmission", 9.077062e-05, "NA", "I_79", "I_65")
  e$add.event("transmission", 1.041083e-04, "NA", "I_42", "I_54")
  result <- get.leaves.names(e)
  expected <- c("I_87","I_79","I_42")
  checkEquals(expected, result)
}
