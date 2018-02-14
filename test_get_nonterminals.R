require(R6)
require(yaml)
source('loadInputs.R')

test.get.nonterminals <- function() {
  e <- EventLogger$new()
  e$add.event("transmission", 4.615755e-05, "NA", "I_95", "I_63")
  e$add.event("transmission", 6.902650e-04, "NA", "I_73", "I_95")
  e$add.event("transmission", 3.345122e-04, "NA", "I_20", "I_49")
  e$add.event("transmission", 3.345122e-04, "NA", "I_94", "I_20")
  result <- get.nonterminals(e)
  expected <- c("I_95","I_20")
  checkEquals(expected, result)
}