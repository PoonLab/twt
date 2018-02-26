require(R6)
require(yaml)
source('loadInputs.R')
settings <- yaml.load_file('test.yaml')
test <- MODEL$new(settings)

# create an EventLogger object for testing
e <- EventLogger$new()
e$add.event("transmission", 4.615755e-05, "NA", "I_95", "I_63")
e$add.event("transmission", 6.902650e-04, "NA", "I_73", "I_95")
e$add.event("transmission", 3.345122e-04, "NA", "I_20", "I_49")
e$add.event("transmission", 3.345122e-04, "NA", "I_94", "I_20")

test.get.leaves.names <- function() {
  result <- get.leaves.names(e)  # return terminal nodes
  expected <- c("I_73","I_94")
  checkEquals(expected, result)
}

test.get.nonterminals <- function() {
  result <- get.nonterminals(e) # return internal nodes of the transmission tree
  expected <- c("I_95","I_20")
  checkEquals(expected, result)
}

test.get.pairs <- function() {
  result <- test$get.pairs() # extract all pairs of pathogen lineages that may coalesce
  expected <- list(`I_1:I_1,I_1:I_2`="I_1", 
                   `I_1:I_1,I_1:I_3`="I_1",
                   `I_1:I_2,I_1:I_3`="I_1",
                   `I_2:I_1,I_2:I_2`="I_2",
                   `I_2:I_1,I_2:I_3`="I_2",
                   `I_2:I_2,I_2:I_3`="I_2")
  checkEquals(expected, result)
}

test.add.pair1 <- function() {
  test$add.pair('I_1:I_2','I_1:I_1',"I_2") 
  result <- test$get.pairs() # when a Lineage is moved from one compartment to another (transmission or migration) or when a Lineage is sampled
  expected <- list(`I_1:I_1,I_1:I_2`="I_2", 
                   `I_1:I_1,I_1:I_3`="I_1",
                   `I_1:I_2,I_1:I_3`="I_1",
                   `I_2:I_1,I_2:I_2`="I_2",
                   `I_2:I_1,I_2:I_3`="I_2",
                   `I_2:I_2,I_2:I_3`="I_2")
  checkEquals(expected, result)
}

test.add.pair2 <- function() {
  test <- MODEL$new(settings)
  test$add.pair('I_1:I_2','I_2:I_1',"I_2") 
  result <- test$get.pairs() # when a Lineage is moved from one compartment to another (transmission or migration) or when a Lineage is sampled
  expected <- list(`I_1:I_1,I_1:I_2`="I_1", 
                   `I_1:I_1,I_1:I_3`="I_1",
                   `I_1:I_2,I_1:I_3`="I_1",
                   `I_2:I_1,I_2:I_2`="I_2",
                   `I_2:I_1,I_2:I_3`="I_2",
                   `I_2:I_2,I_2:I_3`="I_2",
                   `I_1:I_2,I_2:I_1`="I_2")
  checkEquals(expected, result)
}

test.remove.pair <- function() {
  test$remove.pair('I_1:I_3','I_1:I_1') 
  result <- test$get.pairs() # when a Lineage is moved from one compartment to another (transmission or migration) or when a Lineage is sampled
  expected <- list(`I_1:I_1,I_1:I_2`="I_1", 
                   `I_1:I_2,I_1:I_3`="I_1",
                   `I_2:I_1,I_2:I_2`="I_2",
                   `I_2:I_1,I_2:I_3`="I_2",
                   `I_2:I_2,I_2:I_3`="I_2",
                   `I_1:I_2,I_2:I_1`="I_2")
  checkEquals(expected, result)
}

test.get.types <- function(){
  result <- test$get.types()
  
  checkEquals(expected, result)
}

test.get.unsampled.hosts <- function(){
  result <- test$get.unsampled.hosts()
  
  checkEquals(expected, result)
}

test.get.compartments <- function(){
  result <- test$get.compartments()
  
  checkEquals(expected, result)
}

test.get.lineages <- function(){
  result <- test$get.lineages()
  
  checkEquals(expected, result)
}

test.get.extant_lineages <- function(){
  result <- test$get.extant_lineages()
  
  checkEquals(expected, result)
}

test.get.extant_comps <- function(){
  result <- test$get.extant_comps()
  
  checkEquals(expected, result)
}

test.get.non_extant_comps <- function(){
  result <- test$get.non_extant_comps()
  
  checkEquals(expected, result)
}

