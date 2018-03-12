require(R6)
require(RUnit)
require(yaml)
source('pkg/R/loadInputs.R')
settings <- yaml.load_file('tests/fixtures/test.yaml')
test <- MODEL$new(settings)

# create an EventLogger object for testing
e <- EventLogger$new()
e$add.event("transmission", 4, "NA", "I_95", "I_63")
e$add.event("transmission", 6, "NA", "I_73", "I_95")
e$add.event("transmission", 3, "NA", "I_20", "I_73")
e$add.event("transmission", 2, "NA", "I_94", "I_20")
e$add.event("transmission", 1, "NA", "I_97", "I_20")

test.get.leaves.names <- function() {
  result <- get.leaves.names(e)  # return terminal nodes
  expected <- c("I_94","I_97")
  checkEquals(expected, result)
}

test.get.nonterminals <- function() {
  result <- get.nonterminals(e) # return internal nodes of the transmission tree
  expected <- c("I_95","I_73","I_20")
  checkEquals(expected, result)
}

test.get.pairs <- function() {
  result <- test$get.pairs() # extract all pairs of pathogen lineages that may coalesce
  expected <- list(`I_1:I_2,I_1:I_3`="I_1",
                   `I_2:I_2,I_2:I_3`="I_2")
  checkEquals(expected, result)
}

test.add.pair1 <- function() {
  test$add.pair('I_1:I_2','I_1:I_1',"I_2") 
  result <- test$get.pairs() # when a Lineage is moved from one compartment to another (transmission or migration) or when a Lineage is sampled
  expected <- list(`I_1:I_2,I_1:I_3`="I_1",
                   `I_2:I_2,I_2:I_3`="I_2",
                   `I_1:I_1,I_1:I_2`="I_2")
  checkEquals(expected, result)
}

test.add.pair2 <- function() {
  test <- MODEL$new(settings)
  test$add.pair('I_1:I_3','I_1:I_2',"I_2") 
  result <- test$get.pairs() # when a Lineage is moved from one compartment to another (transmission or migration) or when a Lineage is sampled
  expected <- list(`I_1:I_2,I_1:I_3`="I_2",
                   `I_2:I_2,I_2:I_3`="I_2")
  checkEquals(expected, result)
}

test.remove.pair <- function() {
  test$remove.pair('I_2:I_3','I_2:I_2') 
  result <- test$get.pairs() # when a Lineage is moved from one compartment to another (transmission or migration) or when a Lineage is sampled
  expected <- list(`I_1:I_2,I_1:I_3`="I_2")
  checkEquals(expected, result)
}

test.get.types <- function(){
  result <- test$get.types() # retrieves CompartmentType objects in a list
  expected.popn.growth.dynamics <- cbind("start"=list(0,0.75,1.5,2.25), "end"=list(0.75,1.5,2.25,"inf"), "intercept"=list(1,-15,90,30), "slope"=list(15,41,-30,-2))
  
  checkEquals('host', result[[1]]$get.name())
  checkEquals(1, result[[1]]$get.bottleneck.size())
  checkEquals(list('host'=0.1), result[[1]]$get.branching.rates())
  checkEquals(0.1, result[[1]]$get.branching.rate('host'))
  checkEquals(list(), result[[1]]$get.migration.rates())
  checkEquals(NULL, result[[1]]$get.migration.rate('host'))
  checkEquals(expected.popn.growth.dynamics, result[[1]]$get.popn.growth.dynamics())
  checkEquals(1000, result[[1]]$get.susceptible())
  checkEquals(20, result[[1]]$get.unsampled())
  
  result[[1]]$set.bottleneck.size(2)
  checkEquals(2, result[[1]]$get.bottleneck.size())
  
  result[[1]]$set.branching.rate('host',0.2)
  checkEquals(0.2, result[[1]]$get.branching.rate('host'))
  checkEquals(list('host'=0.2), result[[1]]$get.branching.rates())
  
  result[[1]]$set.migration.rate('host',0.01)
  checkEquals(0.01, result[[1]]$get.migration.rate('host'))
  checkEquals(list('host'=0.01), result[[1]]$get.migration.rates())
  
  result[[1]]$set.name('cell')
  checkEquals('cell', result[[1]]$get.name())
  
  result[[1]]$set.susceptible(1010)
  checkEquals(1010, result[[1]]$get.susceptible())
  
  result[[1]]$set.unsampled(25)
  checkEquals(25, result[[1]]$get.unsampled())
}

test.get.unsampled.hosts <- function(){
  test <- MODEL$new(settings)
  result <- test$get.unsampled.hosts() # function creates "blank" Compartment objects for Unsampled Hosts (US), stored in lists for each section within a CompartmentType object
  checkEquals(20, length(result))
  
  result.names <-sapply(result, function(x){
    x$get.name()
  })
  expected.names <- c("US_host_1", "US_host_2", "US_host_3", "US_host_4", "US_host_5",
                      "US_host_6", "US_host_7", "US_host_8", "US_host_9", "US_host_10",
                      "US_host_11", "US_host_12", "US_host_13", "US_host_14", "US_host_15",
                      "US_host_16", "US_host_17", "US_host_18", "US_host_19", "US_host_20")
  checkEquals(expected.names, result.names)
  
  result.types <- sapply(result,function(x){
    x$get.type()$get.name()
  })
  expected.types <- rep('host', 20)
  checkEquals(expected.types, result.types)
}

test.get.compartments <- function(){
  result <- test$get.compartments() # retrieves compartments in a list
  checkEquals(2, length(result))
  result.names <-sapply(result, function(x){
    x$get.name()
  })
  expected.names <- c("I_1","I_2")
  checkEquals(expected.names, result.names)
  
  result.types <- sapply(result,function(x){
    x$get.type()$get.name()
  })
  expected.types <- rep('host', 2)
  checkEquals(expected.types, result.types)
}

test.get.lineages <- function(){
  result <- test$get.lineages() # retrieves lineages in a list
  checkEquals(6, length(result))
  
  result.hosts <- sapply(result,function(x){
    x$get.location()$get.name()
  })
  checkEquals(3, length(which(result.hosts=='I_1')))
  checkEquals(3, length(which(result.hosts=='I_2')))
  
  result.times <- sapply(result,function(x){
    x$get.sampling.time()
  })
  expected.times <- c(0.2,0,0,0.2,0,0)
  checkEquals(expected.times, result.times)
  
  result[[1]]$set.location(test$get.compartments(),'I_2')
  new.location <- result[[1]]$get.location()$get.name()
  checkEquals('I_2', new.location)
}

test.get.extant_lineages <- function(){
  result <- test$get.extant_lineages() # retrieves list of Lineages with sampling.time t=0
  checkEquals(4, length(result))
  result.hosts <- sapply(result,function(x){
    x$get.location()$get.name()
  })
  checkEquals(2, length(which(result.hosts=='I_1')))
  checkEquals(2, length(which(result.hosts=='I_2')))
  result.times <- sapply(result,function(x){
    x$get.sampling.time()
  })
  expected.times <- rep(0, 4)
  checkEquals(expected.times, result.times)
}

test.get.extant_comps <- function(){
  result <- test$get.extant_comps() # retrieves list of Compartments containing Lineages with sampling.time t=0
  checkEquals(2, length(result))
  checkEquals('I_1', result[[1]]$get.name())
  checkEquals('I_2', result[[2]]$get.name())
}

test.get.non_extant_comps <- function(){
  result <- test$get.non_extant_comps() # retrieves list of Compartments containing Lineages with sampling.time t!=0
  checkEquals(NULL, result)
}

