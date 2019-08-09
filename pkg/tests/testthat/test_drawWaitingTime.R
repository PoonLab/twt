require(twt)
setwd('~/git/treeswithintrees')
source('pkg/R/drawWaitingTime.R')
settings <- yaml.load_file('tests/fixtures/test.yaml')
model <- MODEL$new(settings)

wait.time <- function(k, alpha, beta){
  u <- 0.5
  (1-(1-u)^(beta/choose(k,2)))*alpha/beta
}

comps <- model$get.compartments()
compnames <- model$get.names(comps)

# retrieves compartments with multiple extant lineages
extant_comps <- unique(sapply(
  unname(model$get.pairs()),
  function(x) {
    comps[[which(compnames == x)]]
  }
))

# retrieves compartment names for extant lineages
ext.lineages.compnames <- sapply(
  model$get.extant_lineages(0.0),
  function(x) {
    x$get.location()$get.name()
  }
)

# counts the number of extant lineages in one compartment
# calculated for parameter `k` in function `wait.time`
num.ext.lineages <- function(x) {
  length(which(ext.lineages.compnames==x))
}

test.wait.time <- function() {
  result <- wait.time(8, 20, 5)  
  expected <- 0.465690037
  checkEquals(expected, result)
}

test.extant_comps <- function() {
  result.names <-sapply(extant_comps, function(x){
    x$get.name()
  })
  expected.names <- c("I_1","I_2")
  checkEquals(expected.names, result.names)
}

test.ext.lineages.compnames <- function() {
  result <- ext.lineages.compnames
  expected <- c("I_1","I_1","I_2","I_2")
  checkEquals(expected, result)
}

test.num.ext.lineages <- function() {
  result <- num.ext.lineages("I_1")
  expected <- 2
  checkEquals(expected, result)
  
  result <- num.ext.lineages("I_2")
  expected <- 2
  checkEquals(expected, result)
}