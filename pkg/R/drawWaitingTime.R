require(twt)
setwd('~/git/treeswithintrees')
settings <- yaml.load_file('tests/fixtures/test.yaml')
input <- MODEL$new(settings)

comps <- input$get.compartments()
compnames <- input$get.names(comps)
extant_comps <- sapply(unname(input$get.pairs()),function(x){comps[[which(compnames == x)]]})# extant compartments with multiple lineages
compnames.for.extlings <- sapply(input$get.extant_lineages(),function(x){
  x$get.location()$get.name()
})

# count the number of extant lineages in one compartment
num.extlings <- function(x){
  length(which(compnames.for.extlings==x))
}

# draw waiting time for one piece of population growth dynamics
wait.time <- function(k, alpha, beta){
  u <- runif(1, 0, 1)
  (1-(1-u)^(beta/choose(k,2)))*alpha/beta
}

# popn.growth <- sapply(extant_comps,function(x){
#   x$get.type()$get.popn.growth.dynamics()
# })

get.min.waittime <- function(extant_comps){
  waiting.times <- vector()
  for (comp in extant_comps){
    # get popn.growth.dynamics
    popn.growth <- comp$get.type()$get.popn.growth.dynamics()
    # get name of the compartment
    name <- comp$get.name()
    time <- 0
    # iterate through the pieces of population growth dynamics to draw waiting time
    for (i in 1:nrow(popn.growth)){
      piece <- popn.growth[i,]
      wait <- wait.time(num.extlings(name), piece$intercept, piece$slope)
      time = time+wait
      if (time > piece$end){
        time <- piece$end
        # and draw waiting time for next piece
      }
      else {
        waiting.times <- c(waiting.times,time)
        names(waiting.times)[[length(waiting.times)]] <- name
        break
      }
    }
    # take the shortest waiting time as the outcome
    waiting.times
  }
  return(min(waiting.times))
}
