waittimes.for.allextcomps <- function(model){
  # draw waiting times for all compartments that have two or more extant lineages
  # @param model = MODEL object
  # @return vector of waiting times
  comps <- model$get.compartments()
  compnames <- model$get.names(comps)
  extant_comps <- sapply(unname(model$get.pairs()),function(x){comps[[which(compnames == x)]]})# extant compartments with multiple lineages
  compnames.for.extlings <- sapply(model$get.extant_lineages(),function(x){
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
  
  get.waittimes <- function(extant_comps){
    # initialize a vector for waiting times
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
        #if waiting time exceeds the end time of the piece, move time to the end time
        if (time > piece$end){
          time <- piece$end
          # and draw waiting time for next piece
        }
        else {
          # add the waiting time of current compartment to the vector of waiting times
          waiting.times <- c(waiting.times,time)
          # associate the waiting times with their compartments
          names(waiting.times)[[length(waiting.times)]] <- name
          break
        }
      }
    }
    return(waiting.times)
  }
  get.waittimes(extant_comps)
}
