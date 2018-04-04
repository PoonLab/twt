
# draw waiting time for one piece of population growth dynamics
wait.time <- function(k, alpha, beta){
  u <- runif(1, 0, 1)
  (1-(1-u)^(beta/choose(k,2)))*alpha/beta
}

waittimes.for.allextcomps <- function(model, current.time){
  # draw waiting times for all compartments that have two or more extant lineages
  # @param model = MODEL object
  # @return vector of waiting times
  comps <- model$get.compartments()
  compnames <- model$get.names(comps)
  
  # extant compartments with multiple lineages
  extant_comps <- unique(sapply(
    unname(model$get.pairs()),
    function(x) {
      comps[[which(compnames == x)]]
    }
  ))
  
  # compartment names for extant lineages
  compnames.for.extlings <- sapply(
    model$get.extant_lineages(),
    function(x) {
      x$get.location()$get.name()
    }
  )
  
  # count the number of extant lineages in one compartment
  num.extlings <- function(x) {
    length(which(compnames.for.extlings==x))
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
      # get the infection time of the compartment
      infect.time <- comp$get.branching.time()
      # time = compartment infection time - current simulation time
      time <- infect.time-current.time
      x <- which(popn.growth[,'time']< time)
      pieces <- popn.growth[x,]
      # iterate through the pieces of population growth dynamics to draw waiting time
      for (i in 1:nrow(pieces)){
        wait <- wait.time(num.extlings(name), piece['intercept'], piece['slope'])
        time = time-wait
        #if waiting time exceeds the start time of the piece, move time to the start time (end time of the previous piece)
        if (time < piece['time']){
          time <- piece['time']
          waittime.for.piece <- 
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
