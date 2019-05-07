# draw waiting time for one piece of population growth dynamics
# inverse cumulative function (Romero-Severson et al., 2017)
wait.time <- function(k, alpha, beta){
  u <- runif(1, 0, 1)
  (1-(1-u)^(beta/choose(k,2)))*alpha/beta
}

calc.coal.wait.times <- function(model, current.time, dynamic=FALSE){
  # draw waiting times for all compartments that have two or more extant lineages
  # @param model = MODEL object
  # @param current.time = current time of the simulation of inner tree
  # @return waiting.times = vector of waiting times
  comps <- model$get.compartments()
  compnames <- model$get.names(comps)
  
  # retrieves compartment names for extant lineages
  ext.lineages.compnames <- sapply(
    model$get.extant.lineages(current.time),
    function(x) {
      x$get.location()$get.name()
    }
  )
  
  # counts the number of extant lineages in one compartment
  # calculated for parameter `k` in function `wait.time`
  num.ext.lineages <- function(x) {
    length(which(ext.lineages.compnames==x$get.name()))
  }
  
  # retrieves compartments with multiple extant lineages
  ext.comps <- sapply(comps, function(x) {if (num.ext.lineages(x) >= 2) x})
  ext.comps[sapply(ext.comps, is.null)] <- NULL
  
  # calculate waiting times per Compartment
  waiting.times <- vector()
  for (comp in ext.comps) {
    
    compname <- comp$get.name()
    
    # determine whether to use `death.rate.distr` and `popn.growth.dynamics` to populate coalescent waiting times, or uniform coalescent rate
    if (is.null(comp$get.type()$get.death.rate.distr()) == F && is.null(comp$get.type()$get.popn.growth.dynamics()) == F) {
      
      # `death.rate.distr` and `popn.growth.dynamics` are specified in YAML
      # with a death rate, varying coalescent rates are permitted because it should be possible to draw a waiting time to infection of the root compartment
      
      # retrieve user-specified popn.growth.dynamics for this compartment
      # NOTE: `popn.growth.dynamics` is user-specified in FORWARD-TIME
      # RECALL: Compartment `get.branching.time` are specified in COALESCENT (BACKWARDS) TIME
      popn.growth <- comp$get.type()$get.popn.growth.dynamics()       
      
      if (is.null(comp$get.branching.time())) {
        # infection time is unknown for the compartment that started the epidemic -- could be arbitrarily large in COALESCENT (BACKWARDS) TIME
        # assume that initial infection time is at the maximal time before it hits the final piece (when population size becomes constant) issue #46
        infection.time <- popn.growth[nrow(popn.growth), 'startTime']     # Re-implementing after incorporating death dynamics (dated: November 26/2018)
      } else {
        infection.time <- comp$get.branching.time()
      }
      
      # delta time = compartment infection time - current simulation time
      delta.t = overall.t <- infection.time - current.time       
      piece.rows <- which(popn.growth[,'startTime'] < delta.t)
      if (length(piece.rows) == 1) {
        # obtain all the pieces of the popn.growth.dynamics functions that are valid for current simulation time
        pieces <- as.matrix(t(popn.growth[piece.rows, ]))                  
        rownames(pieces) <- 1
      } else {
        pieces <- popn.growth[piece.rows, ]         
      }
      
      for (i in nrow(pieces):1){
        # iterate through the valid pieces of population growth dynamics to draw waiting time
        # in COALESCENT TIME, start at the current piece first (this would be the latest piece in FORWARD TIME)
        # work backwards to earlier pieces in the popn.growth.dynamic functions if needed
        
        piece <- pieces[i,]
        wait <- wait.time(num.ext.lineages(compname), piece['intercept'], piece['slope'])
        delta.t <- delta.t - wait
        
        # if waiting time exceeds the start time of the piece, move delta.t to the start time (end time of the previous piece)
        if (delta.t < piece['startTime']){
          delta.t <- piece['startTime']
          if (delta.t == 0) {
            waiting.times <- c(waiting.times, overall.t)                    # wait time 'maxed out', add total waiting time from current time
            names(waiting.times)[[length(waiting.times)]] <- compname       # associate waiting times w/ their compartment name
          }
          # go to next piece and draw new waiting time
        } else {
          cumul.wait.time <- overall.t - delta.t
          waiting.times <- c(waiting.times, cumul.wait.time)                # add the waiting time of current compartment to the vector of waiting times
          names(waiting.times)[[length(waiting.times)]] <- compname         # associate the waiting times with their compartments
          break
        }
      }
      
    } else {
      
      # either `death.rate.distr` and/or `popn.growth.dynamics` not specified in YAML
      # assume that coalescent rates must be uniform (constant within-host population sizes)
      if (num.ext.lineages(comp) > 0) {
        coal.rate <- comp$get.type()$get.coalescent.rate()
        
        wait.time <- rexp(n=1,rate=coal.rate)
        waiting.times <- c(waiting.times, wait.time)
        names(waiting.times)[[length(waiting.times)]] <- compname
      }
      
    }
    
  }
  
  waiting.times
}

