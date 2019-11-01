#' wait.time
#' 
#' Draw waiting time to next coalescent event for the current interval
#' of a piecewise linear model describing the population 
#' growth dynamics over time.
#' Inverse cumulative function from Romero-Severson et al., 2017 
#' https://doi.org/10.1534/genetics.117.300284
#' 
#' @param k: number of extant lineages that could potentially coalesce
#' @param t1: start time in reverse, relative to intercept
#' @param alpha: intercept, population size at t=0
#' @param beta: slope of linear growth interval
#' 
#' @return
#'   Random variate from waiting time distribution.
wait.time <- function(k, t1, alpha, beta){
  u <- runif(1)
  (1-(1-u)^(beta/choose(k,2)))*(alpha+beta*t1)/beta
}


#' calc.coal.wait.times
#' 
#' Draw waiting times for all compartments that have two or more extant lineages.
#' 
#' @param model: object of class 'Model' or 'Run'
#' @param current.time: current simulation time for inner tree
#' @return Named numeric vector of waiting times per compartment
#' 
#' @examples 
#' # this model specifies 10 compartments with 3 lineages each,
#' # and constant coalescent rate 1.0
#' path <- system.file('extdata', 'SI.yaml', package='twt')
#' settings <- yaml.load_file(path)
#' mod <- MODEL$new(settings)
#' 
#' # these should average 0.333
#' calc.coal.wait.times(mod, 0)
#' 
#' @export
calc.coal.wait.times <- function(model, current.time, dynamic=FALSE){
  # retrieve all compartments
  comps <- c(model$get.compartments(), model$get.unsampled.hosts())
  compnames <- names(comps)
  
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
  counts <- sapply(comps, function(x) num.ext.lineages(x))
  ext.comps <- comps[counts >= 2]
  
  # calculate waiting times per Compartment
  waiting.times <- vector()
  
  for (comp in ext.comps) {
    k <- num.ext.lineages(comp)  # >=2 at this point
    this.type <- comp$get.type()
    
    if (is.null(this.type$get.death.rate.distr()) || 
        is.null(this.type$get.popn.growth.dynamics())) {
      
      # assume constant rate of coalescence
      c.rate <- this.type$get.coalescent.rate()
      waiting.times <- c(waiting.times, rexp(n=1, rate=choose(k,2)*c.rate))
      names(waiting.times)[length(waiting.times)] <- comp$get.name()
    }
    
    else {
      # population growth dynamics with user-specified piecewise linear model
      this.name <- comp$get.name()
      popn.growth <- this.type$get.popn.growth.dynamics()  # the model!
      
      # time is measured in reverse relative to start of simulation at t=0
      if (is.null(comp$get.branchingtime())) {
        # index case, branching time not determined by transmission from other case
        # arbitrarily set br. time to maximum time at limit of population growth model
        infection.time <- popn.growth[nrow(popn.growth), 'startTime']
      } else {
        infection.time <- comp$get.branching.time()
      }
      
      # amount of time until we reach start of infection (backwards)
      delta.t = overall.t <- infection.time - current.time
      
      # index to pieces within the remaining time interval
      piece.rows <- which(popn.growth[, 'startTime'] < delta.t)
      if (length(piece.rows) == 1) {
        # only one interval left - cast row vector as matrix
        pieces <- as.matrix(t(popn.growth[piece.rows, ]))
        row.names(pieces) <- 1
      } else {
        pieces <- popn.growth[piece.rows, ]
      }
      
      for (i in nrow(pieces):1) {
        piece <- pieces[i, ]
        wait <- wait.time(num.ext.lineages(this.name), piece['intercept'], piece['slope'])
        delta.t <- delta.t - wait
        
        if (delta.t < piece['startTime']) {
          # waiting time exceeds limit of this piece
          delta.t <- piece['startTime']
          if (delta.t == 0) {
            # past limit of initial piece - either a different within-host event occurs or 
            # we apply this compartment's bottleneck
            waiting.times <- c(waiting.times, overall.t)
            names(waiting.times)[[length(waiting.times)]] <- this.name
            break
          }
          # otherwise, go to next piece
        } else {
          # waiting time within this piece
          cumul.wait.time <- overall.t - delta.t  # relative to current time
          waiting.times <- c(waiting.times, cumul.wait.time)
          names(waiting.times)[[length(waiting.times)]] <- this.name
          break
        }
      }
      
    }  # end else
    
    # go to next compartment
  }
  
  return(waiting.times)
}

