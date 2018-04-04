# simulation of the inner tree
simulate <- function(model, eventlog) {
  comps <- model$get.compartments()
  compnames <- model$get.names(comps)
  
  extant_l <- model$get.extant_lineages()
  not_extant_l <- list()
  not_yet_sampled_l <- list()
  
  extant_c <- unique(sapply(extant_l,function(x){x$get.location()}))
  not_extant_c <- list()
  not_yet_sampled_c <- list()
}


# inner-tree simulation
inner.tree <- function(model, eventlog) {
  # collect extant lineages at time t=0
  extant.lineages <- model$get.extant_lineages()
  current.time <- 0.0
  
  while (length(extant.lineages) > 1) {
    # calculate waiting times for coalescent events for each compartment with 2 or more lineages
    coal.wait.times <- waittimes.for.allextcomps()
    # calculate total migration rate across all compartments at a given time
    mig.rate <- calc.migration.rates(model)
    
    if (length(coal.wait.times) == 0) {
      if (mig.rate == 0) {
        # if no coalescence or migration events possible at this point in time
        # move up to the earliest transmission event in the EventLogger
        
        next
      }
    }
    
    # retrieve the minimum waiting time of the calculated coalescent event waiting times
    coal.time <- min(coal.wait.times)
    # calculate a waiting time to the next migration event
    mig.time <- rexp(n=1, rate=mig.rate)
    
    # take minimum waiting time
    
    # if minimum waiting time exceeds a transmission event 
    # OR another lineage becomes extant before the waiting time
    # update the current time to add the waiting time and start again
  }
}