# inner-tree simulation
inner.tree <- function(model, eventlog) {
  # collect extant lineages at time t=0
  current.time <- 0.0
  extant.lineages <- model$get.extant_lineages(current.time)
  num.extant <- length(extant.lineages)
  
  while (length(extant.lineages) > 1) {
    # calculate waiting times for coalescent events for each compartment with 2 or more lineages
    coal.wait.times <- calc.coal.wait.times(model, current.time)
    # calculate total migration rate across all compartments at a given time
    mig.rate <- calc.migration.rates(model, current.time)
    # transmission times
    transm.times <- eventlog$get.events('transmission')$time
    # record number of transmission events already included in simulation at this current.time
    num.transm.occurred <- length(transm.times <= current.time)
    
    if (length(coal.wait.times) == 0) {
      if (mig.rate == 0) {
        # if no coalescence or migration events possible at this point in time
        # move up to the earliest transmission event in the EventLogger
        next.time <- min(transm.times)
        current.time <- next.time
        next
      } else {
        # no coalescence possible, but migration events are possible
        
      }
    }
    
    # retrieve the minimum waiting time of the calculated coalescent event waiting times
    coal.time <- min(coal.wait.times)
    # calculate a waiting time to the next migration event
    mig.time <- rexp(n=1, rate=mig.rate)
    
    # take minimum waiting time
    new.time <- min(mig.time, coal.time)
    
    if (length(transm.times <= new.time) > num.transm.occurred) {
      # if minimum waiting time exceeds a transmission event not previously included 
      # (in the count of transmission events that have been recorded to have occurred)
      # update the current time to add the waiting time and start again
      
      all.new.transm <- which(transm.times <= new.time)
      new.time <- transm.times[length(all.new.transm)]
      
    } else if (length(model$get.extant_lineages(new.time)) > num.extant) {
      # OR other lineage(s) become(s) extant before the waiting time
      # update the current time to add the waiting time and start again
      
      new.time <- new.time + current.time
      
    } else {
      
      
      
    }
    
     current.time <- new.time
     
    
    
  }
}