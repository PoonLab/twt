# inner-tree simulation
inner.tree <- function(model, eventlog) {
  # collect extant lineages at time t=0
  current.time <- 0.0
  extant.lineages <- model$get.extant_lineages(current.time)
  num.extant <- length(extant.lineages)
  if (num.extant == 0) {
    stop ('There must be at least one lineage sampled at time t=0.')
  }
  
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
      # calculate a waiting time to the next migration event
      mig.time <- rexp(n=1, rate=mig.rate)
      if (mig.rate == 0) {
        # if no coalescence or migration events possible at this point in time
        # move up to the earliest transmission event in the EventLogger
        current.time <- min(transm.times)
        next
      } else {
        # no coalescence possible, but migration events are possible
        new.time <- mig.time
      }
    } else {
      # retrieve the minimum waiting time of the calculated coalescent event waiting times
      coal.time <- min(coal.wait.times)
      # take minimum waiting time
      new.time <- min(mig.time, coal.time)
    }
    
    # checks for if the minimum waiting time is not exceeded by other fixed events:
    
    if (length(transm.times <= new.time) > num.transm.occurred) {
      # if minimum waiting time exceeds a transmission event not previously included 
      # (in the count of transmission events that have been recorded to have occurred)
      # update the current time to the earlier transmission time (coalesc. time) and start again
      all.new.transm <- which(transm.times <= new.time)
      current.time <- transm.times[length(all.new.transm)]
      next
      
    } else if (length(model$get.extant_lineages(new.time)) > num.extant) {
      # OR other lineage(s) become(s) extant before the waiting time
      # update the current time to add the waiting time and start again
      current.time <- current.time + new.time
      next
      
    } else {
      # checks have been made, move forward with generating a migration or coalescent event with the new time
      if (new.time == mig.time) {
        # next event is a migration event; now draw migrated lineage, and source and recipient compartments
        migrating.lineage <- sample(extant.lineages, 1)              # draw a lineage to be migrated
        l_name <- migrating.lineage$get.name()                       # lineage name
        r_name <- migrating.lineage$get.location()                   # recipient compartment lineage migrated to (forward time)
        s_name <- sample()$get.name()                                # source compartment lineage migrated from (forward time)
        
        migrating.lineage$set.location(s_name)
        eventlog$add.event('migration', next.time, l_name, r_name, s_name)
      } else {
        # next event is a coalescent event; now choose a 'bin'
        eventlog$add.event('coalescent', next.time, obj1, obj2, obj3)
      }
      
      current.time <- current.time + new.time
      
    }
    
    
  }
}