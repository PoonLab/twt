# inner-tree simulation
inner.tree <- function(model, eventlog) {
  # collect extant lineages at time t=0
  transm.events <- eventlog$get.events('transmission')
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
    transm.times <- transm.events$time
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
      # update the current time to the earlier transmission time (coalesc. time) of all newly included transmission events and start again
      # call bottleneck function to mass coalesce lineages in (first) new transmission event newly included
      
      all.new.transm <- which(transm.times <= new.time)
      current.time <- min(setdiff( transm.times[all.new.transm], transm.times )) 
      comp.2.bottle <- transm.events[which(transm.times == current.time),]$compartment1
      generate.bottleneck(model, eventlog, comp.2.bottle, current.time)
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
        
        # add migration event to EventLogger
        eventlog$add.event('migration', next.time, l_name, r_name, s_name)
        
      } else {
        # next event is a coalescent event; now choose a 'bin' from compartment with next.time
        
        coal.comp.name <- names(coal.time)   # name of the compartment with the min coal wait time
        coal.comp.lineages <- sapply(extant.lineages, function(x){
          if (x$get.location()$get.name() == coal.comp.name) {x$get.name()}
          else {NULL}
        })
        lineages.to.coalesce <- sample(coal.comp.lineages, 2)
        
        # create a new ancestral lineage
        ancestral.lineage <- Lineage$new(name = paste(lineages.to.coalesce, sep=';'),
                                         sampling.time = next.time,
                                         location = coal.comp.name)
        model$add.lineage(ancestral.lineage)
        
        # remove pairs containing coalesced lineages from list of pair choices
        # add pairs with new ancestral lineage into list of pair choices
        
        # add coalescent event to EventLogger
        eventlog$add.event('coalescent', next.time, lineages.to.coalesce[1], lineages.to.coalesce[2], ancestral.lineage$get.name(), coal.comp.name)
      }
      
      current.time <- current.time + new.time
      
    }
    
    
  }
}




generate.bottleneck <- function(model, eventlog, compartment.name, current.time) {
  # function coalesces lineages currently extant in given Compartment to the given Compartment's bottleneck size
  # bottleneck size of Compartment is user determined by CompartmentType
  # @param model = MODEL object
  # @param eventlog = EventLogger object
  # @param compartment.name = name of a unique Compartment object
  # @param current.time = time of simulation current to this function call
  # @return
  
  bottleneck.size <-
  extant.lineages <- model$get.extant_lineages(current.time) # the reason not passing extant.lineages directly is b/c could have additional extant lineages updated alongside current.time
  comp.lineages <- sapply(extant.lineages, function(x) {
    if (x$get.location()$get.name() == compartment.name) {x$get.name()}
    else {NULL}
  })
  
  while (length(comp.lineages) > bottleneck.size)
  
}