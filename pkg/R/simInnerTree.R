# inner-tree simulation
sim.inner.tree <- function(model, eventlog) {
  
  inf <- c(model$get.compartments(), model$get.unsampled.hosts())
  inf.names <- model$get.names(inf)
  
  # collect extant lineages at time t=0
  transm.events <- eventlog$get.events('transmission')
  # transmission times
  transm.times <- transm.events$time[!is.na(transm.events$time)]
  current.time <- 0.0
  extant.lineages <- model$get.extant.lineages(current.time)
  num.extant <- length(extant.lineages)
  if (num.extant == 0) {
    stop ('There must be at least one lineage sampled at time t=0.')
  }
  
  
  while (length(extant.lineages) > 1) {
    # calculate waiting times for coalescent events for each compartment with 2 or more lineages
    coal.wait.times <- calc.coal.wait.times(model, current.time)
    # calculate total migration rate across all compartments at a given time
    mig.rate <- calc.migration.rates(model, current.time)
    # record number of transmission events already included in simulation at this current.time
    num.transm.occurred <- length(which(transm.times <= current.time))
    
    if (mig.rate == 0) {
      # no migration events possible
      mig.time <- -1
      if (length(coal.wait.times) == 0) {
        # no coalescence or migration events possible at this point in time
        # move up to the earliest transmission event in the EventLogger
        current.time <- min(transm.times)
        extant.lineages <- model$get.extant.lineages(current.time)
        next
      } else {
        # retrieve the minimum waiting time of the calculated coalescent event waiting times
        new.time <- min(coal.wait.times)
      }
    } else {
      # calculate a waiting time to the next migration event
      mig.time <- rexp(n=1, rate=mig.rate)
      if (length(coal.wait.times) == 0) {
        # no coalescence events possible; accept the migration time as the new time
        new.time <- mig.time
      } else {
        # take the minimum waiting time out of all coalescent wait times and migration time
        new.time <- min(c(mig.time, coal.wait.times))
      }
    }
    
    
    # check to see if the minimum waiting time is not exceeded by other fixed events:
    
    if (length(which(transm.times <= new.time)) > num.transm.occurred) {
      # if minimum waiting time exceeds a transmission event not previously included 
      # (in the count of transmission events that have been recorded to have occurred)
      # update the current time to the earlier transmission time (coalesc. time) of all newly included transmission events and start again
      # call bottleneck function to mass coalesce lineages in (first) new transmission event newly included
      # change the location of the lineages that 'survived' the bottleneck to the source of the transmission
      
      old.transm <- which(transm.times <= current.time)
      all.new.transm <- which(transm.times <= new.time)
      current.time <- min(setdiff( transm.times[all.new.transm], transm.times[old.transm] )) 
      
      transm.event <- transm.events[which(transm.times == current.time),]
      comp.2.bottle <- inf[[which(inf.names == transm.event$compartment1)]]
      
      survivor.lineages <- generate.bottleneck(model, eventlog, comp.2.bottle, current.time)
      sapply(survivor.lineages, function(x) model$add.lineage(x))
      new.comp.location <- transm.event$compartment2
      
      # for each of the survivor lineages, the compartment needs to update its location to the source of the transmission event
      survivor.names <- sapply(survivor.lineages, function(x) {
        x$set.location(inf, new.comp.location)
        x$get.name()
      })
      
      # the transmission event's `lineage` column can now be updated to included the names of the 'survivors'
      eventlog$modify.event(transm.event$time, survivor.names)
      
      # update extant lineages
      extant.lineages <- model$get.extant.lineages(current.time)
      
      next
      
    } else if (length(model$get.extant.lineages(new.time)) > num.extant) {
      # OR other lineage(s) become(s) extant before the waiting time
      # update the current time to add the waiting time and start again
      
      current.time <- current.time + new.time
      # update extant lineages
      extant.lineages <- model$get.extant.lineages(current.time)
      next
      
    } else {
      # checks have been made, move forward with generating a migration or coalescent event with the new time
      
      if (new.time == mig.time) {
        # next event is a migration event; now draw migrated lineage, and source and recipient compartments
        
        migrating.lineage <- sample(extant.lineages, 1)              # draw a lineage to be migrated
        l_name <- migrating.lineage$get.name()                       # lineage name
        r_comp_name <- migrating.lineage$get.location()              # recipient compartment Lineage migrated to (forward time)
        
        r_ind <- which(inf.names == r_comp_name)
        # exclude r_comp, and exclude any source comp currently w/ only one lineage (otherwise it would be considered a transmission)
        filtered.inf <- sapply(inf[-r_ind], function(i) {length(i$get.lineages()) > 1})   
        s_comp <- sample(filtered.inf, 1)                            # source compartment Lineage migrated from (forward time)
        s_name <- s_comp$get.name()
        
        # issue 32: if migration of lineages from US individual not already included in outer tree, have to graft another branch to the outer tree
        outer.tree.events <- eventlog$get.events('transmission')
        outer.tree.comps <- union(outer.tree.events$compartment1, outer.tree.events$compartment2)
        if (s_name %in% outer.tree.comps == FALSE) {
          # source of migration is from US individual not already included in outer tree --> sample a time from stored vector of used & unused times
          t.times <- s_comp$get.type()$get.transmission.times()
          
          # sample an available time
          avail.t.time <- sample(which(t.times[names(t.times)==T] >= mig.time), 1)
          
          # reset the vector to make the sampled time now unavailable
          # CONFIRMED: the type's vector of available t.times is updated for every compartment under this Compartment Type
          sampled.t.time.ind <- which(t.times == avail.t.time)
          names(t.times)[sampled.t.time.ind] <- FALSE
          master.comp.type <- which(model$get.types() == s_comp$get.type())
          model$get.types()[[master.comp.type]]$set.transmission.times(t.times)
          
          # generate new transmission event for branch to be grafted onto tree
          eventlog$add.event('transmission', avail.t.time, l_name, NA, r_comp_name, s_comp_name)
          
        }
        
        migrating.lineage$set.location(inf, s_name)
        
        # add migration event to EventLogger
        eventlog$add.event('migration', mig.time, l_name, NA, r_comp_name, s_comp_name)
        
      } else {
        # next event is a coalescent event; now choose a 'bin' from compartment with new.time
        
        coal.comp.name <- names(coal.wait.times)[which(coal.wait.times == new.time)]   # name of the compartment with the min coal wait time
        coal.comp <- inf[[which(inf.names == coal.comp.name)]]
        coal.comp.lineages <- sapply(extant.lineages, function(x){
          if (x$get.location()$get.name() == coal.comp.name) {x}
          else {NULL}
        })
        coal.comp.lineages[sapply(coal.comp.lineages, is.null)] <- NULL   # cleanup
        lineages.to.coalesce <- sample(coal.comp.lineages, 2)
        names.coal.lineages <- sapply(lineages.to.coalesce, function(x) x$get.name())
        
        # create a new ancestral lineage
        ancestral.lineage <- Lineage$new(name = model$get.node.ident(),           # label internal nodes iteratively from 1..inf by convention
                                         sampling.time = (current.time + new.time),
                                         location = coal.comp)
        model$add.lineage(ancestral.lineage)
        model$update.node.ident()
        
        # remove pairs containing coalesced lineages from list of pair choices
        remove.lineage.pairs(model, lineages.to.coalesce[[1]])
        remove.lineage.pairs(model, lineages.to.coalesce[[2]])
        model$remove.lineage(lineages.to.coalesce[[1]])
        model$remove.lineage(lineages.to.coalesce[[2]])
        
        # add pairs with new ancestral lineage into list of pair choices
        add.lineage.pairs(model, ancestral.lineage)
        
        eventlog$add.event('coalescent', current.time, names.coal.lineages[1], names.coal.lineages[2], ancestral.lineage$get.name(), coal.comp.name)
      }
      
      current.time <- current.time + new.time
      
      
      # issue 40: if coalescent event occurs at a transmission time, force coalescence of all other lineages at this time
      if (length(which(transm.times <= current.time)) > num.transm.occurred) {
        old.transm <- which(transm.times <= current.time)
        all.new.transm <- which(transm.times <= new.time)
        current.time <- min(setdiff( transm.times[all.new.transm], transm.times[old.transm] )) 
        
        transm.event <- transm.events[which(transm.times == current.time),]
        comp.2.bottle <- inf[[which(inf.names == transm.event$compartment1)]]
        new.comp.location <- transm.event$compartment2
        
        survivor.lineages <- generate.bottleneck(model, eventlog, comp.2.bottle, current.time)
        sapply(survivor.lineages, function(x) model$add.lineage(x))
        
        
        # for each of the survivor lineages, the compartment needs to update its location to the source of the transmission event
        survivor.names <- sapply(survivor.lineages, function(x) {
          x$set.location(inf, new.comp.location)
          x$get.name()
        })
        
        # need to remove survivor lineages from comp.2.bottle, and add those same survivor lineages to new.comp.location
        sapply(survivor.lineages, function(x) {
          comp.2.bottle$remove.lineage(x)
          new.comp.location$add.lineage(x)
        })
        
        
        # the transmission event's `lineage` column can now be updated to included the names of the 'survivors'
        eventlog$modify.event(transm.event$time, survivor.names)
      }
      
      
      # update extant lineages
      extant.lineages <- model$get.extant.lineages(current.time)
      
    }
    
  }
  
  eventlog
}



generate.bottleneck <- function(model, eventlog, comp, current.time) {
  # function coalesces lineages currently extant in given Compartment to the given Compartment's bottleneck size
  # bottleneck size of Compartment is user determined by CompartmentType
  # @param model = MODEL object
  # @param eventlog = EventLogger object
  # @param comp = unique Compartment object
  # @param current.time = time of simulation current to this function call
  # @return comp.lineages = lineages that 'survive' the bottleneck on towards source of transmission in coalescent time
  
  bottleneck.size <- comp$get.type()$get.bottleneck.size()
  extant.lineages <- model$get.extant.lineages(current.time) # reason to not pass extant.lineages directly is b/c there could be additional extant lineages updated alongside current.time
  comp.lineages <- unlist(sapply(extant.lineages, function(x) {
    if (x$get.location()$get.name() == comp$get.name()) {x}
    else {NULL}
  }))
  
  while (length(comp.lineages) > bottleneck.size) {
    lineages.to.coalesce <- sample(comp.lineages, 2)
    names.coal.lineages <- sapply(lineages.to.coalesce, function(x) x$get.name())
    
    # create a new ancestral lineage
    ancestral.lineage <- Lineage$new(name = model$get.node.ident(),       # all lineages forced to coalesce in bottleneck have same node identification
                                     sampling.time = current.time,
                                     location = comp)
    model$add.lineage(ancestral.lineage)
    
    # remove pairs containing coalesced lineages from list of pair choices
    remove.lineage.pairs(model, lineages.to.coalesce[[1]])
    remove.lineage.pairs(model, lineages.to.coalesce[[2]])
    model$remove.lineage(lineages.to.coalesce[[1]])
    model$remove.lineage(lineages.to.coalesce[[2]])
   
    # add pairs with new ancestral lineage into list of pair choices
    add.lineage.pairs(model, ancestral.lineage)
    
    eventlog$add.event('coalescent', current.time, names.coal.lineages[1], names.coal.lineages[2], ancestral.lineage$get.name(), comp$get.name())
    
    # update while loop dependent condition: `comp.lineages`
    shortened.list <- unlist(sapply(comp.lineages, function(x) {
      if (x$get.name() %in% names.coal.lineages) {NULL}
      else {x}
    }))
    comp.lineages <- c( shortened.list, ancestral.lineage )
  }
  
  # after all lineages under forced coalescence labeled w/ same node identification, update to new node identity for next internal node
  model$update.node.ident()          
  
  # return list of lineages that 'survive' the bottleneck (coalescent time) onwards to source of transmission 
  comp.lineages
}




remove.lineage.pairs <- function(model, lineage) {
  # this function removes lineage pairs in MODEL obj attr `choices`: list of lineage pair choices to coalesce
  # @param model = MODEL object, one per simulation
  # @param lineage = Lineage object, to be removed from list of lineage pair choices
  
  # extract all the pairs
  current.pairs <- model$get.pairs()
  
  # narrow down to only the compartment with lineage of interest
  compname <- lineage$get.location()$get.name()
  comp.lineage.pairs <- unlist(sapply(1:length(current.pairs), function(x) {
    if (current.pairs[[x]] == compname) {names(current.pairs)[[x]]}
    else {NULL}
  }))
  
  # narrow down to only ones with given lineage in the name
  sapply(comp.lineage.pairs, function(y) {
    # split each pair of names and look for the lineage name
    split.names <- unlist(strsplit(y, '\\,'))
    if (lineage$get.name() %in% split.names) {
      # remove this pair of names from list attr `choices`
      model$remove.pair(split.names[1], split.names[2])
    }
  })
}



add.lineage.pairs <- function(model, lineage) {
  # this function adds lineage pairs in MODEL obj attr `choices`: list of lineage pair choices to coalesce
  # @param model = MODEL object
  # @param lineage = Lineage object to be added to list of lineage pair choices
  
  host.comp <- lineage$get.location()$get.name()
  # find other lineages in this (host.comp) location
  other.lineages <- unlist(sapply(model$get.lineages(), function(x) {
    if (x$get.location()$get.name() == host.comp) {
      if (x$get.name() == lineage$get.name()) {NULL}
      else {x$get.name()}
    } else {NULL}
  }))
  
  # generate pairs with this new lineage and other lineages
  sapply(other.lineages, function(y) {
    model$add.pair(lineage$get.name(), y, host.comp)
  })
}