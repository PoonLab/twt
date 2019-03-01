sim.inner.tree <- function(model, eventlog) {
  # Simulate the coalescence of pathogen lineages within hosts (compartments)
  # AND resolve the migration events that may or may not involve sampled pathogen lineages within hosts
  # EventLogger object updated with inner tree events
  #
  # @param model: R6 object from Model$new()
  # @param eventlog: R6 EventLogger object populated by sim.outer.tree() and sim.migrations()
  
  # vector of all infected Compartments in population at time zero (most recent)
  inf <- c(model$get.compartments(), model$get.unsampled.hosts())
  inf.names <- model$get.names(inf)
  
  # retrieve transmission events and times
  transm.events <- eventlog$get.events('transmission')
  transm.times <- transm.events$time[!is.na(transm.events$time)]
  
  # retrieve migration events **to be resolved** and times
  migration.events <- eventlog$get.migration.events()
  migration.times <- migration.events$time
  
  current.time <- 0.0
  
  # collect extant lineages at time t=0 and store into a vector
  extant.lineages <- model$get.extant.lineages(current.time)
  num.extant <- length(extant.lineages)
  if (num.extant == 0) {
    stop ('There must be at least one lineage sampled at time t=0.')
  }
  
  
  while (length(extant.lineages) > 1) {
    
    # calculate waiting times for coalescent events for each compartment with 2 or more lineages
    coal.wait.times <- calc.coal.wait.times(model, current.time)
    
    # record number of transmission events already included in simulation at this current.time
    num.transm.occurred <- length(which(transm.times <= current.time))
    
    # record number of migration events already included in simulation at this current.time
    num.migrations.occurred <- length(which(migration.times <= current.time))
      
    if (length(coal.wait.times) == 0) {
      # no coalescent events possible at this point in time
      # move up to the earliest event (transmission or migration) in the EventLogger
      
      current.time <- min(c(transm.times, migration.times))
      extant.lineages <- model$get.extant.lineages(current.time)
      next
      
    } else {
      # retrieve the minimum waiting time of the calculated coalescent event waiting times
      chosen.time <- min(coal.wait.times)
      
    }
  
    new.time <- chosen.time + current.time
    
    
    # check to see if the minimum waiting time is not exceeded by other fixed events:
    
    if (length(which(transm.times <= new.time)) > num.transm.occurred) {
      # if minimum waiting time exceeds a transmission event not previously included 
      # (in the count of transmission events that have been recorded to have occurred)
      # update the current time to the earlier transmission time (coalesc. time) of all newly included transmission events and start again
      # call bottleneck function to mass coalesce lineages in (first) new transmission event newly included
      # change the location of the lineages that 'survived' the bottleneck to the source of the transmission
      
      old.transm <- transm.times[which(transm.times <= current.time)]
      all.new.transm <- transm.times[which(transm.times <= new.time)]
      current.time <- min(union( all.new.transm, old.transm )) 
      
      transm.event <- transm.events[which(transm.times == current.time),]
      comp.2.bottle <- inf[[which(inf.names == transm.event$compartment1)]]
      comp.2.receive <- inf[[which(inf.names == transm.event$compartment2)]]
      
      survivor.lineages <- generate.bottleneck(model, eventlog, comp.2.bottle, current.time)
      new.comp.location <- comp.2.receive$get.name()
      
      # TODO: if bottleneck.size > 1, for each of the survivor lineages, the compartment needs to update its location to the source of the transmission event
      #survivor.names <- sapply(survivor.lineages, function(x) {
      #  x$set.location(inf, new.comp.location)
      #  x$get.name()
      #})
      
      # TODO: need to remove survivor lineages from comp.2.bottle, and add those same survivor lineages to new.comp.location
      #sapply(survivor.lineages, function(x) {
      #  comp.2.bottle$remove.lineage(x)
      #  comp.2.receive$add.lineage(x)
      #})
      
      # FIXME: only works for when bottleneck.size > 1
      survivor.lineages$set.location(inf, new.comp.location)
      survivor.names <- survivor.lineages$get.name()
      comp.2.bottle$remove.lineage(survivor.lineages)
      comp.2.receive$add.lineage(survivor.lineages)
      
      # the transmission event's `lineage` column can now be updated to included the names of the 'survivors'
      eventlog$modify.event(transm.event$time, survivor.names)
      extant.lineages <- model$get.extant.lineages(current.time)
      
      next
      
    } else if (length(which(migration.times <= new.time)) > num.migrations.occurred) {
      # next event is a migration event that needs to be resolved      
      
      old.migrations <- migration.times[which(migration.times <= current.time)]
      all.new.migrations <- migration.times[which(migration.times <= new.time)]
      current.time <- min(set.diff(all.new.migrations, old.migrations))
      
      migration.event <- migration.events[which(migration.times == current.time),]
      
      # first, draw a random recipient compartment with a matching recipientType of the migration event to be resolved
      possible.recipients <- sapply(inf, function(x) {
        if (x$get.type()$get.name() == migration.event$r_type) {x}
      })
      chosen.recipient <- sample(possible.recipients, 1)
      
      # does this compartment involve one of our sampled infected compartments?
      if (chosen.recipient$is.unsampled() != TRUE) {
        
        # if so, does this migration event involve one of our sampled lineages?
        # probability of the migration involving a sampled lineage is (number of extant sampled lineages / N ), where N = 1/coal.rate
        prob.samp.lin <- length(chosen.recipient$get.lineages()) * chosen.recipient$get.type()$get.coalescent.rate()
        
        # if so, draw a source compartment with a matching sourceType to the event to be resolved
        # additionally, randomly draw a sampled infected lineage to be the migrating lineage and generate the migration event
        if (runif(1) < prob.samp.lin) {
          
          possible.sources <- sapply(inf, function(y) {
            of (y$get.type()$get.name() == migration.event$s_type) {y}
          })
          chosen.source <- sample(possible.sources, 1)
         
          migrating.lineage <- sample(chosen.recipient$get.lineages(), 1)[[1]]         # draw a lineage to be migrated
          generate.migration(model, eventlog, chosen.source, chosen.recipient, migrating.lineage, migration.event$time)
        } else {
          # TODO : maybe don't "restart this simulation, rather, go straight to the next checks
          current.time <- new.time
          extant.lineages <- model$get.extant.lineages(current.time)
          
          next
        }
        
      } else {
        # TODO : maybe don't "restart this simulation, rather, go straight to the next checks
        current.time <- new.time
        extant.lineages <- model$get.extant.lineages(current.time)
        
        next
      }
      
    } else if (length(model$get.extant.lineages(new.time)) > num.extant) {
      # OR other lineage(s) become(s) extant before the waiting time
      # update the current time to add the waiting time and start again
      
      current.time <- new.time
      extant.lineages <- model$get.extant.lineages(current.time)
      next
      
    } else {
      # checks have been made, move forward with generating a coalescent event with the new time
      
      # choose a 'bin' from compartment with new.time
      coal.comp.name <- names(coal.wait.times)[which(coal.wait.times == chosen.time)]   # name of the compartment with the min coal wait time
      coal.comp <- inf[[which(inf.names == coal.comp.name)]]
      
      coal.comp.lineages <- sapply(extant.lineages, function(x){
        if (x$get.location()$get.name() == coal.comp.name) {x}
        else {NULL}
      })
      
      coal.comp.lineages[sapply(coal.comp.lineages, is.null)] <- NULL   # cleanup
      lineages.to.coalesce <- sample(coal.comp.lineages, 2)
      generate.coalescent(model, eventlog, lineages.to.coalesce, coal.comp, new.time)
      

      # issue 40: if coalescent event occurs at a transmission time, force coalescence of all other lineages at this time
      if (length(which(transm.times <= new.time)) > num.transm.occurred) {
        old.transm <- transm.times[which(transm.times <= current.time)]
        all.new.transm <- which(transm.times <= new.time)
        current.time <- min(union( all.new.transm, old.transm ))
        
        transm.event <- transm.events[which(transm.times == current.time),]
        comp.2.bottle <- inf[[which(inf.names == transm.event$compartment1)]]
        comp.2.receive <- inf[[which(inf.names == transm.event$compartment2)]]
        
        survivor.lineages <- generate.bottleneck(model, eventlog, comp.2.bottle, current.time)
        new.comp.location <- comp.2.receive$get.name()
        
        # for each of the survivor lineages, the compartment needs to update its location to the source of the transmission event
        survivor.names <- sapply(survivor.lineages, function(x) {
          x$set.location(inf, new.comp.location)
          x$get.name()
        })
        
        # need to remove survivor lineages from comp.2.bottle, and add those same survivor lineages to new.comp.location
        sapply(survivor.lineages, function(x) {
          comp.2.bottle$remove.lineage(x)
          comp.2.receive$add.lineage(x)
        })
        
        # the transmission event's `lineage` column can now be updated to included the names of the 'survivors'
        eventlog$modify.event(transm.event$time, survivor.names)
      }
      
      # update current time and extant lineages
      current.time <- new.time
      extant.lineages <- model$get.extant.lineages(current.time)
      
    }
    
  }
  
  #eventlog
}




generate.migration <- function(model, eventlog, source, recipient, migrating.lineage, migration.time) {
  # function records a migration event with a given lineage from a source to recipient migration
  #
  # @param model = MODEl object
  # @param eventlog = EventLogger object
  # @param migrating.lineage = Lineage object; to be migrated
  # @param inf = list of Compartment objects as potential migration sources
  # @param inf.names = list of Compartment object names of type character (coincides with inf)
  # @param current.time = double; time of simulation current to this function call

  l_name <- migrating.lineage$get.name()                       # lineage name
  
  outer.tree.events <- eventlog$get.events('transmission')
  outer.tree.comps <- union(outer.tree.events$compartment1, outer.tree.events$compartment2)

  if (source$get.name() %in% outer.tree.comps == FALSE) {
    # issue 32: if migration of lineages from US individual not already included in outer tree, have to graft another branch to the outer tree
    # source of migration is from US individual not already included in outer tree --> sample a time from stored vector of used & unused times
    t.times <- source$get.type()$get.transmission.times()

    # sample an available time
    avail.t.time <- t.times[ sample(which(t.times[names(t.times)==T] >= migration.time), 1) ]

    # reset the vector to make the sampled time now unavailable
    # CONFIRMED: the type's vector of available t.times is updated for every compartment under this Compartment Type
    sampled.t.time.ind <- which(t.times == avail.t.time)
    names(t.times)[sampled.t.time.ind] <- FALSE
    master.comp.type <- which(model$get.names(model$get.types()) == source$get.type()$get.name())
    model$get.types()[[master.comp.type]]$set.transmission.times(t.times)

    # generate new transmission event for branch to be grafted onto tree
    eventlog$add.event('transmission', avail.t.time, l_name, NA, recipient$get.name(), source$get.name())
  }

  # add lineage to source compartment
  source$add.lineage(migrating.lineage)
  source$add.lineage.pairs(model, migrating.lineage)

  # remove lineage from recipient compartment
  recipient$remove.lineage(migrating.lineage)
  recipient$remove.lineage.pairs(migrating.lineage)

  # set location of migrating lineage to source compartment
  migrating.lineage$set.location(inf, source$get.name())

  # add migration event to EventLogger
  eventlog$add.event('migration', current.time, l_name, NA, recipient$get.name(), source$get.name())

}



generate.coalescent <- function(model, eventlog, lineages.to.coalesce, coalescing.comp, current.time) {
  # function records the coalescent of 2 given lineages currently extant in given Compartment into an ancestral lineage
  #
  # @param model = MODEL object
  # @param eventlog = EventLogger object
  # @param lineages.to.coalesce = vector of Lineage objects of length 2
  # @param coalescing.comp = Compartment object
  # @param current.time = double; time of simulation current to this function call
  
  names.coal.lineages <- sapply(lineages.to.coalesce, function(x) x$get.name())
  
  # create a new ancestral lineage 
  ancestral.lineage <- Lineage$new(name = model$get.node.ident(),           # label internal nodes iteratively from 1..inf by convention
                                   sampling.time = current.time,
                                   location = coalescing.comp)
  
  # remove lineages undergoing coalescence
  model$remove.lineage(lineages.to.coalesce[[1]])
  model$remove.lineage(lineages.to.coalesce[[2]])
  coalescing.comp$remove.lineage(lineages.to.coalesce[[1]])
  coalescing.comp$remove.lineage(lineages.to.coalesce[[2]])
  
  # add ancestral lineage resulting from coalescence
  model$add.lineage(ancestral.lineage)
  coalescing.comp$add.lineage(ancestral.lineage)
  model$update.node.ident()
  
  # remove pairs containing coalesced lineages from list of pair choices
  remove.lineage.pairs(model, lineages.to.coalesce[[1]])
  remove.lineage.pairs(model, lineages.to.coalesce[[2]])
  
  # add pairs with new ancestral lineage into list of pair choices
  add.lineage.pairs(model, ancestral.lineage)
  
  eventlog$add.event('coalescent', current.time, names.coal.lineages[1], names.coal.lineages[2], ancestral.lineage$get.name(), coalescing.comp$get.name())
  
}



generate.bottleneck <- function(model, eventlog, comp, current.time) {
  # function coalesces lineages currently extant in given Compartment to the given Compartment's bottleneck size
  # bottleneck size of Compartment is user determined by CompartmentType
  #
  # @param model = MODEL object
  # @param eventlog = EventLogger object
  # @param comp = Compartment object
  # @param current.time = double; time of simulation current to this function call
  # @return comp.lineages = lineages that 'survive' the bottleneck on towards source of transmission in coalescent time
  
  bottleneck.size <- comp$get.type()$get.bottleneck.size()
  
  #while (length(comp$get.lineages()) > bottleneck.size) {
  #  lineages.to.coalesce <- sample(comp$get.lineages(), 2)
  #  generate.coalescent(model, eventlog, lineages.to.coalesce, comp, current.time)
  #}
  
  # TODO: if bottleneck size > 1, randomly separate lineages into "bottleneck groups"
  bottleneck.lineages <- paste0(sapply(comp$get.lineages(), function(x) {x$get.name()}), collapse=',')
  
  sapply(comp$get.lineages(), function(x) {
    # remove lineages undergoing bottleneck
    model$remove.lineage(x)
    comp$remove.lineage(x)
    # remove pairs containing bottleneck lineages from list of pair choices
    remove.lineage.pairs(model, x)
  })
  
  # create a new ancestral lineage 
  ancestral.lineage <- Lineage$new(name = model$get.node.ident(),           # label internal nodes iteratively from 1..inf by convention
                                   sampling.time = current.time,
                                   location = comp)
  
  # add ancestral lineage resulting from bottleneck
  model$add.lineage(ancestral.lineage)
  comp$add.lineage(ancestral.lineage)
  model$update.node.ident()
  
  # add pairs with new ancestral lineage into list of pair choices (unnecessary for bottleneck.size of 1)
  # add.lineage.pairs(model, ancestral.lineage)
  
  eventlog$add.event('bottleneck', current.time, bottleneck.lineages, NA, ancestral.lineage$get.name(), comp$get.name())
  
  # TODO: if bottleneck size > 1, return list of lineages that 'survive' the bottleneck (coalescent time) onwards to source of transmission 
  ancestral.lineage
}




remove.lineage.pairs <- function(model, lineage) {
  # this function removes lineage pairs in MODEL obj attr `choices`: list of lineage pair choices to coalesce
  #
  # @param model = MODEL object, one per simulation
  # @param lineage = Lineage object, to be removed from list of lineage pair choices
  
  # extract all the pairs
  current.pairs <- model$get.pairs()
  
  # narrow down to only the compartment with lineage of interest
  compname <- lineage$get.location()$get.name()
  
  if (length(current.pairs) != 0) {
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
}



add.lineage.pairs <- function(model, lineage) {
  # this function adds lineage pairs in MODEL obj attr `choices`: list of lineage pair choices to coalesce
  #
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
