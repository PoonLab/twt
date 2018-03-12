library(R6)


# CompartmentType
#  Defines a type of compartment such as a class of host
#  individual (risk group) with a specific transmission rate
CompartmentType  <- R6Class("CompartmentType",
  public = list(
    initialize = function(name=NA, unsampled = NA,
                          susceptible=NA, branching.rates=NA,
                          effective.size=NA, bottleneck.size=NA,
                          migration.rates=NA, popn.growth.dynamics=NA) {
      private$name <- name
      private$unsampled <- unsampled
      private$susceptible <- susceptible
      private$branching.rates <- branching.rates               # named vector of transmission rates corresponding to different Compartment objects
      private$effective.size <- effective.size
      private$bottleneck.size <- bottleneck.size
      private$migration.rates <- migration.rates               # named vector of migration rates of different Compartments
      private$popn.growth.dynamics <- popn.growth.dynamics
    },
    
    get.bottleneck.size = function() {
      private$bottleneck.size
    },
    
    get.name = function() {
      private$name
    },
    
    get.unsampled = function() {
      private$unsampled
    },
    
    get.susceptible = function() {
      private$susceptible
    },
    
    get.branching.rates = function() {
      private$branching.rates
    },
    
    get.branching.rate = function(name.type) {
      private$branching.rates[[name.type]]
    },
    
    get.migration.rates = function() {
      private$migration.rates
    },
    
    get.migration.rate = function(name.type) {
      private$migration.rates[[name.type]]
    },
    
    get.popn.growth.dynamics = function() {
      private$popn.growth.dynamics
    }
    
  ),
  private = list(
    name = NULL,
    unsampled = NULL,
    susceptible = NULL,
    branching.rates = NULL,
    effective.size = NULL,
    bottleneck.size = NULL,
    migration.rates = NULL,
    popn.growth.dynamics = NULL
  )
)




# Compartment
Compartment <- R6Class("Compartment",
  public = list(
    initialize = function(name=NA, type=NA, source=NA, branching.time=NA, sampling.time=NA) {
      private$name <- name
      private$type <- type
      private$source <- source
      private$branching.time <- branching.time
      private$sampling.time <- sampling.time
    },
    
    get.name = function() {
      private$name
    },
    
    get.type = function() {
      private$type
    },
    
    get.source = function() {
      private$source
    },
    
    get.branching.time = function() {
      private$branching.time
    },
    
    get.sampling.time = function() {
      private$sampling.time
    },
    
    set.type = function(new.type) {
      private$type <- new.type
    },
    
    set.source = function(new.source) {
      private$source <- new.source
    },
    
    set.branching.time = function(new.branching.time) {
      private$branching.time <- new.branching.time
    }
  
  ),
  private = list(
    name = NULL,
    type = NULL,          # reference to CompartmentType object
    source = NULL,
    branching.time = NULL,
    sampling.time = NULL
  )
)




# Lineage
Lineage <- R6Class("Lineage",
  public = list(
    initialize = function(name=NA, type=NA, sampling.time=NA, location=NA) {
      private$name <- name
      private$type <- type
      private$sampling.time <- sampling.time
      private$location <- location
    },
    
    get.name = function() {
      private$name
    },
    
    get.type = function() {
      private$type
    },
    
    get.sampling.time = function() {
      private$sampling.time
    },
    
    get.location = function() {
      private$location
    },
    
    set.location = function(locationList, new.locationName) {
      new.locationObj <- locationList[[ which(sapply(locationList, function(x) {x$get.name()}) == new.locationName) ]]
      private$location <- new.locationObj
    }
    
  ),
  private = list(
    name = NULL,
    type = NULL,            # potential reference to LineageType object
    sampling.time = NULL,
    location = NULL
  )
)




# Event Logger (tracks migration, transmission, and coalescent events) (bottleneck events are logged as coalescent events)
EventLogger <- R6Class("EventLogger",
  public = list(

    initialize = function(events = data.frame(event.type=character(),
                                              time=numeric(),
                                              lineage1=character(),
                                              lineage2=character(),
                                              compartment1=character(),
                                              compartment2=character()
                                             ),
    events.noncumul = data.frame(event.type=character(),
                                 time=numeric(),
                                 lineage1=character(),
                                 lineage2=character(),
                                 compartment1=character(),
                                 compartment2=character()
                                )
    ){
      private$events <- events
      private$events.noncumul <- events.noncumul
    },
   
   
    get.all.events = function(cumulative=TRUE) {
      if (nrow(private$events.noncumul) == 0) {cat('No events to display.')}
      else {
        if (cumulative) { 
          private$.generate.cumul.eventlog(private$events.noncumul)                        # default eventlog shows cumulative time b/c more user friendly
        } else {
          private$events.noncumul
        }
      }
    },
   
   
    get.events = function(event.type, cumulative=TRUE) {
      if (cumulative) {
        eventList <- private$.generate.cumul.eventlog(private$events.noncumul)
      } else {
        eventList <- private$events.noncumul
      }
      indices <- which(eventList$event.type == event.type)
      if (length(indices) != 0) {
        as.data.frame(t(sapply(indices, function(x) {eventList[x,]})))
      } else {
        cat('No events of type "', event.type, '".')
      }
     
    },
   
   
    add.event = function(name, time, obj1, obj2, obj3) {
      # @param name = event type of either transmission, migration, or coalescent
      # @param time = NON-CUMULATIVE time that event has occurred between two compartments in a transmission/migration event,
      # or CUMULATIVE time that event has occurred in a coalescent event   -------------MAY CHANGE IN FUTURE
      # @param obj1 = name of a lineage object in mode character
      if (tolower(name) == 'transmission' || tolower(name) == 'migration') {
        nonCumulEvent <- list(event.type=name, time=time, lineage1=obj1, lineage2=NA, compartment1=obj2, compartment2=obj3)
      } else if (tolower(name) == 'coalescent') {
        nonCumulEvent <- list(event.type=name, time=time, lineage1=obj1, lineage2=obj2, compartment1=obj3, compartment2=NA)
      }
      private$events.noncumul <- rbind(private$events.noncumul, nonCumulEvent, stringsAsFactors=F)
    },
   
   
    clear.events = function() {
      private$events = private$events.noncumul <- data.frame(event.type=character(),
                                                       time=numeric(),
                                                       lineage1=character(),
                                                       lineage2=character(),
                                                       compartment1=character(),
                                                       compartment2=character()
      )
    }
   
  ),
  private = list(
    events = NULL,
    events.noncumul = NULL,
    
    .generate.cumul.eventlog = function(noncumul.eventlog) {
      private$events <- data.frame(event.type=character(),
                                time=numeric(),
                                lineage1=character(),
                                lineage2=character(),
                                compartment1=character(),
                                compartment2=character()
      )
      event.types <- unique(noncumul.eventlog$event.type)        # up to 3 different event types
      
      sapply(event.types, function(event.name) {
        # for each type of event (transmission, migration, and/or coalescent)
        # retrieve set of events of that event.name type
        events <- noncumul.eventlog[ which(noncumul.eventlog$event.type == event.name), ]
        
        if (event.name == 'transmission') {
          # find the longest path from root to tip
          root <- setdiff(events$compartment2, events$compartment1)
          tips <- setdiff(events$compartment1, events$compartment2)
          maxTime <- private$find.max.time(events, root, tips)
          
          # trace from root to tips and calculate all subsequent cumulative times based on maxTime (end time of transmission tree simulation)
          private$events <- rbind(private$events, private$generate.events(events, maxTime, root, tips), stringsAsFactors=F)
          
        } else {
          # set of events and their times remain the same, since they are inputted as cumulative time already
          NULL
        }
      })
      private$events
      
    },
   
    find.max.time = function(events, root, tips) {
      # @oaram events = set of events of a particular type (ie. transmission, migration, coalescent)
      # @return maxTime = maximum cumulative time of longest path from root to tip; the 'time span' of the transmission tree simulation
      maxTime <- 0
     
      # recursive function to calculate cumulative time
      calc.cumul.time <- function(node, node_time) {
        parent <- events[ which(events$compartment1 == node), 'compartment2']
        if (parent == root) {
          return (node_time)
        } else {
          parent_time <- as.numeric(events[ which(events$compartment1 == parent), 'time'])
          return (node_time + calc.cumul.time(parent, parent_time))
        }
      } 
     
      # for each tip, trace back and calculate resulting cumulative time at the root
      for (tip in tips) {
        tip_time <- as.numeric(events[ which(events$compartment1 == tip), 'time'])
        cumul.time <- calc.cumul.time(tip, tip_time)
        if (cumul.time > maxTime) {
          maxTime <- cumul.time
        }
      }
      maxTime
    },
   
   
   
    generate.events = function(events, maxTime, root, tips) {
      
      generate.indiv.event <- function(node, parent_time) {
        # recursive function to generate cumulative times for each individual event
        # returns a childEvent or NULL to be added to the eventlog data frame
        if (node %in% tips) {
          return (NULL)
        } else {
          nodeEvents <- events[ which(events$compartment2 == node), ]
          if (length(row.names(nodeEvents)) == 0) {
            return(NULL)
          } else {
            for (x in seq_along(nodeEvents)) {
              childEvent <- nodeEvents[x,]
              childEvent['time'] <- parent_time - as.numeric(childEvent['time'])
              private$events <- rbind(private$events, childEvent, stringsAsFactors=F)
              # same for the descendants
              generate.indiv.event(as.character(childEvent['compartment1']), childEvent['time'])
            }
            return(private$events)
          }
         
        }
      }
     
      rootEvents <- events[ which(events$compartment2 == root), ] 
      for (x in seq_along(rootEvents)) {
        parentEvent <- rootEvents[x,]
        parentEvent['time'] <- maxTime - as.numeric(parentEvent['time'])
        private$events <- rbind(private$events, parentEvent, stringsAsFactors=F)
        # do the same for all of the descendants
        generate.indiv.event(as.character(parentEvent['compartment1']), parentEvent['time'])
      }
     
      indices <- grep('NA', row.names(private$events), ignore.case=T, invert=T)
      match.noncumul.ordering <- order(as.numeric(row.names(private$events[indices,])))
      private$events[indices,][match.noncumul.ordering,]
     
    }
                         
 )
)

