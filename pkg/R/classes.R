library(R6)


# CompartmentType
#  Defines a type of compartment such as a class of host
#  individual (risk group) with a specific transmission rate
CompartmentType  <- R6Class("CompartmentType",
  public = list(
    name = NULL,
    unsampled = NULL,
    susceptible = NULL,
    branching.rates = NULL,
    effective.size = NULL,
    bottleneck.size = NULL,
    migration.rates = NULL,
    popn.growth.dynamics = NULL,
    
    initialize = function(name=NA, unsampled = NA,
                          susceptible=NA, branching.rates=NA,
                          effective.size=NA, bottleneck.size=NA,
                          migration.rates=NA, popn.growth.dynamics=NA) {
      self$name <- name
      self$unsampled <- unsampled
      self$susceptible <- susceptible
      self$branching.rates <- branching.rates               # named vector of transmission rates corresponding to different Compartment objects
      self$effective.size <- effective.size
      self$bottleneck.size <- bottleneck.size
      self$migration.rates <- migration.rates               # named vector of migration rates of different Compartments
      self$popn.growth.dynamics <- popn.growth.dynamics
    },
    
    get.bottleneck.size = function() {
      self$bottleneck.size
    },
    
    get.name = function() {
      self$name
    },
    
    get.unsampled = function() {
      self$unsampled
    },
    
    get.susceptible = function() {
      self$susceptible
    },
    
    get.branching.rates = function() {
      self$branching.rates
    },
    
    get.branching.rate = function(name.type) {
      self$branching.rates[[name.type]]
    },
    
    get.migration.rates = function() {
      self$migration.rates
    },
    
    get.migration.rate = function(name.type) {
      self$migration.rates[[name.type]]
    },
    
    get.popn.growth.dynamics = function() {
      self$popn.growth.dynamics
    },
    
    set.bottleneck.size = function(new.size) {
      self$bottleneck.size <- new.size
    }, 
    
    set.name = function(new.name) {
      self$name <- new.name
    },
    
    set.unsampled = function(new.no) {
      self$unsampled <- new.no
    },
    
    set.susceptible = function(new.no) {
      self$susceptible <- new.no
    },
    
    set.branching.rate = function(name.type, new.rate) {
      self$branching.rates[[name.type]] <- new.rate
    },
    
    set.migration.rate = function(name.type, new.rate) {
      self$migration.rates[[name.type]] <- new.rate
    }
    
  ),
  private = list()
)




# Compartment
Compartment <- R6Class("Compartment",
  public = list(
    name = NULL,
    type = NULL,          # reference to CompartmentType object
    source = NULL,
    branching.time = NULL,
    sampling.time = NULL,
    
    initialize = function(name=NA, type=NA, source=NA, branching.time=NA, sampling.time=NA) {
      self$name <- name
      self$type <- type
      self$source <- source
      self$branching.time <- branching.time
      self$sampling.time <- sampling.time
    },
    
    get.name = function() {
      self$name
    },
    
    get.type = function() {
      self$type
    },
    
    get.source = function() {
      self$source
    },
    
    get.branching.time = function() {
      self$branching.time
    },
    
    get.sampling.time = function() {
      self$sampling.time
    },
    
    set.type = function(new.type) {
      self$type <- new.type
    },
    
    set.source = function(new.source) {
      self$source <- new.source
    },
    
    set.branching.time = function(new.branching.time) {
      self$branching.time <- new.branching.time
    }
  
  ),
  private = list()
)




# Lineage
Lineage <- R6Class("Lineage",
  public = list(
    name = NULL,
    type = NULL,            # potential reference to LineageType object
    sampling.time = NULL,
    location = NULL,
    initialize = function(name=NA, type=NA, sampling.time=NA, location=NA) {
      self$name <- name
      self$type <- type
      self$sampling.time <- sampling.time
      self$location <- location
    },
    
    get.name = function() {
      self$name
    },
    
    get.type = function() {
      self$type
    },
    
    get.sampling.time = function() {
      self$sampling.time
    },
    
    get.location = function() {
      self$location
    },
    
    set.sampling.time = function(new.sampling.time) {
      self$sampling.time <- new.sampling.time
    },
    
    set.location = function(locationList, new.locationName) {
      new.locationObj <- locationList[[ which(sapply(locationList, function(x) {x$get.name()}) == new.locationName) ]]
      self$location <- new.locationObj
    }
    
  ),
  private = list()
)




# Event Logger (tracks migration, transmission, coalescent, and bottleneck events)
EventLogger <- R6Class("EventLogger",
  public = list(
    events = NULL,
    events.noncumul = NULL,
    
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
                          ) 
    {
      self$events <- events
      self$events.noncumul <- events.noncumul
    },
    
    
    
    get.all.events = function(cumulative=TRUE) {
      if (nrow(self$events.noncumul) == 0) {cat('No events to display.')}
      else {
        if (cumulative) { 
          private$generate.cumul.eventlog(self$events.noncumul)                        # default eventlog shows cumulative time b/c more user friendly
          self$events
        } else {
          self$events.noncumul
        }
      }
    },
    
    
    
    get.events = function(event.type) {
      eventList <- self$events.noncumul
      indices <- which(eventList$event.type == event.type)
      if (length(indices) != 0) {
        as.data.frame(t(sapply(indices, function(x) {eventList[x,]})))
      } else {NULL}
      
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
      self$events.noncumul <- rbind(self$events.noncumul, nonCumulEvent, stringsAsFactors=F)
    },
    
    
    
    clear.events = function() {
      self$events = self$events.noncumul <- data.frame(event.type=character(),
                                                       time=numeric(),
                                                       lineage1=character(),
                                                       lineage2=character(),
                                                       compartment1=character(),
                                                       compartment2=character()
                                                       )
    }
    
  ),
  private = list(

    generate.cumul.eventlog = function(noncumul.eventlog) {
      self$events <- data.frame(event.type=character(),
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
        events <- noncumul.eventlog$get.events(event.name)
        
        if (event.name == 'coalescent') {
          # set of events and their times remain the same, since they are inputted as cumulative time already
          
        } else {
          # find the longest path from root to tip
          root <- setdiff(events$compartment2, events$compartment1)
          tips <- setdiff(events$compartment1, events$compartment2)
          maxTime <- private$find.max.time(events, root, tips)
          
          # trace from root to tips and calculate all subsequent cumulative times based on maxTime (end time of transmission tree simulation)
          private$generate.events(events, maxTime, root)
        }
      })
      # return newly populated eventlog containing cumulative times
      self$events
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
          parent_time <- events[ which(events$compartment1 == parent), 'time']
          return (node_time + calc.cumul.time(parent, parent_time))
        }
      } 
      
      # for each tip, trace back and calculate resulting cumulative time at the root
      for (tip in tips) {
        tip_time <- events[ which(events$compartment1 == tip), 'time']
        cumul.time <- calc.cumul.time(tip, tip_time)
        if (cumul.time > maxTime) {
          maxTime <- cumul.time
        }
      }
      maxTime
    },
    
    
    
    generate.events = function(events, maxTime, root, tips) {
      
      # recursive function to generate cumulative times for each individual event
      generate.indiv.event <- function(node, parent_time) {
        if (node %in% tips) {
          return (NULL)
        } else {
          nodeEvents <- events[ which(events$compartment2 == node), ]
          if (length(nodeEvents) == 0) {
            return(NULL)
          } else {
            for (x in seq_along(nodeEvents)) {
              childEvent <- nodeEvents[x,]
              #delta_time <- events[ which(events$compartment1 == node), 'time' ]
              childEvent['time'] <- parent_time - childEvent['time']
              self$events <- rbind(self$events, childEvent, stringsAsFactors=F)
              # same for the descendants
              generate.indiv.event(as.character(childEvent['compartment1']), childEvent['time'])
            }
          }
          
        }
      }
      
      
      rootEvents <- events[ which(events$compartment2 == root), ] 
      for (x in seq_along(rootEvents)) {
        parentEvent <- rootEvents[x,]
        parentEvent['time'] <- maxTime - parentEvent['time']
        self$events <- rbind(self$events, parentEvent, stringsAsFactors=F)
        # do the same for the descendants
        generate.indiv.event(as.character(parentEvent['compartment1']), parentEvent['time'])
      }
    
    }
    
    
  )
)

