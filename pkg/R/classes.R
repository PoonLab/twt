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
    initialize = function(name=NA, type=NA, source=NA, branching.time=NA) {
      private$name <- name
      private$type <- type
      private$source <- source
      private$branching.time <- branching.time
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
    branching.time = NULL
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
      if (nrow(private$events) == 0) {cat('No events to display.')}
      else {
        if (cumulative) {
          private$events                    # default eventlog shows cumulative time b/c more user friendly
        } else {
          private$generate.noncumul.eventlog(private$events)
          private$events.noncumul
        }
      }
    },
   
   
    get.events = function(event.type, cumulative=TRUE) {
      if (cumulative) {
        eventList <- private$events[ which(private$events$event.type == event.type), ]
      } else {
        private$generate.noncumul.eventlog(private$events)
        eventList <- private$events.noncumul[ which(private$events.noncumul$event.type == event.type), ]
      }
      if (nrow(eventList) != 0) {
        eventList
      } else {
        cat('No events of type "', event.type, '".')
      }
     
    },
   
   
    add.event = function(name, time, obj1, obj2, obj3) {
      # @param name = event type of either transmission, migration, or coalescent
      # @param time = CUMULATIVE time that event has occurred between two compartments in a transmission/migration/coalescent event
      # @param obj1 = name of a lineage object in mode character
      # @param obj2 = name of a recipient compartment object in mode character for a transmission or migration event
                 # OR name of a second lineage object in mode character for a coalescent event
      # @param obj3 = name of a source compartment object in mode character for a transmission or migration event
                 # OR name of a compartment object in mode character for a coalescent event
      if (tolower(name) == 'transmission' || tolower(name) == 'migration') {
        CumulEvent <- list(event.type=name, time=time, lineage1=obj1, lineage2=NA, compartment1=obj2, compartment2=obj3)
      } else if (tolower(name) == 'coalescent') {
        CumulEvent <- list(event.type=name, time=time, lineage1=obj1, lineage2=obj2, compartment1=obj3, compartment2=NA)
      }
      private$events <- rbind(private$events, CumulEvent, stringsAsFactors=F)
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
    
    generate.noncumul.eventlog = function(cumul.eventlog) {
      # generates an event log with non-cumulative times of events
      # @param cumul.eventlog = an event log with cumulative times of events
      private$events.noncumul <- data.frame(event.type=character(),
                                            time=numeric(),
                                            lineage1=character(),
                                            lineage2=character(),
                                            compartment1=character(),
                                            compartment2=character()
                                            )
      event.types <- unique(cumul.eventlog$event.type)   # up to 3 different event types
      
      sapply(event.types, function(event.name) {
        # for each type of event (transmission, migration, and/or coalescent)
        # retrieve set of events of that event.name type
        events <- self$get.events(event.name)
        if (nrow(events) == 0) {
          NULL
        } else {
          if (event.name == 'transmission' || event.name == 'migration') {
            root <- setdiff(events$compartment2, events$compartment1)
            tips <- setdiff(events$compartment1, events$compartment2)
            
            # trace from root to tips and calculate all subsequent non-cumulative times based on maxTime (cumulative time of root)
            private$events.noncumul <- rbind(private$events.noncumul, private$generate.events(events, root, tips), stringsAsFactors=F)
          } else if (event.name == 'coalescent') {
            
          }
        }
      })
      
    },
   
   
   
    generate.events = function(events, root, tips) {
      
      # inner recursive helper function
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
            for (x in 1:nrow(nodeEvents)) {
              childEvent <- nodeEvents[x,]
              # traverse descendants
              generate.indiv.event(as.character(childEvent['compartment1']), as.numeric(childEvent['time']))
              childEvent['time'] <- parent_time - as.numeric(childEvent['time'])
              private$events.noncumul <- rbind(private$events.noncumul, childEvent, stringsAsFactors=F)
            }
            return(private$events.noncumul)
          }
         
        }
      }
      
      # beginning of function generate.events()
      rootEvents <- events[ which(events$compartment2 == root), ] 
      for (x in 1:nrow(rootEvents)) {
        parentEvent <- rootEvents[x,]
        # traverse descendants
        generate.indiv.event(as.character(parentEvent['compartment1']), as.numeric(parentEvent['time']))
        # root's individualt delta t from when it was infected to when it made its first transmission is 'undefined'
        parentEvent['time'] <- NaN
        private$events.noncumul <- rbind(private$events.noncumul, parentEvent, stringsAsFactors=F)
      }
     
      indices <- grep('NA', row.names(private$events.noncumul), ignore.case=T, invert=T)
      match.cumul.ordering <- order(as.numeric(row.names(private$events.noncumul[indices,])))
      private$events.noncumul[indices,][match.cumul.ordering,]
     
    }
                         
 )
)

