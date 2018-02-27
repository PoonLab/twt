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
      self$migrations.rates[[name.type]] <- new.rate
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
    
    set.location = function(new.location) {
      self$location <- new.location
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
      if (nrow(self$events) == 0) {cat('No events to display.')}
      else {
        if (cumulative) { 
          self$events                       # default eventlog shows cumulative time b/c user can input cumulative inf.times manually, and those must be cumulative
        } else {
          self$events.noncumul
        }
      }
    },
    
    get.events = function(event.type, cumulative=TRUE) {
      if (cumulative) {
        eventList <- self$events
      } else {
        eventList <- self$events.noncumul
      }
      indices <- which(eventList$event.type == event.type)
      if (length(indices) != 0) {
        as.data.frame(t(sapply(indices, function(x) {eventList[x,]})))
      } else {NULL}
      
    },
    
    add.event = function(name, time, obj1, obj2, obj3, cumulative=TRUE) {
      current.events <- self$get.events(name)                  # retrieve all events with that particular name
      
      if (cumulative) {
        CumulTime <- time
      } else {
        nonCumulTime <- time
      }
      
      if (tolower(name) == 'transmission' || tolower(name) == 'migration') {
        if (is.null(current.events)) {
          if (cumulative) {
            nonCumulTime <- time
          } else {
            CumulTime <- time
          }
        } else if (obj2 %in% sapply(1:nrow(current.events), function(x){current.events[x,]$compartment2}) == F) {
          if (cumulative) {
            nonCumulTime <- time
          } else {
            CumulTime <- time
          }
          # if (obj3 %in% sapply(1:nrow(current.events), function(x){current.events[x,]$compartment1})) {
          #   private$modify.times(name, time, CumulTime, obj3, cumulative)
          # }
          
        } else {
          eventsAsSource <- which(sapply(1:nrow(current.events), function(x){current.events[x,]$compartment2}) == obj2)
          maxTime <- max(unlist(current.events[eventsAsSource, 'time']))
          if (cumulative) {
            nonCumulTime <- time - maxTime
          } else {
            CumulTime <- time + maxTime
          }
          # if (obj3 %in% sapply(1:nrow(current.events), function(x){current.events[x,]$compartment1})) {
          #   private$modify.times(name, time, CumulTime, obj3, cumulative)
          # }
        }
        
        CumulEvent <- list(event.type=name, time=CumulTime, lineage1=obj1, lineage2=NA, compartment1=obj2, compartment2=obj3)
        nonCumulEvent <- list(event.type=name, time=nonCumulTime, lineage1=obj1, lineage2=NA, compartment1=obj2, compartment2=obj3)
        
      } else if (tolower(name) == 'coalescent') {
        
        tCumul.lin1 <- current.events[ which(current.events$time == obj1), 'time']
        if (cumulative) {
          nonCumulTime <- time - tCumul.lin1
        } else {
          CumulTime <- time + tCumul.lin1
        }
        
        CumulEvent <- list(event.type=name, time=CumulTime, lineage1=obj1, lineage2=obj2, compartment1=obj3, compartment2=NA)
        nonCumulEvent <- list(event.type=name, time=nonCumulTime, lineage1=obj1, lineage2=obj2, compartment1=obj3, compartment2=NA)
        
      } else {
        stop(name, 'is not an event name.')
      }
      
      self$events <- rbind(self$events, CumulEvent, stringsAsFactors=F)
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
    # 
    # modify.times = function(name, time, CumulTime, sourceName, cumulative) {
    #   event.type.indices <- which(self$events$event.type == name)
    #   # now check if the source is present in the recipient list, and re-modify CumulTime or nonCumulTime of those entries
    #   
    #   entry <- which(sapply(1:nrow(self$get.events(name)), function(x){current.events[x,]$compartment1 == sourceName}))
    #   current.cumulTime <- as.numeric(self$get.events(name)[entry, 'time'])
    #   current.nonCumulTime <- as.numeric(self$get.events(name, cumulative=F)[entry, 'time'])
    #   
    #   if (cumulative) {
    #     if (current.cumulTime > time) {
    #       modifiedTime <- current.cumulTime - time
    #       if (current.nonCumulTime > modifiedTime) {
    #         self$events.noncumul[event.type.indices[entry],'time'] <- modifiedTime
    #       }
    #     }
    #   } else {
    #     if (self$get.events(name)[entry, 'compartment2'] %in% sapply(1:nrow(current.events), function(x){self$get.events(name)[x,]$compartment1}) ) {
    #       name <- self$get.events(name)[entry, 'compartment2']
    #       time <- self$get.events(name, cumulative=F)[entry, 'time']
    #       CumulTime <- self$get.events(name)[entry, 'time']
    #       sourceName <- self$get.events(name)[which(sapply(1:nrow(current.events), function(x){self$get.events(name)[x,]$compartment1})), 'compartment2']
    #       CumulTime <- modify.times(name, time, CumulTime, sourceName, cumulative)      ## will need to be RECURSIVE
    #     }
    #     modifiedTime <- CumulTime + current.nonCumulTime
    #     if (modifiedTime > current.cumulTime ) {
    #       self$events[event.type.indices[entry], 'time'] <- modifiedTime
    #     }
    #     
    #   }
    #   
    #  
    # }
    
  )
)

