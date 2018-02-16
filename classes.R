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
    
    initialize = function(name=NA, unsampled = NA,
                          susceptible=NA, branching.rates=NA,
                          effective.size=NA, bottleneck.size=NA,
                          migration.rates=NA) {
      self$name <- name
      self$unsampled <- unsampled
      self$susceptible <- susceptible
      self$branching.rates <- branching.rates               # named vector of transmission rates corresponding to different Compartment objects
      self$effective.size <- effective.size
      self$bottleneck.size <- bottleneck.size
      self$migration.rates <- migration.rates               # named vector of migration rates of different Compartments
    },
    
    get.bottleneck.size = function() {
      self$bottleneck.size
    },
    
    get.effective.size = function() {
      self$effective.size
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
      rate <- self$branching.rates[[name.type]]
      rate
    },
    
    get.migration.rates = function() {
      self$migration.rates
    },
    
    get.migration.rate = function(name.type) {
      rate <- self$migration.rates[[name.type]]
      rate
    },
    
    set.bottleneck.size = function(new.size) {
      self$bottleneck.size <- new.size
    }, 
    
    set.effective.size = function(new.size) {            # could also set Ne as a univariate function in mode character()
      self$effective.size <- new.size
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
    initialize = function(events=data.frame(event.type=character(),
                                            time=numeric(),
                                            lineage1=character(),
                                            lineage2=character(),
                                            compartment1=character(),
                                            compartment2=character()
                                            ),
                          event.noncumul=data.frame(event.type=character(),
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
      as.data.frame(t(sapply(indices, function(x) {eventList[x,]})))
    },
    
    add.event = function(name, time, obj1, obj2, obj3, cumulative=TRUE) {
      if (cumulative) {
        # adding an event with cumulative time
        if (tolower(name) == 'transmission' || tolower(name) == 'migration') {
          new.event <- list(event.type=name, time=time, lineage1=obj1, lineage2=NA, compartment1=obj2, compartment2=obj3)
        } else if (tolower(name) == 'coalescent' || tolower(name) == 'bottleneck') {
          new.event <- list(event.type=name, time=time, lineage1=obj1, lineage2=obj2, compartment1=obj3, compartment2=NA)
        } else {
          stop(name, 'is not an event name.')
        }
        self$events <- rbind(self$events, new.event, stringsAsFactors=F)
        
        # need to break down into individual delta t and store in self$events.noncumul
        
      } else {
        # adding an event with only an individual delta t
        if (tolower(name) == 'transmission' || tolower(name) == 'migration') {
          new.event <- list(event.type=name, time=time, lineage1=obj1, lineage2=NA, compartment1=obj2, compartment2=obj3)
        } else if (tolower(name) == 'coalescent' || tolower(name) == 'bottleneck') {
          new.event <- list(event.type=name, time=time, lineage1=obj1, lineage2=obj2, compartment1=obj3, compartment2=NA)
        } else {
          stop(name, 'is not an event name.')
        }
        self$events.noncumul <- rbind(self$events, new.event, stringsAsFactors=F)
        
        # need to add up individual delta t and convert into cumulative time and store into self$events
        
      }
      
    },
    
    clear.events = function() {
      self$events <- data.frame(event.type=character(),
                                time=numeric(),
                                lineage1=character(),
                                lineage2=character(),
                                compartment1=character(),
                                compartment2=character()
                                )
    }
    
  ),
  private = list()
)

