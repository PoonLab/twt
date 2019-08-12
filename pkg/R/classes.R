#' CompartmentType
#' 
#' \code{CompartmentType} is an R6 class that defines a type of compartment, 
#' such as a class of host individual (risk group) with a specific transmission 
#' rate.
#' 
#' @param name: a character string that uniquely identifies the class
#' @param unsampled: if TRUE, no Compartments of this Type can contain sampled
#' lineages (directly observed tips of the inner tree).
#' @param susceptible: if TRUE, no Compartments of this Type contain either
#' sampled or unsampled lineages.
#' @param branching.rates: a named vector of transmission rates *to* other
#' CompartmentTypes.
#' @param migration.rates: a named vector of migration rates *to* other 
#' CompartmentTypes.
#' @param bottleneck.size: the maximum number of lineages that can be transmitted
#' to a Compartment of this Type.
#' @param coalescent.rate: the rate at which lineages coalesce within a Compartment
#' of this Type.
#' @param death.rate.distr: a text expression for the waiting time distribution
#' to a death event. Required for `popn.growth.dynamics` to override 
#' `coalescent.rate`.
#' @param wait.time.distr: a text expression for the waiting time distribution
#' to a transmission event.
#' @param popn.growth.dynamics: a text expression for population growth dynamics
#' in forward time. If not NULL, can override `coalescent.rate`.
#' @param transmission.times: numeric vector of transmission event times, 
#' populated by `outer.tree.sim` from class parameters.
#' 
#' @export
CompartmentType  <- R6Class("CompartmentType",
  public = list(
    initialize = function(name=NA, unsampled = NA,
                          susceptible=NA, branching.rates=NA,
                          migration.rates=NA, bottleneck.size=NA,
                          coalescent.rate=NA, death.rate.distr=NA, wait.time.distr=NA,
                          popn.growth.dynamics=NA, transmission.times=NA) {
      private$name <- name
      private$unsampled <- unsampled
      private$susceptible <- susceptible
      private$branching.rates <- branching.rates               # named vector of transmission rates corresponding to different Compartment objects
      private$migration.rates <- migration.rates
      private$bottleneck.size <- bottleneck.size
      private$coalescent.rate <- coalescent.rate               # named vector of migration rates of different Compartments
      private$death.rate.distr <- death.rate.distr
      private$wait.time.distr <- wait.time.distr
      private$popn.growth.dynamics <- popn.growth.dynamics
      private$transmission.times <- transmission.times         # populated after outer.tree.sim, tracked used and unused for migration events in inner.tree.sim
    },
    
    # accessor functions
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
    
    get.coalescent.rate = function() {
      private$coalescent.rate
    },
    
    get.death.rate.distr = function() {
      private$death.rate.distr
    },
    
    get.wait.time.distr = function() {
      private$wait.time.distr
    },
    
    get.popn.growth.dynamics = function() {
      private$popn.growth.dynamics
    },
    
    get.transmission.times = function() {
      private$transmission.times
    },
    
    set.transmission.times = function(vector.transm.times) {
      private$transmission.times <- vector.transm.times
    },
    
    set.migration.rate = function(recipient.type, new.migr.rate) {
      private$migration.rates[[recipient.type]] <- new.migr.rate
    }
    
  ),
  private = list(
    name = NULL,
    unsampled = NULL,
    susceptible = NULL,
    branching.rates = NULL,
    migration.rates = NULL,
    bottleneck.size = NULL,
    coalescent.rate = NULL,
    death.rate.distr = NULL,
    wait.time.distr = NULL,
    popn.growth.dynamics = NULL,
    transmission.times = NULL
  )
)




#' Compartment
#' 
#' \code{Compartment} is an R6 class for objects that represent the units of an
#' outer tree simulation, such as a host individual or deme.
#' 
#' @param name: a character string that uniquely identifies the Compartment
#' @param type: a reference to a CompartmentType object
#' @param source: a reference to another Compartment from which a Lineage was 
#' transmitted to this Compartment
#' @param branching.time: stores the origin time of this Compartment, which 
#' corresponds to a branching event in the "outer" tree.
#' @param unsampled: if TRUE, then any Lineage carried by this Compartment is not
#' directly observed, i.e., it does not represent a tip in the "inner" tree.
#' 
#' @export 
Compartment <- R6Class("Compartment",
  public = list(
    initialize = function(name=NA, type=NA, source=NA, branching.time=NA, unsampled=FALSE, lineages=list()) {
      private$name <- name
      private$type <- type
      private$source <- source
      private$branching.time <- branching.time
      private$unsampled <- unsampled                   # attr req later when identifying new US Comps to be promoted in mig events
      private$lineages <- lineages
    },
    
    # accessor functions
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
    },
    
    is.unsampled = function() {
      private$unsampled
    },
    
    get.lineages = function() {
      private$lineages
    },
    
    add.lineage = function(new.lineage) {
      private$lineages[[length(private$lineages)+1]] <- new.lineage
    },
    
    remove.lineage = function(ex.lineage) {
      lin.ind <- which(sapply(private$lineages, function(x){x$get.name() == ex.lineage$get.name()}))
      private$lineages <- private$lineages[-lin.ind]
    }
  
  ),
  private = list(
    name = NULL,
    type = NULL,          # reference to CompartmentType object
    source = NULL,
    branching.time = NULL,
    unsampled = NULL,
    lineages = NULL
  )
)




#' Lineage
#' 
#' \code{Lineage} is an R6 class for objects that represent pathogen lineages
#' that are carried by Compartments and which comprise the "inner" tree of the 
#' simulation.
#' 
#' @param name: a character string that uniquely identifies the Lineage
#' @param type: a reference to an object of class LineageType (not yet implemented)
#' @param sampling.time: the time that the Lineage was sampled; left to NA for 
#' unsampled Lineages
#' @param location: a reference to a Compartment object
#' 
#' @export
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
    
    get.type = function() {                                     # in the future, will be a pointer to a LineageType object
      private$type
    },
    
    get.sampling.time = function() {
      private$sampling.time
    },
    
    get.location = function() {
      private$location
    },
    
    set.location = function(locationList, new.locationName) {
      new.locationObj <- locationList[[ 
        which(sapply(locationList, function(x) {x$get.name()}) == new.locationName) 
        ]]
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




#' EventLogger
#' 
#' \code{EventLogger} is an R6 class for an object that tracks migration, 
#' transmission, and coalescent events.  Note that bottleneck events are logged as 
#' coalescent events.
#' 
#' @param events: a data frame where each row represents an event with a time 
#' stamp in forward time.
#' @param events.noncumul: optionally, a data frame where events are stamped 
#' with the waiting time since the preceding event
#' @param migration.events.storage: a data frame of migration events, returned by
#' `simMigrations.R:.calc.migration.events()`.
#' 
#' @export
EventLogger <- R6Class("EventLogger",
  public = list(
    # FIXME: it would be more sensible to have one "events" argument and a 
    #        boolean for whether the data frame is cumulative time or not.
    initialize = function(events = data.frame(event.type=character(),
                                              time=numeric(),
                                              lineage1=character(),
                                              lineage2=character(),
                                              compartment1=character(),
                                              compartment2=character(),
                                              stringsAsFactors = FALSE
                                             ),
                          events.noncumul = data.frame(event.type=character(),
                                                       time=numeric(),
                                                       lineage1=character(),
                                                       lineage2=character(),
                                                       compartment1=character(),
                                                       compartment2=character(),
                                                       stringsAsFactors = FALSE
                                                      ),
                          migration.events.storage = data.frame(stringsAsFactors = FALSE)
    ){
      private$events <- events
      private$events.noncumul <- events.noncumul
      private$migration.events.storage <- migration.events.storage
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
        # cat('No events of type "', event.type, '".\n')
        NULL
      }
     
    },
    
    add.event = function(type, time, line1, line2, comp1, comp2) {
      # @param type: event type, one of 'transmission', 'migration', 'coalescence',
      # or 'bottleneck'.
      # @param time: CUMULATIVE time that event has occurred between two compartments 
      # in a transmission/migration/coalescent event
      
      if (is.element(type, c('transmission', 'migration'))) {
        e <- list(event.type=type, time=time, lineage1=line1, lineage2=NA,
                  compartment1=comp1, compartment2=comp2)
      } else if (is.element(type, c('coalescent', 'bottleneck'))) {
        e <- list(event.type=type, time=time, lineage1=line1, lineage2=line2,
                  compartment1=comp1, compartment2=NA)
      } else {
        stop("Error, unrecognized type argument in add.event()")
      }
      
      private$events <- rbind(private$events, e, stringsAsFactors=F)
    } 
   
   
    clear.events = function() {
      private$events = private$events.noncumul <- data.frame(event.type=character(),
                                                             time=numeric(),
                                                             lineage1=character(),
                                                             lineage2=character(),
                                                             compartment1=character(),
                                                             compartment2=character(),
                                                             stringsAsFactors = FALSE
                                                            )
    },
    
    
    modify.event = function(transmission.time, lineages) {
      # when inner tree simulation has reached a transmission event, need to fill in the lineage column w/ the lineages that are present at transmission time
      transmission.events <- self$get.events('transmission')                 # in the case of bottleneck events, will have same time so have to isolate transmission event times
      index <- which(transmission.events$time == transmission.time)
      rowname <- rownames(transmission.events)[index]
      eventlog.index <- which(rownames(self$get.all.events()) == rowname)
      private$events[eventlog.index, 'lineage1'] <- lineages
    },
    
    get.migration.events = function() {
      private$migration.events.storage
    },
    
    store.migration.events = function(migration.events) {
      private$migration.events.storage <- migration.events
    },
    
    get.fixed.samplings = function() {
      # retrieves the fixed sampling times of the tips of a MODEL object
      private$fixed.samplings.storage
    },
    
    store.fixed.samplings = function(model.fixed.samplings) {
      # stores the fixed sampling times of the tips of a MODEL object
      private$fixed.samplings.storage <- model.fixed.samplings
    }
    
  ),
  private = list(
    events = NULL,
    events.noncumul = NULL,
    migration.events.storage = NULL,
    fixed.samplings.storage = NULL,
    
    generate.noncumul.eventlog = function(cumul.eventlog) {
      # generates an event log with non-cumulative times of events
      # @param cumul.eventlog = an event log with cumulative times of events
      private$events.noncumul <- data.frame(event.type=character(),
                                            time=numeric(),
                                            lineage1=character(),
                                            lineage2=character(),
                                            compartment1=character(),
                                            compartment2=character(),
                                            stringsAsFactors = FALSE
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
      maxRootTime <- max(rootEvents$time)
      for (x in 1:nrow(rootEvents)) {
        parentEvent <- rootEvents[x,]
        # traverse descendants
        generate.indiv.event(as.character(parentEvent['compartment1']), as.numeric(parentEvent['time']))
        # root's individualt delta t from when it was infected to when it made its first transmission is 'undefined'
        parentEvent['time'] <- maxRootTime - parentEvent['time']      # 0 or 1 by convention (see treeswithintrees closed issue #29)
        private$events.noncumul <- rbind(private$events.noncumul, parentEvent, stringsAsFactors=F)
      }
     
      indices <- grep('NA', row.names(private$events.noncumul), ignore.case=T, invert=T)
      match.cumul.ordering <- order(as.numeric(row.names(private$events.noncumul[indices,])))
      private$events.noncumul[indices,][match.cumul.ordering,]
     
    }
                         
 )
)

