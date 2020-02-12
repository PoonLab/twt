#' CompartmentType
#' 
#' \code{CompartmentType} is an R6 class that defines a type of compartment, 
#' such as a class of host individual (risk group) with a specific transmission 
#' rate.
#' 
#' @param name: a character string that uniquely identifies the class
#' @param unsampled: if TRUE, no Compartments of this Type can contain sampled
#'        lineages (directly observed tips of the inner tree).
#' @param branching.rates: a named vector of transmission rates *to* other
#'        CompartmentTypes.
#' @param tranistion.rates: a named vector of transition rates to other 
#'        CompartmentTypes.
#' @param migration.rates: a named vector of migration rates *to* other 
#'        CompartmentTypes.
#' @param bottleneck.size: the maximum number of lineages that can be transmitted
#'        to a Compartment of this Type.
#' @param effective.size: reciprocal of the rate at which lineages coalesce within 
#'        a Compartment of this Type.
#' @param popn.growth.dynamics: a text expression for population growth dynamics
#'        in forward time. If not NULL, can override `effective.size`.
#' @param generation.time:  scales coalescent events to the time scale of 
#'        outer events (including transmission).
#' @param transmission.times: numeric vector of transmission event times, 
#'        populated by `outer.tree.sim` from class parameters.
#' 
#' @examples 
#' 
#' # load CompartmentTypes from a YAML object
#' path <- system.file('extdata', 'SI.yaml', package='twt')
#' settings <- yaml.load_file(path)
#' mod <- Model$new(settings)
#' mod$get.types()
#' 
#' # manually specify a CompartmentType object (usually done by YAML)
#' host <- CompartmentType$new(name='host', branching.rates=c(host=0.1), 
#' bottleneck.size=1, coalescent.rate=1)
#' 
#' @export
CompartmentType  <- R6Class(
  "CompartmentType",
  public = list(
    initialize = function(name=NA, unsampled = NA, branching.rates=NA, 
                          transition.rates=NA, migration.rates=NA, bottleneck.size=NA,
                          effective.size=NA, popn.growth.dynamics=NA, 
                          generation.time=NA, transmission.times=NA) {
      private$name <- name
      private$unsampled <- unsampled
      
      # named vector of transmission rates corresponding to different Compartment 
      # objects
      private$branching.rates <- branching.rates
      private$transition.rates <- transition.rates
      private$migration.rates <- migration.rates
      private$bottleneck.size <- bottleneck.size
      
      # named vector of migration rates of different Compartments
      private$effective.size <- effective.size
      private$popn.growth.dynamics <- popn.growth.dynamics
      private$generation.time <- generation.time
      
      # populated after outer.tree.sim, tracked used and unused for migration 
      # events in inner.tree.sim
      private$transmission.times <- transmission.times
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
    
    get.branching.rates = function() {
      private$branching.rates
    },
    
    get.branching.rate = function(current.time, name.type) {
      
      # Get timing of rate changes
      rate.changes <- as.numeric(names(private$branching.rates))
      rate.changes <- rate.changes[order(rate.changes, decreasing = T)]
      
      # Does this CompartmentType have multiple rates?
      if (length(rate.changes) > 1) {
        
        # Is there a rate for the current time?
        if (current.time %in% rate.changes) {
          private$branching.rates[[as.character(current.time)]][[name.type]]
        }
        else {
          index <- max(which(rate.changes > current.time))
          private$branching.rates[[as.character(index)]][[name.type]]
        }
      }
      else {
        private$branching.rates[[1]][[name.type]]
      }
    },
    
    get.migration.rates = function() {
      private$migration.rates
    },
    
    get.migration.rate = function(name.type) {
      private$migration.rates[[name.type]]
    },
    
    get.transition.rates = function() {
      private$transition.rates
    },
    
    get.transition.rate = function(name.type) {
      private$transition.rates[[name.type]]
    },
    
    get.effective.size = function() {
      private$effective.size
    },
    
    get.popn.growth.dynamics = function() {
      private$popn.growth.dynamics
    },
    
    get.generation.time = function() {
      private$generation.time
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
    
    branching.rates = NULL,
    transition.rates = NULL,
    migration.rates = NULL,
    
    bottleneck.size = NULL,
    effective.size = NULL,
    popn.growth.dynamics = NULL,
    generation.time = NULL,
    transmission.times = NULL
  )
)


print.CompartmentType <- function(obj) {
  cat(paste(obj$name, ":"))
}


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
#' @examples 
#' # load Compartments from a YAML object
#' path <- system.file('extdata', 'SI.yaml', package='twt')
#' settings <- yaml.load_file(path)
#' mod <- MODEL$new(settings)
#' 
#' # display first Compartment object
#' host1 <- mod$get.compartments()[[1]]
#' host1
#' 
#' # manually initialize a new Compartment object
#' hostN <- Compartment$new(name='newHost', unsampled=TRUE)
#' hostN$set.type(host1$get.type())
#' 
#' @export 
Compartment <- R6Class("Compartment",
  public = list(
    initialize = function(name=NA, type=NA, source=NA, branching.time=NA, 
                          sampling.time=NA, unsampled=FALSE, lineages=list()) {
      private$name <- name
      private$type <- type
      private$source <- source
      private$branching.time <- branching.time
      private$sampling.time <- sampling.time
      # attr req later when identifying new US Comps to be promoted in mig events
      private$unsampled <- unsampled
      private$lineages <- lineages
    },
    
    copy = function(deep=FALSE) {
      # see https://github.com/r-lib/R6/issues/110
      cloned <- self$clone(deep)  # calls deep_clone method
      if (deep) {
        # attach new Lineages to new Compartment
        for (lineage in cloned$get.lineages()) {
          lineage$set.location(cloned)
        }
      }
      cloned  # return
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
    },
    
    set.sampling.time = function() {
      private$sampling.time <- max(sapply(private$lineages, function(line) {
        line$get.sampling.time()
      }))
    },
    
    set.unsampled = function(is.unsampled) {
      private$unsampled <- is.unsampled
    },
    
    is.unsampled = function() {
      private$unsampled
    },
    
    get.lineages = function() {
      private$lineages
    },
    
    add.lineage = function(new.lineage) {
      private$lineages[[length(private$lineages)+1]] <- new.lineage
      self$set.sampling.time()
    },
    
    remove.lineage = function(ex.lineage) {
      private$lineages[[ex.lineage$get.name()]] <- NULL
      self$set.sampling.time()
    }
    
  ),
  private = list(
    name = NULL,
    type = NULL,  # reference to CompartmentType object
    source = NULL,
    branching.time = NULL,
    sampling.time = NULL,
    unsampled = NULL,
    lineages = NULL,
    
    deep_clone = function(name, value) {
      if (name == 'lineages') {
        # map deep clone to Lineage copy() method
        lapply(value, function(lineage) lineage$copy(deep=TRUE))
      } 
      else {
        value
      }
    }
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
#' @examples
#' # load Compartments from a YAML object
#' path <- system.file('extdata', 'SI.yaml', package='twt')
#' settings <- yaml.load_file(path)
#' mod <- MODEL$new(settings)
#' 
#' # display first Lineage in first Compartment
#' comp <- mod$get.compartments()[[1]]
#' comp$get.lineages()  # display all 3 lineages
#' 
#' # manually add an unsampled Lineage
#' lin <- Lineage$new(name="L0", location=comp)
#' comp$add.lineage(lin)
#' 
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
    
    parent = NULL,
   
    copy = function(deep=FALSE) {
      # see https://github.com/r-lib/R6/issues/110
      if (deep) {
        parent <- private$location
        private$location <- NULL  # temporarily erase before cloning!
      }
     
      cloned <- self$clone(deep)
     
      if (deep) {
        private$location <- parent  # restore original reference
      }
     
      cloned
    },
   
    get.name = function() {
      private$name
    },
   
    get.type = function() {  # in the future, will be a pointer to a LineageType object
      private$type
    },
   
    get.sampling.time = function() {
      private$sampling.time
    },
   
    get.location = function() {
      private$location
    },
     
    set.location = function(comp) {
      private$location <- comp
    },
   
    set.location.by.name = function(locationList, new.locationName) {
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


