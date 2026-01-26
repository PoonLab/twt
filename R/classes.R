#' Compartment
#' \code{Compartment} is an R6 class that defines a population of host 
#' individuals that share the same epidemiological state, e.g., susceptible,
#' infected and recovered.  A compartmental model defines rates of transit 
#' between compartments.  A compartment may also represent individuals with 
#' one or more shared characteristics such as a high transmission risk group, 
#' or a geographic location.
#' 
#' @param name:  character, uniquely identifies the class
#' @param rates:  list, transition rates expressed as expressions.  Dependence 
#'                on other Compartments is identified by name.  May be a 
#'                numeric value when rate is constant.
#' @param bottleneck.size:  str, an expression for the bottleneck size 
#'        distribution affecting lineages within a member of this 
#'        Compartment.
#' @param coalescence.rate:  str, an expression for the coalescent rate for 
#'        lineages within a member of this Compartment.
#' @export
Compartment <- R6Class(
  "Compartment",
  public = list(
    initialize = function(name=NA, size=NA, rates=NA, bottleneck.size=NA,
                          coalescent.rate=NA, generation.time=NA) {
      private$name <- name
      private$size <- size
      private$rates <- rates
      private$bottleneck.size <- bottleneck.size
      private$coalescent.rate <- coalescent.rate
    },
    
    # accessor functions
    get.name = function() { private$name },
    get.size = function() { private$size },
    get.rates = function() { private$rates },
    get.bottleneck.size = function() { private$bottleneck.size },
    get.coalescent.rate = function() { private$coalescent.rate }
  ),
  
  private = list(
    name = NULL,
    size = NULL,
    rates = NULL,
    bottleneck.size=NULL,
    coalescent.rate=NULL
  )
)


#' Host
#' 
#' \code{Host} is an R6 class representing a host individual carrying 
#' one or more sampled Pathogen lineages.  The transmission history among
#' hosts comprise the "outer" tree of the simulation.
#' 
#' @param name:  character, unique identifier for the host
#' @param compartment:  character, reference to Compartment the host currently
#'        belongs to.
#' @param source:  character, reference to another Host from which a Pathogen 
#'        lineage was transmitted to this Host.  May store more than one source
#'        in the case of superinfection.
#' @param transmission.time:  numeric, stores the time(s) that a Pathogen 
#'        lineage weas transmitted to this Host.  Uses a named vector to 
#'        make the association of time to source Host less ambiguous.
#' @param sampling.time:  numeric, the earliest time that a Pathogen lineage 
#'        was directly sampled from this Host.  If NA, then this is an unsampled
#'        Host, i.e., its lineages are ancestral to lineages in subsequent 
#'        hosts where they are sampled.
#' @param pathogens:  list, store objects of class Pathogen associated with 
#'        this Host
#' @export
Host <- R6Class(
  "Host",
  public = list(
    initialize = function(name=NA, compartment=NA, source=NA, 
                          transmission.time=NA, sampling.time=NA, 
                          pathogens=list()) {
      private$name <- name
      private$compartment <- compartment
      private$source <- source
      private$transmission.time <- transmission.time
      private$sampling.time <- sampling.time
      private$pathogens <- pathogens
    },
    
    copy = function(deep=FALSE) {
      # see https://github.com/r-lib/R6/issues/110
      cloned <- self$clone(deep)  # calls deep_clone method
      if (deep) {
        # attach new Pathogens to new Host
        for (pathogen in cloned$get.pathogens()) {
          pathogen$set.location(cloned)
        }
      }
      cloned  # return
    },
    
    # accessor functions
    get.name = function() { private$name },
    get.compartment = function() { private$compartment },
    
    get.source = function() { private$source },
    set.source = function(new.source) {
      private$source <- new.source
    },
    
    get.transmission.time = function() { private$transmission.time },
    set.transmission.time = function(new.transmission.time) {
      private$transmission.time <- new.transmission.time
    },
    
    get.sampling.time = function() { private$sampling.time },
    set.sampling.time = function() {
      private$sampling.time <- min(sapply(private$pathogens, function(line) {
        line$get.sampling.time()
      }))
    },
    
    get.pathogens = function() { private$pathogens },
    add.pathogen = function(new.pathogen) {
      private$pathogens[[length(private$pathogens)+1]] <- new.pathogen
      self$set.sampling.time()
    },
    remove.pathogen = function(ex.pathogen) {
      private$pathogens[[ex.pathogen$get.name()]] <- NULL
      self$set.sampling.time()
    },
    
    is.sampled = function() {
      !all(is.na(private$sampling.time))
    }
),

private = list(
  name = NULL,
  compartment = NULL,
  source = NULL,
  transmission.time = NULL,
  sampling.time = NULL,
  pathogens = NULL,
  
  deep_clone = function(name, value) {
    if (name == 'pathogens') {
      # map deep clone to Pathogen copy() method
      lapply(value, function(pathogen) pathogen$copy(deep=TRUE))
    } 
    else {
      value
    }
  }
)
)


#' Pathogen
#' 
#' \code{Pathogen} is an R6 class for objects that represent pathogen lineages
#' that are carried by Hosts and which comprise the "inner" tree of the 
#' simulation.
#' 
#' @param name: character, unique identifier for Pathogen object
#' @param sampling.time: numeric, time that the Pathogen was sampled; left to 
#'                       NA for unsampled Pathogens
#' @param host: a reference to a Host object
#' 
#' @export
Pathogen <- R6Class(
  "Pathogen",
  public = list(
    initialize = function(name=NA, type=NA, sampling.time=NA, location=NA) {
      private$name <- name
      private$sampling.time <- sampling.time
      private$host <- host
    },
    
    parent = NULL,
    
    copy = function(deep=FALSE) {
      # see https://github.com/r-lib/R6/issues/110
      if (deep) {
        parent <- private$host
        private$host <- NULL  # temporarily erase before cloning!
      }
      cloned <- self$clone(deep)
      if (deep) {
        private$location <- parent  # restore original reference
      }
      
      cloned
    },
    
    get.name = function() { private$name },
    get.sampling.time = function() { private$sampling.time },
    
    get.host = function() { private$host },
    set.host = function(host) {
      private$host <- host
    },
    
    set.host.by.name = function(hostList, new.hostName) {
      # FIXME: unused, deprecated?
      new.hostObj <- hostList[[ 
        which(sapply(hostList, function(x) {x$get.name()}) == new.hostName) 
      ]]
      private$host <- new.hostObj
    }
  ),
  private = list(
    name = NULL,
    sampling.time = NULL,
    host = NULL
  )
)

