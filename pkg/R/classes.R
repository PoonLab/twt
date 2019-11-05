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
#' @examples 
#' 
#' # load CompartmentTypes from a YAML object
#' path <- system.file('extdata', 'SI.yaml', package='twt')
#' settings <- yaml.load_file(path)
#' mod <- MODEL$new(settings)
#' mod$get.types()
#' 
#' # manually specify a CompartmentType object (usually done by YAML)
#' host <- CompartmentType$new(name='host', branching.rates=c(host=0.1), bottleneck.size=1, coalescent.rate=1, wait.time.distr='rexp(1,1)')
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
    
    # FIXME: I'd rather use the actual Compartment object as a single argument
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


