library(R6)


# CompartmentType
#  Defines a type of compartment such as a class of host
#  individual (risk group) with specific transmission rate.
CompartmentType <- R6Class("CompartmentType",
  public = list(
    name = NULL,
    transmission.rates = NULL,
    coalescent.rate = NULL,
    bottleneck.size = NULL,
    migration.rates = NULL,
    initialize = function(name=NA, transmission.rates=NA,
                          coalescent.rate=NA, bottleneck.size=NA, 
                          migration.rates=NA) {
      self$name <- name
      self$transmission.rates  <- transmission.rates   # named vector of transmission rates corresponding to different Compartment objects
      self$coalescent.rate <- coalescent.rate
      self$bottleneck.size <- bottleneck.size
      self$migration.rates <- migration.rates          # named vector of migrations rates of different Compartments
    }
  ),
  private = list()
)  

Compartment <- R6Class("Compartment",
  public = list(
    type = NULL,  # reference to CompartmentType object
    source = NULL,
    inf.time = NULL,
    sampling.time = NULL,
    initialize = function(type=NA, source=NA, inf.time=NA, sampling.time=NA) {
      self$type <- type
      self$source <- source
      self$inf.time <- inf.time
      self$sampling.time <- sampling.time
    }
  )
)

Lineage <- R6Class("Lineage",
  public = list(
    type = NULL,           # potential reference to LineageType object
    sampling.time = NULL,
    location = NULL,
    initialize = function(type=NA, sampling.time=NA, location=NA) {
      self$type <- type
      self$sampling.time <- sampling.time
      self$location <- location
    }
  ),
  private = list()
)


