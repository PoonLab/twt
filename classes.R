library(R6)


# CompartmentType
#  Defines a type of compartment such as a class of host
#  individual (risk group) with a specific transmission rate
CompartmentType  <- R6Class("CompartmentType",
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
      self$transmission.rates <- transmission.rates               # named vector of transmission rates corresponding to different Compartment objects
      self$coalescent.rate <- coalescent.rate
      self$bottleneck.size <- bottleneck.size
      self$migration.rates <- migration.rates                     # named vector of migration rates of different Compartments
    },
    
    get.bottleneck.size = function() {
      self$bottleneck.size
    },
    
    get.coalescent.rate = function() {
      self$coalescent.rate
    },
    
    get.name = function() {
      self$name
    },
    
    get.transmission.rate = function(name.type) {
      rate <- self$transmission.rates[[name.type]]
      rate
    },
    
    get.migration.rate = function(name.type) {
      rate <- self$migration.rates[[name.type]]
      rate
    },
    
    set.bottleneck.size = function(new.size) {
      self$bottleneck.size <- new.size
    }, 
    
    set.coalescent.rate = function(new.rate) {
      self$coalescent.rate <- new.rate
    },
    
    set.name = function(new.name) {
      self$name <- new.name
    },
    
    set.transmission.rate = function(name.type, new.rate) {
      self$transmission.rates[[name.type]] <- new.rate
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
    type = NULL,          # reference to CompartmentType object
    source = NULL,
    inf.time = NULL,
    sampling.time = NULL,
    
    initialize = function(type=NA, source=NA, inf.time=NA, sampling.time=NA) {
      self$type <- type
      self$source <- source
      self$inf.time <- inf.time
      self$sampling.time <- sampling.time
    },
    
    get.type = function() {
      self$type
    },
    
    get.source = function() {
      self$source
    },
    
    get.inf.time = function() {
      self$inf.time
    },
    
    get.sampling.time = function() {
      self$sampling.time
    },
    
    set.type = function(new.type) {
      self$type <- new.type
    },
    
    set.source = function(new.source) {
      self$source <- new.source
    }
  
  ),
  private = list()
)




# Lineage
Lineage <- R6Class("Lineage",
  public = list(
    type = NULL,            # potential reference to LineageType object
    sampling.time = NULL,
    location = NULL,
    initialize = function(type=NA, sampling.time=NA, location=NA) {
      self$type <- type
      self$sampling.time <- sampling.time
      self$location <- location
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
    
    set.location = function(new.location) {
      self$location <- new.location
    }
    
  ),
  private = list()
)


