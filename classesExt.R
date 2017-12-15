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
      self$transmission.rates <- transmission.rates
      self$coalescent.rate <- coalescent.rate
      self$bottleneck.size <- bottleneck.size
      self$migration.rates <- migration.rates
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
    type = NULL,
    source = NULL,
    inf.time = NULL,
    sampling.time = NULL,
    
    initialize = function(type=NA, source=NA, ing.time=NA, sampling.time=NA) {
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
    }
  
  ),
  private = list()
)




# Lineage
Lineage <- R6Class("Lineage",
  public = list(
    type = NULL,
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





# Bring them all together
Mastermind <- R6Class("Mastermind",
  public = list(
    settings = NULL,
    types = NULL,
    compartments = NULL,
    lineages = NULL,
    
    initialize = function(settings=NA) {
      #self$types <- settings
      #self$compartments <- settings
      #self$lineages <- settings
      self$load.types(settings)
      self$load.compartments(settings)
      self$load.lineages(settings)
    },
    
    
    load.types = function(settings) {
      types <- sapply(names(settings$CompartmentType), function(x) {
                 params <- settings$CompartmentType[[x]]
                 x <- CompartmentType$new(name = x,
                                          transmission.rates = params$transmission.rates,
                                          migration.rates = params$migration.rates,
                                          coalescent.rate = params$coalescent.rate,
                                          bottleneck.size = params$bottleneck.size
                                          )
      })
      self$types <- types
    },
    
    
    load.compartments = function(settings) {
      compartments <- sapply(names(settings$Compartments), function(x) {
        compartX <- list()
        params <- settings$Compartments[[x]]
        nIndiv <- params$pop.size
        for(obj in 1:nIndiv) {
          x <- Compartment$new(type = params$type,
                               source = params$source,
                               inf.time = params$inf.time,
                               sampling.time = params$sampling.time
                               )
          compartX[[obj]] <- x
        }
        compartX
      })
      self$compartments <- compartments
    },
    
    
    load.lineages = function(settings) {
      lineages <- sapply(names(settings$Lineages), function(x) {
        lineageX <- list()
        params <- settings$Lineages[[x]]
        nIndiv <- params$pop.size
        for (obj in 1:nIndiv) {
          x <- Lineage$new(type = params$type,
                           sampling.time = params$sampling.time,
                           location = params$location
                           )
          lineageX[[obj]] <- x
        }
        lineageX
      })
      self$lineages <- lineages
    },
    
    
    get.types = function() {
      self$types
    },
    
    get.compartments = function() {
      self$compartments
    },
    
    get.lineages = function() {
      self$lineages
    }
    
  ),
  private = list()
)

