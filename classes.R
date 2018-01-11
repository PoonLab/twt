library(R6)


# CompartmentType
#  Defines a type of compartment such as a class of host
#  individual (risk group) with a specific transmission rate
CompartmentType  <- R6Class("CompartmentType",
  public = list(
    name = NULL,
    no.unsampled = NULL,
    no.susceptible = NULL,
    transmission.rates = NULL,
    coalescent.rate = NULL,
    bottleneck.size = NULL,
    migration.rates = NULL,
    
    initialize = function(name=NA, no.unsampled = NA,
                          no.susceptible=NA, transmission.rates=NA,
                          coalescent.rate=NA, bottleneck.size=NA,
                          migration.rates=NA) {
      self$name <- name
      self$no.unsampled = no.unsampled
      self$no.susceptible = no.susceptible
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
    
    get.no.unsampled = function(name.type) {
      num <- self$no.unsampled[[name.type]]
      num
    },
    
    get.no.susceptible = function(name.type) {
      num <- self$no.susceptible[[name.type]]
      num
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
    
    set.no.unsampled = function(name.type, new.no) {
      self$no.unsampled[[name.type]] <- new.no
    },
    
    set.no.susceptible = function(name.type, new.no) {
      self$no.susceptible[[name.type]] <- new.no
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
    name = NULL,
    type = NULL,          # reference to CompartmentType object
    source = NULL,
    inf.time = NULL,
    sampling.time = NULL,
    
    initialize = function(name=NA, type=NA, source=NA, inf.time=NA, sampling.time=NA) {
      self$name <- name
      self$type <- type
      self$source <- source
      self$inf.time <- inf.time
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
    
    set.location = function(new.location) {
      self$location <- new.location
    }
    
  ),
  private = list()
)




# Event Logger (tracks migration, transmission, coalescent, and bottleneck events)
EventLogger <- R6Class("EventLogger",
  public = list(
    migrations = NULL,
    transmissions = NULL,
    coalescences = NULL,
    bottlenecks = NULL,
    initialize = function(migrations=NA, transmissions=NA, coalescences=NA, bottlenecks=NA) {
      self$migrations <- migrations
      self$transmissions <- transmissions
      self$coalescences <- coalescences
      self$bottlenecks <- bottlenecks
    },
    
    add.migration = function(time, lineage1, compartment1, compartment2) {
      new.event <- list(time=time, L1=lineage1, C1=compartment1, C2=compartment2)
      self$migrations <- c(self$migrations, new.event)
    },
    
    add.transmission = function(time, lineage1, compartment1, compartment2) {
     new.event <- list(time=time, L1=lineage1, C1=compartment1, C2=compartment2)
     self$transmissions <- c(self$transmissions, new.event)
    },
    
    add.coalescence = function(time, lineage1, lineage2, compartment1) {
      new.event <- list(time=time, L1=lineage1, L2=lineage2, C1=compartment1)
      self$coalescences <- c(self$coalescences, new.event)
    },
    
    add.bottleneck = function(time, lineage1, lineage2, compartment1) {
      new.event <- list(time=time, L1=lineage1, L2=lineage2, C1=compartment1)
      self$bottlenecks <- c(self$bottlenecks, new.event)
    }
    
  ),
  private = list()
)
