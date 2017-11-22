library(R6)

# TODO: how are we going to handle varying probabilities of 
# a lineage being transmitted from one type to another?

# CompartmentType
#  Defines a type of compartment such as a class of host
#  individual (risk group) with specific transmission rate.
CompartmentType <- R6Class("CompartmentType",
  public = list(
    name = NULL,
    transmission.rate = NULL,
    coalescent.rate = NULL,
    bottleneck.size = NULL,
    migration.rate = NULL,
    initialize = function(name=NA, transmission.rate=NA,
                          coalescent.rate=NA, bottleneck.size=NA, 
                          migration.rate=NA) {
      self$name <- name
      self$transmission.rate  <- transmission.rate
      self$coalescent.rate <- coalescent.rate
      self$bottleneck.size <- bottleneck.size
      self$migration.rate <- migration.rate
    }
  ),
  private = list()
)  

Compartment <- R6Class("Compartment",
  public = list(
    type = NULL,  # reference to CompartmentType object
    lineages = NULL,
    initialize = function(type=NA, lineages=list()) {
      self.type <- type
      self.lineages <- lineages
    }
  )
)

Lineage <- R6Class("Lineage",
  public = list(
    sampling.time = NULL,
    parent = NULL
  ),
  private = list()
)


