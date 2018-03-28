# calculate migration rates across all compartments at a givent time
calc.migration.rates <- function(model) {
  # retrieve all lineages that are currently extant
  extant.lineages <- model$get.extant_lineages()
  possibleSourceTypes <- list() # list of different types of source compartment that each lineage could possibly have migrated from 
  
  # for each lineage
  sapply(extant.lineages, function(x) {
    recipientType <- x$get.location()$get.type()$get.name()
    recipientRates - sapply(model$get.types(), function(a) {
      if (a$get.migration.rate(recipientType) == 0) { NULL }
    })
  })
  # store all possible compartments that lineage can be migrated from 
  # retrieve the type of each possible compartment
  # store the different possibilities of source -> recipient pair types migration rates
  # eliminate those with a migration rate of 0
  # for each resulting migration rate, multiply the rate by the number of source -> recipient pair types
  # add up all of the individual rates and quantities to give the total migration rate
  # overall migration rate is the rate than any migration event happens across all existing compartments
  
}