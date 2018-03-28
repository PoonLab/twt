# calculate migration rates across all compartments at a givent time
calc.migration.rates <- function(model) {
  # calculate migration rates across all compartments at a given time
  # @param model = MODEL object
  # @result = total rate of ANY mgiration event occurring from the list of extant lineages
  # retrieve all lineages that are currently extant
  extant.lineages <- model$get.extant_lineages()
  possibleSourceTypes <- list() # list of different types of source compartment that each lineage could possibly have migrated from 
  
  # for each lineage
  for (x in extant.lineages) {
    recipientType <- x$get.location()$get.type()$get.name()
    recipientRates <- sapply(model$get.types(), function(a) {
      if (length(a$get.migration.rates()) == 0) {   # TODO: boot this set out completely (every rate will be zero, no such migration from this source)
        NULL
      } else {
        if (a$get.migration.rate(recipientType) == 0) { NULL }
        else { a$get.migration.rate(recipientType) }
      } 
    })
    # eliminate those with a migration rate of 0
    recipientRates[sapply(recipientRates, is.null)] <- NULL    # remove source -> recipient pairs w/ migration rates of 0
    
    if (length(recipientRates) == 0) {
      # means that this lineage will never be a recipient in a migration event
      next
    } else {
      # retrieve the type of each possible source compartment
      # store all possible compartments that a lineage could have migrated from 
      possibleSourceTypes <- c(possibleSourceTypes, list(names(recipientRates)))
      names(possibleSourcetypes)[[length(possibleSourceTypes)]] <- recipientType
    }
  }
  
  # store the different possibilities of source -> recipient pair types migration rates
  indiv.rates.n.quantities <- sapply(model$get.types(), function(b) {
    sourceType <- b$get.name()
    recipientTypes <- names(b$get.migration.rates())
    sapply(recipientTypes, function(c) {
      # for each resulting migration rate, multiply the rate by the number of source -> recipient pair types
      # rate = ?
      pairRate <- b$get.migration.rate(c)
      
      qualified.r <- which(names(possibleSourceTypes) == c)
      if (length(qualified.r) == 0) {
        nPairs <- 0
      } else {
        qualified.sr <- which(sapply(qualified.r, function(d) {
          sourcetype %in% possibleSourceTypes[d]
        }))
        if (length(qualified.sr) == 0) {
          nPairs <- 0
        } else {
          nPairs <- length(qualified.sr)                # the number of pairs that have this migration pair type
        }
      }
    })
    pairRate * nPairs
  })
  
  # total rate of ANY migration event occurring is the weighted sum of these rates in the dictionary
  sum(indiv.rates.n.quantities)
  
}