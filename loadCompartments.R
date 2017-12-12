## parsing the YAML file
# takes in a given tree 
# or takes in params to be able to construct a tree

require(yaml)

load.compartments <- function(file) {
  settings <- yaml.load_file(file)
  
  # parse settings and load into CompartmentType and Compartment objects
  compartmentTypes <- names(settings$CompartmentType)
  compartments <- names(settings$Compartments)
  lineages <- names(settings$Lineages)
  
  init.objects <- sapply(compartmentTypes, function(x) {
    params <- settings$CompartmentType[[x]]
    x <- CompartmentType$new(name = x,
                             transmission.rates = params$transmission.rates,
                             migration.rates = params$migration.rates,
                             coalescent.rate = params$coalescent.rate,
                             bottleneck.size = params$bottleneck.size
                             )
  })
  
  
  init.compartments <- sapply(compartments, function(x) {
    compartX <- list()
    params <- settings$Compartments[[x]]
    nIndiv <- params$pop.size                                # potential for multiple generic objects w/ same params
    for (obj in 1:nIndiv) {                                  
      x <- Compartment$new(type = params$type,
                           source = params$source,
                           inf.time = params$inf.time,
                           sampling.time = params$sampling.time
                           )
      compartX[[obj]] <- x
    }
    compartX
  })
  
  
  init.lineages <- sapply(lineages, function(x) {
    lineageX <- list()
    params <- settings$Lineages[[x]]
    nIndiv <- params$pop.size
    for(obj in 1:nIndiv) {
      x <- Lineage$new(type = params$type,
                       sampling.time = params$sampling.time,
                       location = params$location
                      )
      lineageX[[obj]] <- x
    }
    lineageX
  })
  
  
  return (c(init.objects, init.compartments, init.lineages))
}




#--------------------------------------------------------------------------------------------
## Are these getters and setters needed? They can be easily accessed directly within each object

# retrieve a transmission rate for a given Compartment object
get.transmission.rate <- function(CompartmentType, Compartment) {
  if (!exists(CompartmentType)) { stop('CompartmentType "', CompartmentType, '" is not found') }     #could also check if present specifically in the Global env
  if (!is.element(Compartment, names(CompartmentType$transmission.rates))) { stop('Compartment "', Compartment, '" is not found') }
  index <- which(names(CompartmentType$transmission.rates) == Compartment)
  return(CompartmentType$transmission.rates[index])
}

# set a new transmission rate for a given Compartment object 
set.transmission.rate <- function(CompartmentType, Compartment, new.transmission.rate) {
  if (!exists(CompartmentType)) { stop('CompartmentType "', CompartmentType, '" is not found') }     
  if (!is.element(Compartment, names(CompartmentType$transmission.rates))) { stop('Compartment "', Compartment, '" is not found') }
  index <- which(names(CompartmentType$transmission.rates) == Compartment)
  CompartmentType$transmission.rates[index] <- new.transmission.rate
}

# retrieve a migration rate for a given Compartment object
get.migration.rate <- function(CompartmentType, Compartment) {
  if (!exists(CompartmentType)) { stop('CompartmentType "', CompartmentType, '" is not found') }
  if (!is.element(Compartment, names(CompartmentType$migration.rates))) { stop('Compartment "', Compartment, '" is not found') }
  index <- which(names(CompartmentType$migration.rates) == Compartment)
  return(CompartmentType$migration.rates[index])
}

# set a new migration rate for a given Compartment object
set.migration.rate <- function(CompartmentType, Compartment, new.migration.rate) {
  if (!exists(CompartmentType)) { stop('CompartmentType "', CompartmentType, '" is not found') }
  if (!is.element(Compartment, names(CompartmentType$migration.rates))) { stop('Compartment "', Compartment, '" is not found') }
  index <- which(names(CompartmentType$migration.rates) == Compartment)
  CompartmentType$migration.rates[index] <- as.numeric(new.migration.rate)
}