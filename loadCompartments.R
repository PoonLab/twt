## parsing the YAML file
# takes in a given tree 
# or takes in params to be able to construct a tree

require(yaml)

load.compartments <- function(file) {
  settings <- yaml.load_file(file)
  
  # parse settings and load into CompartmentType and Compartment objects
  
}






#--------------------------------------------------------------------------------------------

# retrieve a transmission rate for a given Compartment object
get.transmission.rate <- function(CompartmentType, Compartment) {
  if (!is.element(Compartment, names(CompartmentType$transmission.rates))) {
    stop('Compartment "', Compartment, '" is not found')
  }
  index <- which(names(CompartmentType$transmission.rates) == Compartment)
  return(CompartmentType$transmission.rates[index])
}

# set a new transmission rate for a given Compartment object 
set.transmission.rate <- function(CompartmentType, Compartment, new.transmission.rate) {
  if (!is.element(Compartment, names(CompartmentType$transmission.rates))) {
    stop('Compartment "', Compartment, '" is not found')
  }
  index <- which(names(CompartmentType$transmission.rates) == Compartment)
  CompartmentType$transmission.rates[index] <- new.transmission.rate
}

# retrieve a migration rate for a given Compartment object
get.migration.rate <- function(CompartmentType, Compartment) {
  if (!is.element(Compartment, names(CompartmentType$migration.rates))) {
    stop('Compartment "', Compartment, '" is not found')
  }
  index <- which(names(CompartmentType$migration.rates) == Compartment)
  return(CompartmentType$migration.rates[index])
}

# set a new migration rate for a given Compartment object
set.migration.rate <- function(CompartmentType, Compartment, new.migration.rate) {
  if(!is.element(Compartment, names(CompartmentType$migration.rates))) {
    stop('Compartment "', Compartment, '" is not found')
  }
  index <- which(names(CompartmentType$migration.rates) == Compartment)
  CompartmentType$migration.rates[index] <- new.migration.rate
}