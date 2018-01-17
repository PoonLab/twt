## test script
require(R6)
require(yaml)
setwd('~/git/treeswithintrees')
settings <- yaml.load_file('example1.yaml')
test <- NestedCoalescent$new(settings)
e <- EventLogger$new()
tips.n.heights <- init.fixed.samplings(test)
init.fixed.transmissions(test, e)


# Load all of the different objects into one larger class
NestedCoalescent <- R6Class("NestedCoalescent",
  public = list(
    settings = NULL,
    types = NULL,
    unsampled.hosts = NULL,
    compartments = NULL,
    lineages = NULL,

    extant = NULL,
    choices = NULL,
    
    initialize = function(settings=NA) {
      private$load.types(settings)
      private$load.unsampled.hosts()
      private$load.compartments(settings)
      private$set.sources()
      private$load.lineages(settings)
      # TODO: populate extant with extant lineages
    },
    
    get.types = function() {self$types},
    get.unsampled.hosts = function() {self$unsampled.hosts},
    get.compartments = function() {self$compartments},
    get.lineages = function() {self$lineages},
    
    
    ## collect host locations of all extant lineages into named list of host1:[pathogen1, pathogen2, ...] name-value pairs
    get.locations = function() {
      locations <- list()  # reset the list
      for (node in self$extant) {    # sapply does not work here... will create list of `$1.host`, `$2.host`, etc
        # TODO: check that lineage is extant
        my.comp.type <- node$get.location()$get.type()$get.name()
        # append this lineage to the vector
        locations[[my.comp.type]] <- c(locations[[my.comp.type]], node)
      }
      locations
    },
    
    
    ## function to extract all pairs of lineages that may coalesce
    get.pairs = function() {
      locations <- self$get.locations()
      choices <- list()
      sapply(locations, function(x) {
        if (length(x) > 1) {
          pairs <- t(combn(1:length(x), 2))
          for (row in 1:nrow(pairs)) {
            pair <- pairs[row,]
            choices[[pair]] <- c(choices[[pair]], x)    # update list of pathogen pairs in same host
          } # TODO: store not in tuples, but in some other data structure
        }
      })
      self$choices <- choices
    }
    
    
  ),
  
  
  private = list(
    
    ## function creates CompartmentType objects
    ## within each CompartmentType, there are distinct compartments with individual transmission & migration rates, and unsampled host & susceptible populations 
    load.types = function(settings) {
      types <- sapply(names(settings$CompartmentType), function(x) {
        params <- settings$CompartmentType[[x]]
        x <- CompartmentType$new(name = x,
                                 no.unsampled = eval(parse(text=paste('list', params$no.unsampled))),
                                 no.susceptible = eval(parse(text=paste('list', params$no.susceptible))),
                                 transmission.rates = eval(parse(text=paste('list', params$transmission.rates))),
                                 migration.rates = eval(parse(text=paste('list', params$migration.rates))),
                                 coalescent.rate = params$coalescent.rate,
                                 bottleneck.size = params$bottleneck.size
        )
      })
      self$types <- types
    },
    
    
    ## function creates "blank" Compartment objects for Unsampled Hosts (US) 
    ## stored in lists for each section within a CompartmentType object
    load.unsampled.hosts = function() {
      us.hosts <- sapply(self$types, function(x) {
        types.unsampled <- names(x$no.unsampled)           # accessing a private variable here; maybe add another public method in CompartmentType instead
        indiv <- lapply(types.unsampled, function(y) {
          compartY <- list()
          nBlanks <- x$get.no.unsampled(y)
          for(blank in 1:nBlanks) {
            x <- Compartment$new()
            compartY[[blank]] <- x
          }
          compartY
        })
        indiv
      })
      self$unsampled.hosts <- us.hosts
    },
    
    
    ## function creates Compartment objects
    ## `type` attr points directly back to a CompartmentType object, and `name` attr is a unique identifier
    load.compartments = function(settings) {
      compartments <- sapply(names(settings$Compartments), function(comp) {
        compartX <- list()
        params <- settings$Compartments[[comp]]
        if (params$type %in% names(self$types)) {
          typeObj <- self$types[[ which(names(self$types) == params$type) ]]    # pointer to CompartmentType object
        } else {
          stop(params$type, ' of Compartment ', comp, ' is not a specified Compartment Type object')
        }
        nIndiv <- params$pop.size
        for(obj in 1:nIndiv) {
          x <- Compartment$new(name = paste0(comp,'_', obj),                     # unique identifier
                               type = typeObj,
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
    
    ## re-iterates over generated Compartment objects and populates `source` attr with R6 objects
    ## sets 'pointers' to other Compartment objects after all have been generated with function call private$load.compartments()
    set.sources = function() {
      compNames <- sapply(self$compartments, function(n){n$get.name()})
      compartments <- sapply(self$compartments, function(x) {
        if (paste0(x$get.source(),'_1') %in% compNames) {                                        # FIXME: arbitrary assignment of source
          sourceObj <- self$compartments[[ which(compNames == paste0(x$get.source(),'_1')) ]]    # FIXME: arbitrarily assigning source to first object in Compartment with $name == x$get.source()
          x$set.source(sourceObj)
        } # TODO: else statement { if source is 'undefined' or not in the list, must be assigned to an unsampled host (US) }
        x
      })
      self$compartments <- compartments
    },
    
    
    ## function creates Lineage objects
    ## `location` attr points directly to a Compartment object, and `name` attr is unique identifier
    load.lineages = function(settings) {
      lineages <- sapply(names(settings$Lineages), function(label) {
        lineageX <- list()
        params <- settings$Lineages[[label]]
        
        if (params$location %in% names(self$compartments)) {                            # set 'pointer' to Compartment object for location
          locationObj <- self$compartments[[ which( sapply(self$compartments, function(y){which(y$get.name() == paste0(params$location,'_1'))}) ==1) ]]  # FIXME: arbitrarily assigns location to first object in Compartment with $name == params$location
        } else {
          stop(params$location, ' of Lineage ', label, ' is not a specified Compartment object')
        }
        nIndiv <- params$pop.size
        if (is.character(params$sampling.time)) {
          vec <- unlist(strsplit(params$sampling.time, split='[`(`|,|`)`]'))
          sampleTimes <- as.double(vec[nzchar(x=vec)])
          if (length(sampleTimes) != nIndiv) { stop(paste('Lineage', label, 'sampling.time does not match pop.size specified for respective Lineage.'))}
        } else {
          sampleTimes <- params$sampling.time
        }
        for (obj in 1:nIndiv) {
          x <- Lineage$new(name = paste0(label,'_',obj),                                 # unique identifier
                           type = params$type,
                           sampling.time = sampleTimes[obj],
                           location = locationObj
          )
          lineageX[[obj]] <- x
        }
        lineageX
      })
      self$lineages <- lineages
    }
    
  )
)
