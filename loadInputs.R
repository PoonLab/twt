## test script
require(R6)
require(yaml)
setwd('~/git/treeswithintrees')
settings <- yaml.load_file('example2.yaml')
test <- NestedCoalescent$new(settings)
e <- EventLogger$new()
tips.n.heights <- init.fixed.samplings(test)
init.fixed.transmissions(test, e)    # applies only to example1.yaml for now
e <- generate.transmission.events(test, e)


# Load all of the different objects into one larger class
NestedCoalescent <- R6Class("NestedCoalescent",
  public = list(
    
    settings = NULL,
    types = NULL,
    unsampled.hosts = NULL,
    compartments = NULL,
    lineages = NULL,

    extant_lineages = NULL,
    extant_comps = NULL,
    
    locations = NULL,         
    choices = NULL,
    
    initialize = function(settings=NA) {
      private$load.types(settings)
      private$load.unsampled.hosts()
      private$load.compartments(settings)
      private$set.sources()
      private$load.lineages(settings)
      
      private$init.extant.lineages()
      self$extant_comps <- self$compartments
    },
    
    get.types = function() {self$types},
    get.unsampled.hosts = function() {self$unsampled.hosts},
    get.compartments = function() {self$compartments},
    get.lineages = function() {self$lineages},
    get.extant_lineages = function() {self$extant_lineages},
    get.extant_comps = function() {self$extant_comps},
    
    get.pairs = function() {
      # extract all pairs of pathogen lineages that may coalesce
      self$get.locations()
      self$choices <- list()
      for (host in self$locations) {
        if (length(host) > 1) {
          for (pair in combn(host,2)) {
            self$choices[[paste(pair[1], pair[2], sep=',')]] <- host
          }
        }
      }
    },
    
    get.locations = function() {
      # collect host locations of all extant pathogen lineages into dict of host1: [path1, path2, path3, ...]
      self$locations <- list()      # reset list
      for (node in self$extant_lineages) {
        if (node$get.location()$get.name() %in% names(self$locations)) {
          self$locations[[node$get.location()$get.name()]] <- list()
        }
        self$locations[[node$get.location()$get.name()]] <- c(self$locations[[node$get.location()$get.name()]], node)
      }
    }
    
  ),
  
  
  private = list(
    
    ## function creates CompartmentType objects
    ## within each CompartmentType, there are distinct compartments with: 
        # individual transmission & migration rates
        # unsampled host & susceptible populations 
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
      self$types <- unlist(types)
    },
    
    
    ## function creates "blank" Compartment objects for Unsampled Hosts (US) 
    ## stored in lists for each section within a CompartmentType object
    load.unsampled.hosts = function() {
      us.hosts <- sapply(self$types, function(x) {
        types.unsampled <- names(x$no.unsampled)           # FIXME: accessing a private variable here; maybe add another public method in CompartmentType instead
        
        lapply(types.unsampled, function(y) {
          nBlanks <- x$get.no.unsampled(y)
          sapply(1:nBlanks, function(blank) {
            Compartment$new(name = paste0('US_', x$get.name(), '_', blank),
                            type = x)
          })
        })
        
      })
      self$unsampled.hosts <- unlist(us.hosts)
    },
    
    
    ## function creates Compartment objects
    ## `type` attr points directly back to a CompartmentType object, and `name` attr is a unique identifier
    load.compartments = function(settings) {
      compartments <- sapply(names(settings$Compartments), function(comp) {
        params <- settings$Compartments[[comp]]
        
        if (params$type %in% names(settings$CompartmentType)) {
          searchTypes <- which(names(settings$CompartmentType) == params$type)
          typeObj <- self$types[[ searchTypes ]]                        # pointer to CompartmentType object
        } else {
          stop(params$type, ' of Compartment ', comp, ' is not a specified Compartment Type object')
        }
        
        nIndiv <- params$pop.size
        
        sapply(1:nIndiv, function(obj) {
          Compartment$new(name = paste0(comp,'_', obj),            # unique identifier
                          type = typeObj,
                          source = params$source,        
                          inf.time = params$inf.time,
                          sampling.time = params$sampling.time
          )
        })
        
      })
      self$compartments <- unlist(compartments)
    },
    
    ## re-iterates over generated Compartment objects and populates `source` attr with R6 objects
    ## sets 'pointers' to other Compartment objects after generated w/ private$load.compartments()
    set.sources = function() {
      compNames <- sapply(self$compartments, function(n){n$get.name()})
      compartments <- sapply(self$compartments, function(x) {
        if (paste0(x$get.source(),'_1') %in% compNames) {                         # FIXME: arbitrary assignment of source
          searchComps <- which(compNames == paste0(x$get.source(),'_1'))          # FIXME: arbitrarily assigning source to first object in Compartment with $name == x$get.source()
          sourceObj <- self$compartments[[ searchComps ]]    
          x$set.source(sourceObj)
        } # TODO: else statement { if source is 'undefined' or not in the list, must be assigned to an unsampled host (US) }
        x
      })
      self$compartments <- compartments
    },
    
    
    ## function creates Lineage objects
    ## `location` attr points directly to a Compartment object, and `name` attr is unique identifier
    ## identifiers create unique Lineages for each Compartment, 
    ## but Compartment A_1 could have Lineage cell_1 and Compartment B_1 also have a separate Lineage cell_1
    load.lineages = function(settings) {
      lineages <- sapply(names(settings$Lineages), function(label) {
        params <- settings$Lineages[[label]]
        
        if (params$location %in% names(settings$Compartments)) {
          searchComps <- which(names(settings$Compartments) == params$location)
          nlocationObj <- settings$Compartments[[ searchComps ]]$pop.size
        } else {
          stop(params$location, ' of Lineage ', label, ' is not a specified Compartment object')
        }
        
        nIndiv <- params$pop.size
        
        if (is.numeric(params$sampling.time)) {
          sampleTimes <- rep.int(params$sampling.time, times=nIndiv)
        } else {
          vec <- unlist(strsplit(params$sampling.time, split='[`(`|,|`)`]'))
          sampleTimes <- as.double(vec[nzchar(x=vec)])
          if (any(is.na(sampleTimes))) {
              stop('\nSampling times must all be specified as type `numeric` or type `double`. ', 
                   'Sampling time "', vec, '" is of type `', typeof(vec), '`.\n')
          }
          if (length(sampleTimes) != nIndiv) {
            stop('attribute `sampling.time` of Lineage ', label, ' does not match pop.size specified for respective Lineage.')
          }
        }
        
        sapply(1:nlocationObj, function(compNum) {
          sapply(1:nIndiv, function(obj) {
            # set 'pointer' to Compartment object for location
            searchComps <- sapply(self$compartments, function(y){which(y$get.name() == paste0(params$location, '_', compNum))})
            locationObj <- self$compartments[[ which( searchComps == 1) ]]
            Lineage$new(name = paste0(label,'_',obj),                          # unique identifier
                        type = params$type,
                        sampling.time = sampleTimes[obj],
                        location = locationObj
            )
          })
        })
      
      })
      self$lineages <- unlist(lineages)
    },
    
    
    # intializes list of Lineages with sampling.time t=0
    init.extant.lineages = function() {
      res <- sapply(self$lineages, function(b){
        if (b$get.sampling.time() == 0) {b}
      })
      self$extant_lineages <- unlist(res)
    }
    
  )
)

