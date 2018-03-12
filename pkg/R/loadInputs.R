require(R6)

# Load all of the different objects into one larger class
MODEL <- R6Class("MODEL",
  public = list(
    initialize = function(settings=NA) {
      private$load.types(settings)
      private$load.unsampled.hosts()
      private$load.compartments(settings)
      private$set.sources()
      private$load.lineages(settings)
      
      private$init.extant.lineages()
      private$init.extant.comps()
      private$init.non.extant.comps()
      
      private$init.locations()
      private$init.pairs()
    },
    
    get.types = function() {private$types},
    get.unsampled.hosts = function() {private$unsampled.hosts},
    get.compartments = function() {private$compartments},
    get.lineages = function() {private$lineages},
    get.extant_lineages = function() {private$extant_lineages},
    get.extant_comps = function() {private$extant_comps},
    get.non_extant_comps = function() {private$non_extant_comps},
    
    get.leaves.names = function(e) {
      # return a vector of Compartment object names that are terminal nodes (only a recipient)
      # @param e = EventLogger object
      t_events <- e$get.events('transmission', cumulative=F)
      setdiff(unlist(t_events$compartment1), unlist(t_events$compartment2))
    },
    
    get.nonterminals = function(e) {
      # return an iterator over all names of internal nodes of the transmission tree
      # @param e = EventLogger object
      t_events <- e$get.events('transmission', cumulative=F)
      intersect(unlist(t_events$compartment1), unlist(t_events$compartment2))
    },
    
    get.node.heights = function() {
      # calculate node heights for all nodes of the tree
      # annotate nodes with heights in place
    },
    
    get.pairs = function() {
      # function extracts and returns all the current pairs of pathogen lineages that may coalesce
      private$choices
    },
    
    add.pair = function(L1, L2, host) {
      # when a Lineage is moved from one compartment to another (transmission or migration)
      # when a Lineage is sampled
      # can also be used to update the location of a pair
      pair <- sort(c(L1, L2))
      private$choices[[paste(pair[1], pair[2], sep=',')]] <- host
    },
    
    remove.pair = function(L1, L2) {
      # when a coalescence occurs
      # when Lineages reach a tranmission bottleneck, forcing coalescence
      pair <- sort(c(L1, L2))
      private$choices[[paste(pair[1], pair[2], sep=',')]] <- NULL
    }
    
  ),
  
  
  private = list(
    types = NULL,
    unsampled.hosts = NULL,
    compartments = NULL,
    lineages = NULL,
    
    extant_lineages = NULL,
    extant_comps = NULL,
    non_extant_comps = NULL, 
    
    locations = NULL,         
    choices = NULL,
    

    load.types = function(settings) {
      ## function creates CompartmentType objects
      ## within each CompartmentType, there are distinct compartments with: 
      # individual transmission & migration rates
      # unsampled host & susceptible populations 
      types <- sapply(names(settings$CompartmentTypes), function(x) {
        params <- settings$CompartmentTypes[[x]]
        x <- CompartmentType$new(name = x,
                                 unsampled = params$unsampled,
                                 susceptible = params$susceptible,
                                 branching.rates = eval(parse(text=paste('list', params$branching.rates))),
                                 migration.rates = eval(parse(text=paste('list', params$migration.rates))),
                                 bottleneck.size = params$bottleneck.size,
                                 popn.growth.dynamics = private$init.popn.growth.dynamics(params$popn.growth.dynamics)
        )
      })
      private$types <- unlist(types)
    },
    
    
    
    load.unsampled.hosts = function() {
      ## function creates "blank" Compartment objects for Unsampled Hosts (US) 
      ## stored in lists for each section within a CompartmentType object
      us.hosts <- sapply(private$types, function(x) {
          nBlanks <- x$get.unsampled()
          sapply(1:nBlanks, function(blank) {
            Compartment$new(name = paste0('US_', x$get.name(), '_', blank),
                            type = x)
          })
      })
      private$unsampled.hosts <- unlist(us.hosts)
    },
    
    
   
    load.compartments = function(settings) {
      ## function creates Compartment objects
      ## `type` attr points directly back to a CompartmentType object, and `name` attr is a unique identifier
      compartments <- sapply(names(settings$Compartments), function(comp) {
        params <- settings$Compartments[[comp]]
        
        if (params$type %in% names(settings$CompartmentType)) {
          searchTypes <- which(names(settings$CompartmentType) == params$type)
          typeObj <- private$types[[ searchTypes ]]                  # pointer to CompartmentType object
        } else {
          stop(params$type, ' of Compartment ', comp, ' is not a specified Compartment Type object')
        }
        
        nIndiv <- params$replicates
        
        sapply(1:nIndiv, function(obj) {
          Compartment$new(name = paste0(comp,'_', obj),            # unique identifier
                          type = typeObj,
                          source = params$source,        
                          branching.time = params$branching.time,
                          sampling.time = params$sampling.time
          )
        })
        
      })
      private$compartments <- unlist(compartments)
    },
    
    
    set.sources = function() {
      ## re-iterates over generated Compartment objects and populates `source` attr with R6 objects
      ## sets 'pointers' to other Compartment objects after generated w/ private$load.compartments()
      compNames <- sapply(private$compartments, function(n){n$get.name()})
      compartments <- sapply(private$compartments, function(x) {
        if (paste0(x$get.source(),'_1') %in% compNames) {                         # FIXME: arbitrary assignment of source
          searchComps <- which(compNames == paste0(x$get.source(),'_1'))          # FIXME: arbitrarily assigning source to first object in Compartment with $name == x$get.source()
          sourceObj <- private$compartments[[ searchComps ]]    
          x$set.source(sourceObj)
        } # TODO: else statement { if source is 'undefined' or not in the list, must be assigned to an unsampled host (US) }
        x
      })
      private$compartments <- compartments
    },
    
    

    load.lineages = function(settings) {
      ## function creates Lineage objects
      ## `location` attr points directly to a Compartment object, and `name` attr is unique identifier
      ## identifiers create unique Lineages for each Compartment, 
      ## but Compartment A_1 could have Lineage cell_1 and Compartment B_1 also have a separate Lineage cell_1
      lineages <- sapply(names(settings$Lineages), function(label) {
        params <- settings$Lineages[[label]]
        
        if (params$location %in% names(settings$Compartments)) {
          searchComps <- which(names(settings$Compartments) == params$location)
          nlocationObj <- settings$Compartments[[ searchComps ]]$replicates
        } else {
          stop(params$location, ' of Lineage ', label, ' is not a specified Compartment object')
        }
        
        nIndiv <- params$replicates
        
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
            stop('attribute `sampling.time` of Lineage ', label, ' does not match number of replicates specified for respective Lineage.')
          }
        }
        
        sapply(1:nlocationObj, function(compNum) {
          sapply(1:nIndiv, function(obj) {
            # set 'pointer' to Compartment object for location
            searchComps <- sapply(private$compartments, function(y){which(y$get.name() == paste0(params$location, '_', compNum))})
            locationObj <- private$compartments[[ which( searchComps == 1) ]]
            Lineage$new(name = paste0(locationObj$get.name(),':',label,'_',obj),                          # unique identifier
                        type = params$type,
                        sampling.time = sampleTimes[obj],
                        location = locationObj
            )
          })
        })
      
      })
      private$lineages <- unlist(lineages)
    },
    
    
   
    init.extant.lineages = function() {
      # intializes list of Lineages with sampling.time t=0
      res <- sapply(private$lineages, function(b){
        if (b$get.sampling.time() == 0) {b}
      })
      private$extant_lineages <- unlist(res)
    },
    
    
    
    init.extant.comps = function() {
      # initializes list of Compartments containing Lineages with sampling.time t=0
      res <- sapply(private$extant_lineages, function(b){
        b$get.location()
      })
      private$extant_comps <- unique(res)
    },
    
    
    
    init.non.extant.comps = function() {
      # initializes list of Compartments containing Lineages with sampling.time t!=0
      extant.names <- sapply(private$extant_comps, function(a){a$get.name()})
      res <- sapply(private$compartments, function(b) {
        if (b$get.name() %in% extant.names == F) {b}
      })
      private$non_extant_comps <- unlist(res)
    },
    
    
    
    init.popn.growth.dynamics = function(pieces) {
      # @param pieces, list of linear pieces of a given CompartmentType
      # @return matrix, where each row represents a linear piece and comprises the following columns
        # start time
        # end time
        # intercept
        # slope
      mat <- matrix(nrow=length(pieces), ncol=4)
      res <- t(sapply(seq_along(pieces), function(x) {
        sapply(seq_along(pieces[[x]]), function(y){
          if (names(pieces[[x]][[y]]) == 'end') {
            if (is.character(pieces[[x]][[y]])) {         # deals with single case where end time infinity is of mode character
              entry <- NA
            } else {
              entry <- pieces[[x]][[y]]
            }
          } else {
            entry <- pieces[[x]][[y]]
          }
          mat[x,y] <- entry
        })
      }))
      
      # checks for the following:
      
      # only one piece must have a start time of zero
      if ( length(which(res[,'start'] == 0)) != 1 ) {
        stop ('One and only one linear piece in the population growth dynamics functions is allowed to have a start time of zero (forward in time).')
      }
      # only one piece must have an end time of infinity (undefined)
      if ( length( which(sapply(res[,'end'], function(x){is.character(x)})) ) != 1 ) {
        stop ('One and only one linear piece in the population growth dynamics functions is allowed to have an undefined end time (forward in time).')
      }
      # every other start time must have an equal end time and vice versa ( no gaps in time for piecewise growth function)
      if ( length(intersect(res[,'start'], res[,'end'])) != (nrow(res)-1) ) {
        stop ('There are time gaps and/or overlaps within the population growth dynamics functions.')
      }
      
      res
    },
    
    
    
    init.locations = function() {
      # helper function for private$init.pairs()
      # collect host locations of all extant pathogen lineages into dict of host1: [path1, path2, path3, ...]
      private$locations <- list()      # reset list
      for (node in private$extant_lineages) {
        if (node$get.location()$get.name() %in% names(private$locations) == F) {
          private$locations[[node$get.location()$get.name()]] <- list()
        }
        private$locations[[node$get.location()$get.name()]] <- c(private$locations[[node$get.location()$get.name()]], node$get.name())
      }
    },
    
    
    
    init.pairs = function() {
      # extract all pairs of pathogen lineages that may coalesce
      private$choices <- list()
      for (hostNum in seq_along(private$locations)) {
        host <- private$locations[[hostNum]]
        if (length(host) > 1) {
          combns <- combn(host,2)
          for (column in 1:ncol(combns)) {
            pair <- sort(unlist(combns[,column]))
            private$choices[[paste(pair[1], pair[2], sep=',')]] <- names(private$locations)[[hostNum]]
          }
        }
      }
    }
    
    
  )
)

