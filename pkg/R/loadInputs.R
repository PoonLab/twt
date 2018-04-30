# Load all of the different objects into one larger class
MODEL <- R6Class("MODEL",
  public = list(
    initialize = function(settings=NA) {
      private$types <- private$load.types(settings)
      private$unsampled.hosts <- private$load.unsampled.hosts()
      private$compartments <- private$load.compartments(settings)
      private$compartments <- private$set.sources()
      private$lineages <- private$load.lineages(settings)
      
      private$extant_lineages <- private$get.extant.lineages(0)
      private$locations <- private$init.locations()
      private$choices <- private$init.pairs()
    },
    
    get.types = function() {private$types},
    get.unsampled.hosts = function() {private$unsampled.hosts},
    get.compartments = function() {private$compartments},
    get.lineages = function() {private$lineages},
    
    get.extant_lineages = function(time) {
      private$extant_lineages <- private$get.extant.lineages(time)
      private$extant_lineages
    },
    
    get.names = function(listR6obj) {
      # @param listR6obj = list of R6 objects of class CompartmentType, Compartment, or Lineage
      unname(sapply(listR6obj, function(x){x$get.name()}))
    },
    
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
    }, 
    
    add.lineage = function(ancestral.lineage) {
      # at a coalescent event, an ancestral lineage must be created
      private$lineages[[length(private$lineages)+1]] <- ancestral.lineage
    },
    
    remove.lineage = function(lineage) {
      # at a coalescent event, lineages that coalesce must be removed
      lin.index <- which(sapply(private$lineages, function(x){x$get.name() == lineage$get.name()}))
      private$lineages[[lin.index]] <- NULL
    }
    
  ),
  
  
  private = list(
    types = NULL,
    unsampled.hosts = NULL,
    compartments = NULL,
    lineages = NULL,
    
    extant_lineages = NULL,
    
    locations = NULL,         
    choices = NULL,
    

    load.types = function(settings) {
      ## function creates CompartmentType objects
      ## within each CompartmentType, there are distinct compartments with: 
      # individual transmission & migration rates
      # unsampled host & susceptible populations 
      unlist(sapply(names(settings$CompartmentTypes), function(x) {
        params <- settings$CompartmentTypes[[x]]
        x <- CompartmentType$new(name = x,
                                 unsampled = params$unsampled,
                                 susceptible = params$susceptible,
                                 branching.rates = eval(parse(text=paste('list', params$branching.rates))),
                                 migration.rates = eval(parse(text=paste('list', params$migration.rates))),
                                 bottleneck.size = params$bottleneck.size,
                                 popn.growth.dynamics = private$init.popn.growth.dynamics(params$popn.growth.dynamics)
        )
      }))
    },
    
    
    
    load.unsampled.hosts = function() {
      ## function creates "blank" Compartment objects for Unsampled Hosts (US) 
      ## stored in lists for each section within a CompartmentType object
      unlist(sapply(private$types, function(x) {
          nBlanks <- x$get.unsampled()
          sapply(1:nBlanks, function(blank) {
            Compartment$new(name = paste0('US_', x$get.name(), '_', blank),
                            type = x)
          })
      }))
    },
    
    
   
    load.compartments = function(settings) {
      ## function creates Compartment objects
      ## `type` attr points directly back to a CompartmentType object, and `name` attr is a unique identifier
      unlist(sapply(names(settings$Compartments), function(comp) {
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
                          branching.time = params$branching.time
          )
        })
        
      }))
    },
    
    
    set.sources = function() {
      ## re-iterates over generated Compartment objects and populates `source` attr with R6 objects
      ## sets 'pointers' to other Compartment objects after generated w/ private$load.compartments()
      compNames <- sapply(private$compartments, function(n){n$get.name()})
      sapply(private$compartments, function(x) {
        if (paste0(x$get.source(),'_1') %in% compNames) {                         # FIXME: arbitrary assignment of source
          searchComps <- which(compNames == paste0(x$get.source(),'_1'))          # FIXME: arbitrarily assigning source to first object in Compartment with $name == x$get.source()
          sourceObj <- private$compartments[[ searchComps ]]    
          x$set.source(sourceObj)
        } # TODO: else statement { if source is 'undefined' or not in the list, must be assigned to an unsampled host (US) }
        x
      })
    },
    
    

    load.lineages = function(settings) {
      ## function creates Lineage objects
      ## `location` attr points directly to a Compartment object, and `name` attr is unique identifier
      ## identifiers create unique Lineages for each Compartment, 
      ## but Compartment A_1 could have Lineage cell_1 and Compartment B_1 also have a separate Lineage cell_1
      unlist(sapply(names(settings$Lineages), function(label) {
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
      
      }))
    },
    
    
   
    get.extant.lineages = function(time) {
      # intializes list of Lineages with sampling.time t=0
      unlist(sapply(private$lineages, function(b){
        if (b$get.sampling.time() <= time) {b}
      }))
    },
    
    
    
    init.popn.growth.dynamics = function(pieces) {
      # @param pieces, list of linear pieces of a given CompartmentType
      # @return matrix, where each row represents a linear piece and comprises the following columns:
        # time = inputted time of the piece
        # popn = inputted popn size of the piece at time `time`
        # slope = slope of the piece 
        # intercept = y-intercept of the piece
      mat <- matrix(nrow=length(pieces), 
                    ncol=6, 
                    dimnames=list(1:length(pieces), c('startTime', 'startPopn', 'endTime', 'endPopn', 'slope', 'intercept')))
      for (x in seq_along(pieces)) {
        if ('startTime' %in% names(unlist(pieces[[x]])) == F) {
          stop ('Parameter "startTime" not defined for piece "', names(pieces)[[x]], '".')
        } else {
          mat[x,1] <- unlist(pieces[[x]])['startTime']
        }
        
        if ('startPopn' %in% names(unlist(pieces[[x]])) == F) {
          stop ('Parameter "startPopn" not defined for piece "', names(pieces)[[x]], '".')
        } else {
          mat[x,2] <- unlist(pieces[[x]])['startPopn']
        }
        
        if ('endTime' %in% names(unlist(pieces[[x]])) == F) {
          mat[x,3] <- NA                                            # for the final piece with `inf` time, assumed that population stays constant from startPopn
          mat[x,4] <- unlist(pieces[[x]])['startPopn']
        } else {
          mat[x,3] <- unlist(pieces[[x]])['endTime']
          if ('endPopn' %in% names(unlist(pieces[[x]])) == F) {
            stop ('Parameter "endPopn" not defined for piece "', names(pieces)[[x]], '".')
          } else {
            mat[x,4] <- unlist(pieces[[x]])['endPopn']
          }
        }
        
      }
      
      # checks for the following:
      
      # only one piece must have a start time of zero
      if ( length(which(mat[,'startTime'] == 0)) != 1 ) {
        stop ('One and only one linear piece in the population growth 
              dynamics functions must have a start time of zero (forward in time).')
      }
      if ( length(which(is.na(mat[,'endTime']))) != 1 ) {
        stop ('One and only one linear piece in the population growth
              dynamics functions must have an end time of infinity (forward in time).
              The rest of the pieces must have attr `endTime` specified for each piece.')
      }
      # all times must be of mode numeric
      if ( is.numeric(mat[,'startTime']) == F || is.numeric(mat[,'endTime']) == F ) {
        stop ('All linear pieces must have start and end times in mode `numeric`.')
      }
      # all times must be unique
      if ( length(unique(mat[,'startTime'])) < length(mat[,'startTime']) ) {
        stop ('All linear pieces must have unique start times.')
      }
      # start and end times of middle pieces must be continuous (with no gaps or overlaps)
      if ( length(intersect(mat[,'startTime'], mat[,'endTime'])) != (nrow(mat)-1) ) {
        stop ('There are time gaps and/or overlapping times within the population growth dynamics pieces.')
      }
      
      # order pieces sequentially based on times (in order to calculate slopes)
      if (nrow(mat) == 1) {
        mat <- as.matrix(t(mat[order(mat[,'startTime']), ]))
        rownames(mat) <- c(1:length(nrow(mat)))
      } else {
        mat <- mat[order(mat[,'startTime']), ]
      }
      
      for (x in 1:nrow(mat)) {
        # calculate slope and intercept for each piece and populate the matrix
        if (x == nrow(mat)) {
          # it is assumed that the final piece will continue towards infinite time at the same population
          # size as the final piece's, with a slope of 0
          slope <- 0
        } else {
          slope <- (mat[x, 'endPopn'] - mat[x, 'startPopn']) / (mat[x, 'endTime'] - mat[x, 'startTime'])
        }
        yInt <- mat[x, 'startPopn'] - slope * mat[x, 'startTime']
        mat[x, 'slope'] <- slope
        mat[x, 'intercept'] <- yInt
      }
      mat
    },
    
    
    
    init.locations = function() {
      # helper function for private$init.pairs()
      # collect host locations of all extant pathogen lineages into dict of host1: [path1, path2, path3, ...]
      private$locations <- list()      # reset list
      for (node in private$extant_lineages) {
        compName <- node$get.location()$get.name()
        if (compName %in% names(private$locations) == F) {
          private$locations[[compName]] <- list()
        }
        private$locations[[compName]] <- c(private$locations[[compName]], node$get.name())
      }
      private$locations
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
      private$choices
    }
    
    
  )
)

