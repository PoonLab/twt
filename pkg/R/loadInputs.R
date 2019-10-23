#' MODEL
#'
#' \code{MODEL} is an R6 class that defines an object that generates all of the 
#' objects of the various classes (e.g., Compartment, CompartmentType) that 
#' define a simulation model.
#' 
#' @param settings: a named list returned by `yaml.load_file()` that contains 
#' user specifications of the simulation model.
#' 
#' @field initial.conds  Initial Conditions
#' @field types  vector of CompartmentType objects
#' @field compartments  vector of Compartment objects
#' @field lineages  vector of Lineage objects
#' @field extant.lineages  list of Lineage objects with sampling time t=0 (most recent)
#' @field locations  a named list of Lineage objects keyed by Compartment
#' @field choices  list of Lineage pairs that may coalesce per Compartment, 
#' keyed by index to `locations`
#' @field fixed.samplings  list of Lineage names and sampling times for plotting
#' 
#' @examples 
#' require(twt)
#' # get path to example YAML file
#' path <- system.file('extdata', 'SI.yaml', package='twt')
#' 
#' # load file and parse to construct MODEL object
#' settings <- yaml.load_file(path)
#' mod <- MODEL$new(settings)
#' 
#' # display the first compartment
#' mod$get.compartments()[1]
#' 
#' @export
MODEL <- R6Class("MODEL",
  public = list(
    initialize = function(settings=NA) {
      private$initial.conds <- private$load.initial.conds(settings)
      private$types <- private$load.types(settings)
      private$compartments <- private$load.compartments(settings)
      private$compartments <- private$set.sources()
      private$lineages <- private$load.lineages(settings)
      
      private$extant.lineages <- private$retrieve.extant.lineages(0)
      private$locations <- private$init.locations()
      private$choices <- private$init.pairs()
      private$fixed.samplings <- private$init.fixed.samplings()
    },
    
    get.initial.conds = function() {private$initial.conds},
    get.types = function() {private$types},
    get.unsampled.hosts = function() {private$unsampled.hosts},  # populated by sim.outer.tree
    get.compartments = function() {private$compartments},
    get.lineages = function() {private$lineages},
    
    get.extant.lineages = function(time) {
      # returns lineages extant at a given time
      # @param time = coalescent (cumulative time) of the simulation
      private$extant.lineages <- private$retrieve.extant.lineages(time)
      private$extant.lineages
    },
    
    get.names = function(listR6obj) {
      # returns names of a given list of R6 objects
      # @param listR6obj = list of R6 objects of class CompartmentType, Compartment, or Lineage
      unname(sapply(listR6obj, function(x){x$get.name()}))
    },
    
    get.pairs = function() {
      # function extracts and returns all the current pairs of pathogen lineages that may coalesce
      private$choices
    },
    
    add.pair = function(L1, L2, host) {
      # FIXME: this and subsequent functions involve manipulation of private data
      # which means that a MODEL changes state after simulation.  Make a derived class?
      
      # adds a pair of pathogen lineages that may coalesce within a given compartment into list of `choices`
        # A. when a Lineage is moved from one compartment to another (transmission or migration)
        # B. when a Lineage is sampled
        # C. can also be used to update the location of a pair
      # @param L1,L2 = Lineage objects
      # @param host = Compartment object
      pair <- sort(c(L1, L2))
      private$choices[[paste(pair[1], pair[2], sep=',')]] <- host
    },
    
    remove.pair = function(L1, L2) {
      # removes a pair of pathogen lineages that can no longer coalesce from the list of `choices`
        # A. when a coalescence occurs
        # B. when Lineages reach a transmission bottleneck, forcing coalescence
      # @param L1,L2 = Lineage objects
      pair <- sort(c(L1, L2))
      private$choices[[paste(pair[1], pair[2], sep=',')]] <- NULL
    }, 
    
    add.lineage = function(ancestral.lineage) {
      # at a coalescent event, an ancestral lineage must be created
      private$lineages[[length(private$lineages)+1]] <- ancestral.lineage
    },
    
    remove.lineage = function(lineage) {
      # at a coalescent event, lineages that coalesce must be removed
      lin.indices <- which(sapply(private$lineages, function(x){x$get.name() == lineage$get.name()}))
      private$lineages <- private$lineages[-lin.indices]
    },
    
    get.node.ident = function() {
      # returns unique identity for internal nodes (inner tree sim, ancestral lineages)
      private$node.ident
    },
    
    update.node.ident = function() {
      # generates new unique identity for next `$get.node.ident()` call (inner tree sim, internal ancestral lineages)
      private$node.ident <- private$node.ident + 1
    }, 
    
    generate.unsampled = function(num.unsampled, t) {
      # function creates "blank" Compartment objects for Unsampled Hosts (US)
      # @param num.unsampled = number of unsampled
      # @param t = CompartmentType object
      private$unsampled.hosts <- c(private$unsampled.hosts, 
                                   unlist(sapply(1:num.unsampled, function(blank) {
        Compartment$new(name=paste0('US_', t$get.name(), '_', blank),
                        type=t,
                        unsampled=TRUE)
      }))) 
    },
    
    clear.unsampled = function() {
      private$unsampled.hosts <- NULL
    },
    
    get.fixed.samplings = function() {
      # retrieves fixed sampling times of tips
      private$fixed.samplings
    }
    
  ),
  
  
  private = list(
    initial.conds = NULL,
    types = NULL,
    unsampled.hosts = NULL,
    compartments = NULL,
    lineages = NULL,
    
    extant.lineages = NULL,
    
    locations = NULL,         
    choices = NULL,
    node.ident = 1,  # used in simulation of inner tree for generating 
                     # unique idents for internal nodes of ancestral lineages
    fixed.samplings = NULL,

        
    load.initial.conds = function(settings) {
      # Loads initial conditions
      # @return list object
      result <- list()
      
      # origin time is measured in reverse (prior to most recent sampled lineage)
      params <- settings$InitialConditions #[[x]]
      if (!is.numeric(params$originTime)) {
        stop("InitialConditions:originTime must be numeric; is this key missing?")
      }
      if (params$originTime <= 0) {
        stop("InitialConditions:originTime must be >0")
      }
      result['originTime'] <- params$originTime
      
      
      # load initial numbers of Compartments per Type
      if (is.null(settings$CompartmentType)) {
        stop("InitialConditions:CompartmentType key missing")
      }
      comp.types <- names(settings$CompartmentTypes)
      result$size <- list()
      for (i in 1:length(params$size)) {
        typename <- names(params$size)[i]
        if (!is.element(typename, comp.types)) {
          stop(paste("InitialConditions:size key", typename, 
                     "must match CompartmentType"))
        }
        size <- params$size[[typename]]
        if (!is.numeric(size)) {
          stop(paste("InitialConditions:size:", typename, " not numeric"))
        }
        result$size[typename] <- size
      }
       
      
      # specify CompartmentType of index case
      if (is.null(params$indexType)) {
        stop("Settings must specify 'indexType' CompartmentType under InitialConditions")
      }
      if (!is.element(params$indexType, comp.types)) {
        stop("'indexType' in InitialConditions does not match any CompartmentType in settings")
      }
      result['indexType'] <- params$index
      
      result
    },
    
    
    load.types = function(settings) {
      # function creates CompartmentType objects
      # within each CompartmentType, there are distinct compartments with
      # individual transmission & migration rates
      # @return: a named vector of CompartmentType R6 objects
      
      if (is.null(settings$CompartmentTypes)) {
        stop("Missing 'CompartmentTypes' field in settings")
      }
      if (length(settings$CompartmentTypes) == 0) {
        stop("Empty 'CompartmentTypes' field in settings.")
      }
      
      required <- c('branching.rates', 'migration.rates', 'bottleneck.size',
                    'coalescent.rate', 'wait.time.distr')
      
      unlist(sapply(names(settings$CompartmentTypes), function(x) {
        params <- settings$CompartmentTypes[[x]]
        
        missing <- which( !is.element(required, names(params)) )
        if (length(missing) > 0) {
          stop(paste("CompartmentTypes:", x, " missing required field(s): ", required[missing]))
        }
        
        # Generate wait time distribution between Compartment infection time and 
        # its first sampling time
        # code below is based directly from Poonlab/Kaphi/pkg/R/smcConfig.R
        sublist <- params$wait.time.distr
        rng.call <- paste('d', sublist$dist, '(x,', sep='')
        args <- sapply(sublist[['hyperparameters']], function(x) paste(names(x), x, sep='='))
        rng.call <- paste(rng.call, paste(args, collapse=','), ')', sep='')
        
        CompartmentType$new(
          name = x,
          branching.rates = eval(parse(text=paste('list', params$branching.rates))),
          migration.rates = eval(parse(text=paste('list', params$migration.rates))),
          
          bottleneck.size = params$bottleneck.size,
          coalescent.rate = params$coalescent.rate,
          
          death.rate.distr = params$death.rate.distr,  # not implemented yet
          wait.time.distr = rng.call,
          
          popn.growth.dynamics = private$init.popn.growth.dynamics(
            params$popn.growth.dynamics
            )
          # transmission.times parameter populated later when simulating outer tree 
          # (simOuterTree.R) to be used when simulating migration events in inner 
          # tree (simInnerTree.R)
        )
      }))
    },
    
    
   
    load.compartments = function(settings) {
      ## function creates Compartment objects
      ## `type` attr points directly back to a CompartmentType object, and 
      ## `name` attr is a unique identifier
      if (is.null(settings$Compartments)) {
        stop("Missing 'Compartments' field in settings.")
      }
      unlist(sapply(names(settings$Compartments), function(comp) {
        if (grepl("_", comp)) {
          stop("Error: Underscore characters are reserved, please modify Compartment name", comp)
        }
        
        params <- settings$Compartments[[comp]]
        if (is.null(params$type)) {
          stop(paste("Compartment", comp, "missing required field 'type'."))
        }
        
        if (params$type %in% names(settings$CompartmentType)) {
          searchTypes <- which(names(settings$CompartmentType) == params$type)
          # reference to CompartmentType object
          typeObj <- private$types[[ searchTypes ]]
        }
        else {
          stop(params$type, ' of Compartment ', comp, 
               ' is not a specified Compartment Type object')
        }
        
        if (is.null(params$replicates)) {
          nIndiv <- 1
        } else {
          nIndiv <- params$replicates
        }
        
        sapply(1:nIndiv, function(obj) {
          Compartment$new(
            name = paste0(comp,'_', obj),  # unique identifier
            type = typeObj,
            source = params$source,        
            branching.time = params$branching.time
          )
        })
        
      }))
    },
    
    
    set.sources = function() {
      ## Sets `source` attribute to other Compartment objects after generated 
      ## w/ private$load.compartments().  Note assignment of Compartments
      ## generated through `replicates` setting is not supported.
      
      compNames <- sapply(private$compartments, function(n){n$get.name()})
      
      sapply(private$compartments, function(x) {
        matches <- which(compNames == paste0(x$get.source(), '_1'))
        if (length(matches) > 0) {
          sourceObj <- private$compartments[[ matches ]]    
          x$set.source(sourceObj)
        } 
        else {
          # TODO: if source is 'undefined' or not in the list, must be 
          # assigned to an unsampled host (US)
        }
        x
      })
    },
    
    

    load.lineages = function(settings) {
      ## Parse `Lineages` field of settings
      
      if (is.null(settings$Lineages)) {
        stop("'Lineages' is a required field in settings.")
      }
      
      unlist(sapply(names(settings$Lineages), function(label) {
        if (grepl("_", label)) {
          stop("Error: underscore characters are reserved, please modify Lineage name",
               label)
        }
         
        params <- settings$Lineages[[label]]
        if (is.null(params$location)) {
          stop("Lineage", label, "missing required field 'location'")
        }
        
        
        # find all Compartments that match the specified location
        if (length(private$compartments) == 0) {
          stop("Compartments must be parsed before Lineages.")
        }
        
        # FIXME: this should work..
        #compNames <- self$get.names(Compartment)
        compNames <- sapply(private$compartments, function(n){n$get.name()})
        
        searchComps <- which(grepl(paste0("^", params$location, "_"), compNames))
        if (length(searchComps) == 0) {
          stop(params$location, ' of Lineage ', label, 
               ' is not a specified Compartment object.\n',
               'Available compartments:\n', 
               paste(compNames))
        }
        
        
        if (is.null(params$replicates)) {
          nIndiv <- 1
        } else {
          nIndiv <- params$replicates  
        }
        
        
        if (is.null(params$sampling.time)) {
          stop("Lineage", label, "missing required field 'sampling.time'")
        }
        if (is.numeric(params$sampling.time)) {
          sampleTimes <- rep.int(params$sampling.time, times=nIndiv)
        } 
        else {
          vec <- unlist(strsplit(params$sampling.time, split='[`(`|,|`)`]'))
          # exclude empty strings
          sampleTimes <- as.double(vec[nzchar(x=vec)])
          if (any(is.na(sampleTimes))) {
            # failed to parse string as double
            stop('\nSampling times must all be specified as type `numeric` or type `double`. ', 
                 'Sampling time "', vec, '" is of type `', typeof(vec), '`.\n')
          }
          if (length(sampleTimes) != nIndiv) {
            stop('attribute `sampling.time` of Lineage ', label, 
                 ' does not match number of replicates specified for respective Lineage.')
          }
        }

                
        sapply(searchComps, function(compIdx) {
          comp <- private$compartments[[compIdx]]
          
          sapply(1:nIndiv, function(i) {            
            # generate unique identifier
            new.lineage <- Lineage$new(
              # for example, "host_3__virus_2"
              name = paste0(comp$get.name(),'__',label,'_',i),
              type = params$type,  # LineageType, UNUSED
              sampling.time = sampleTimes[i],
              location = comp  # reference to Compartment object
            )
            # add Lineage to Compartment
            comp$add.lineage(new.lineage)
            
            new.lineage  # return object
          })
        })
      
      }))
    },
    
    
   
    retrieve.extant.lineages = function(time) {
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
      
      if (is.null(pieces)) {
        return (NULL)
      }
      mat <- matrix(
        nrow=length(pieces), ncol=6, 
        dimnames=list(1:length(pieces), 
                      c('startTime', 'startPopn', 'endTime', 'endPopn', 
                        'slope', 'intercept'))
        )
      
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
          # for the final piece with `inf` time, assumed that population stays constant from startPopn
          mat[x,3] <- NA
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
      # collect host locations of all extant pathogen lineages into dict 
      # of host1: [path1, path2, path3, ...]
      private$locations <- list()      # reset list
      for (node in private$extant.lineages) {
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
      
      # reset container
      private$choices <- list()
      
      for (hostNum in seq_along(private$locations)) {
        # retrieve vector of Lineages in this host
        lineages <- private$locations[[hostNum]]
        hostname <- names(private$locations)[[hostNum]]
        
        if (length(lineages) > 1) {
          combns <- combn(lineages, 2)
          
          for (col in 1:ncol(combns)) {
            pair <- sort(unlist(combns[,col]))
            key <- paste(pair, collapse=',')
            private$choices[[key]] <- hostname
          }
          #. <- apply(combns, 2, function(pair) {
          #  key <- paste(sort(unlist(pair)), collapse=',')
          #  private$choices[[key]] <- hostname
          #})
        }
      }
      private$choices
    },
    
    
    
    init.fixed.samplings = function() {
      # retrieve sampling time and populate tip labels / times in ape:: phylo object 
      # (for plotting Eventlogger function)
      list(
        # store label w/ corresponding tip height in new ape::phylo object (not casted into `phylo` yet)
        tip.label = sapply(private$lineages, function(x) x$get.name()),
        
        # only used for calculating edge length
        tip.height = sapply(private$lineages, function(x) x$get.sampling.time())
      )
    }
    
    
  )
)


print.MODEL <- function(obj) {
  cat("Initial conditions:\n")
  cat("  originTime: ", obj$origin, "\n")
  
}
