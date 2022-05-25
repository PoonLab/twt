#' Model
#'
#' \code{Model} is an R6 class that defines an object that generates all of the 
#' objects of the various classes (e.g., Compartment, CompartmentType) that 
#' define a simulation model.  A \code{Model} object is immutable - it should
#' not change over the course of simulation.  Instead, we derive a \code{Run}
#' object class that inherits from a \code{Model}.
#' 
#' @param settings: a named list returned by `yaml.load_file()` that contains 
#' user specifications of the simulation model.
#' 
#' @field initial.conds  Initial Conditions
#' @field types  vector of CompartmentType objects
#' @field compartments  vector of Compartment objects
#' @field lineages  vector of Lineage objects
#' @field fixed.samplings  list of Lineage names and sampling times for plotting
#' 
#' @examples 
#' require(twt)
#' # get path to example YAML file
#' path <- system.file('extdata', 'SI.yaml', package='twt')
#' 
#' # load file and parse to construct Model object
#' settings <- yaml.load_file(path)
#' mod <- Model$new(settings)
#' 
#' # display summary information (calls S3 print method)
#' mod
#' 
#' @export
Model <- R6Class("Model",
  #lock_objects = FALSE,
  
  public = list(
    initialize = function(settings=NA, name=NA) {
      self$name <- name
      
      private$initial.conds <- private$load.initial.conds(settings)
      private$types <- private$load.types(settings)
      private$compartments <- private$load.compartments(settings)
      private$compartments <- private$set.sources()
      private$lineages <- private$load.lineages(settings)
      
      private$fixed.samplings <- private$init.fixed.samplings()
    },
    
    # public variables
    name = NULL,
    
    # ACCESSOR FUNCTIONS
    get.initial.conds = function() { private$initial.conds },
    
    get.types = function() { private$types },
    get.compartments = function() { private$compartments },
    get.lineages = function() { private$lineages },
    
    get.names = function(listR6obj) {
      # returns names of a given list of R6 objects
      # @param listR6obj = list of R6 objects of class CompartmentType, Compartment, or Lineage
      unname(sapply(listR6obj, function(x){x$get.name()}))
    },
    
    get.fixed.samplings = function() {
      # retrieves fixed sampling times of tips
      private$fixed.samplings
    }
    
  ),
  
  
  private = list(
    initial.conds = NULL,
    
    types = NULL,
    compartments = NULL,
    lineages = NULL,
    
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
        result$size[typename] <- floor(size)
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
      
      required <- c('branching.rates', 'transition.rates', 'migration.rates', 
                    'bottleneck.size', 'generation.time')
      
      unlist(sapply(names(settings$CompartmentTypes), function(x) {
        params <- settings$CompartmentTypes[[x]]
        
        missing <- which( !is.element(required, names(params)) )
        if (length(missing) > 0) {
          stop(paste("CompartmentTypes:", x, " missing required field(s): ", required[missing]))
        }
        
        # CompartmentType must specify EITHER one overall rate or a rate for the origin time
        if (!is.null(names(params$branching.rates))) {
          if (!is.element(settings$InitialConditions$originTime, names(params$branching.rates))) {
            stop("Must declare branching rate for origin time ",
                 settings$InitialConditions$originTime,
                 " in CompartmentType ", x)
          }
          
          if (!all(!is.na(as.numeric(names(params$branching.rates))))) {
            wrong_labs <- names(params$branching.rates)[is.na(as.numeric(names(params$branching.rates)))]
            stop("Time-hetetrogeneous branching rates must be declared with numeric time labels '",
                 wrong_labs,
                 "' in CompartmentType ", x, " failed coercion to numeric")
          }
        }
        
        # CompartmentType must specify EITHER effective.size or piecewise linear model
        if (!is.element('effective.size', names(params)) &&
            !is.element('popn.growth.dynamics', names(params))) {
          stop("Either `effective.size` or `popn.growth.dynamics` must be",
               "declared in CoalescentType", x)
        }
        
        # handle epochs (rate changes at specified points in time)
        rate.changes <- lapply(params$branching.rates, function(x) {
          eval(parse(text=paste('list', x))) 
          })
        if (length(rate.changes) > 1) {
          rate.changes <- rate.changes[order(as.numeric(names(rate.changes)), decreasing=T)]
        }
        
        rate.changes2 <- lapply(params$transition.rates, function(x) {
          eval(parse(text=paste('list', x))) 
        })
        
        CompartmentType$new(
          name = x,
          branching.rates = rate.changes,
          transition.rates = rate.changes2,
          migration.rates = eval(parse(text=paste('list', params$migration.rates))),
          
          bottleneck.size = params$bottleneck.size,
          bottleneck.theta = ifelse(is.null(params$bottleneck.theta), 0, 
                                params$bottleneck.theta),
          effective.size = params$effective.size,
          
          popn.growth.dynamics = private$init.popn.growth.dynamics(
            params$popn.growth.dynamics
            ),
          
          generation.time = params$generation.time
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
      
      res <- unlist(sapply(names(settings$Compartments), function(comp) {
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
        
        # parse unsampled field
        if (is.null(params$unsampled)) {
          unsampled <- FALSE  # default to sampled
        }
        else {
          unsampled <- as.logical(params$unsampled)
          if (is.na(unsampled)) {
            warning("Coercing ", params$unsampled, " to NA")
            unsampled <- FALSE
          }
        }
        
        sapply(1:nIndiv, function(obj) {
          Compartment$new(
            name = paste0(comp,'_', obj),  # unique identifier
            type = typeObj,
            source = params$source,
            branching.time = params$branching.time,
            unsampled = unsampled
          )
        })
        
      }))
      
      # return as named vector
      names(res) <- sapply(res, function(comp) comp$get.name())
      res
    },
    
    
    set.sources = function() {
      ## Sets `source` attribute to other Compartment objects after generated 
      ## w/ private$load.compartments().  Note assignment of Compartments
      ## generated through `replicates` setting is not supported.
      ## This is only used when user is manually specifying the outer tree.
      
      compNames <- sapply(private$compartments, function(n){n$get.name()})
      
      sapply(private$compartments, function(x) {
        matches <- which(compNames == paste0(x$get.source(), '_1'))
        if (length(matches) > 0) {
          sourceObj <- private$compartments[[ matches ]]    
          x$set.source(sourceObj)
        } 
        else {
          # source will be set by outer tree simulation
        }
        x
      })
    },
    

    load.lineages = function(settings) {
      ## Parse `Lineages` field of settings
      
      if (is.null(settings$Lineages)) {
        stop("'Lineages' is a required field in settings.")
      }
      
      # return object will be concatenated vector of Lineage objects
      result <- unlist(sapply(names(settings$Lineages), function(label) {
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
            stop('\nSampling times must all be specified as type `numeric`', 
                 ' or type `double`. Sampling time "', vec, '" is of type `', 
                 typeof(vec), '`.\n')
          }
          if (length(sampleTimes) != nIndiv) {
            stop('attribute `sampling.time` of Lineage ', label, 
                 ' does not match number of replicates specified for',
                 'respective Lineage.')
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
            # this also updates Compartment sampling.time attribute
            comp$add.lineage(new.lineage)
            
            new.lineage  # return object
          })
        })
      }))
      
      names(result) <- sapply(result, function(l) l$get.name())
      result
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
        nrow=length(pieces), ncol=5, 
        dimnames=list(1:length(pieces), 
                      c('startTime', 'startPopn', 'endTime', 'endPopn', 'slope'))
        )
      
      for (x in seq_along(pieces)) {
        vec <- unlist(pieces[[x]])  # a named vector
        
        if (!is.element('startTime', names(vec))) {
          stop ('Parameter "startTime" not defined for piece "', names(pieces)[[x]], '".')
        } else {
          mat[x,1] <- vec['startTime']
        }
        
        if (!is.element('startPopn', names(vec))) {
          stop ('Parameter "startPopn" not defined for piece "', names(pieces)[[x]], '".')
        } else {
          mat[x,2] <- vec['startPopn']
        }
        
        if (!is.element('endTime', names(vec))) {
          # for the final piece with `inf` time, assumed that population stays constant from startPopn
          mat[x,3] <- NA
          mat[x,4] <- vec['startPopn']
        } 
        else {
          mat[x,3] <- vec['endTime']
          if (!is.element('endPopn', names(vec))) {
            stop ('Parameter "endPopn" not defined for piece "', names(pieces)[[x]], '".')
          } else {
            mat[x,4] <- vec['endPopn']
          }
        }
        
      }
      
      # checks for the following:
      
      # only one piece must have a start time of zero
      if ( length(which(mat[,'startTime'] == 0)) != 1 ) {
        stop ('One and only one linear piece in the population growth ',
              'dynamics functions must have a start time of zero (forward in time).')
      }
      if ( length(which(is.na(mat[,'endTime']))) != 1 ) {
        stop ('One and only one linear piece in the population growth ', 
              'dynamics functions must have an end time of infinity (forward in time). ',
              'The rest of the pieces must have attr `endTime` specified for each piece.')
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
        stop ('There are time gaps and/or overlapping times within the population ', 
              'growth dynamics pieces.')
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
          slope <- (mat[x, 'startPopn'] - mat[x, 'endPopn']) / 
            (mat[x, 'endTime'] - mat[x, 'startTime'])
        }
        #yInt <- mat[x, 'startPopn'] - slope * mat[x, 'startTime']
        mat[x, 'slope'] <- slope
        #mat[x, 'intercept'] <- yInt
      }
      mat
    },
    
    
    init.fixed.samplings = function() {
      # retrieve sampling time and populate tip labels / times in ape:: phylo object 
      # (for plotting Eventlogger function)
      sapply(private$lineages, function(x) x$get.sampling.time())
    }
    
  )
)


#' print.Model
#' S3 class function to display contents of a Model object
#' @export
print.Model <- function(obj) {
  cat("twt Model", ifelse(is.na(obj$name), '<unnamed>', obj$name), "\n\n")
  
  cat("Initial conditions:\n")
  init <- obj$get.initial.conds()
  cat("  originTime: ", init$originTime, "\n")
  
  cat("  size:\n")
  for (i in 1:length(init$size)) {
    cat(paste0("    ", names(init$size)[i], ": ", init$size[i]), "\n")
  }
  
  cat("  indexType: ", init$indexType, "\n")
  
  cat("CompartmentTypes:\n")
  types <- obj$get.types()
  for (i in 1:length(types)) {
    ct <- types[[i]]
    cat( sprintf("  %10s ", c(names(types)[i])) )
  }
}

#' summary.Model
#' S3 class function to summarize a Model object - simply a wrapper 
#' around print()
#' @export
summary.Model <- function(obj) {
  print(obj)
}
