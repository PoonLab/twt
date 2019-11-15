#' Run
#' 
#' \code{Run} is an R6 class that is derived from \code{Model}.  It inherits
#' all simulation model features including CompartmentTypes, Compartments, 
#' and Lineages.  Unlike \code{Model}, a \code{Run} object is used to track 
#' an outcome of the simulation model and is therefore mutable where \code{Model}
#' is immutable.  For example, \code{Run} tracks which lineages are extant and 
#' the numbers of pairs of Lineages that can coalesce over time.  In addition, 
#' the number of unsampled compartments is stochastic outcome of the model that
#' is tracked at this level.
#' 
#' @param model: an object of class Model
#' 
#' @field extant.lineages  list of Lineage objects with sampling time t=0 (most recent)
#' @field locations  a named list of Lineage objects keyed by Compartment
#' @field choices  list of Lineage pairs that may coalesce per Compartment, 
#' keyed by index to `locations`
#' 
#' @export
Run <- R6Class(
  "Run",
  lock_objects = FALSE,
  
  public = list(
    
    initialize = function(model) {
      private$eventlog <- EventLogger$new()
      
      # transfer objects from Model
      private$initial.conds <- model$get.initial.conds()
      private$types <- model$get.types()
      
      
      # deep copy of Compartments and Lineages
      private$compartments <- lapply(model$get.compartments(), function(comp) {
        comp$copy(deep=TRUE)
      })
      
      private$lineages <- unlist(
        lapply(private$compartments, function(x) x$get.lineages())
        )
      
      
      private$fixed.samplings <- model$get.fixed.samplings()
      
      #private$extant.lineages <- private$get.extant.lineages(0)
      private$extant.lineages <- private$retrieve.extant.lineages(0)
      private$locations <- private$init.locations()
      private$choices <- private$init.pairs()
      
      # placeholder to be populated by sim.outer.tree()
      # not used if user supplies transmission history
      private$unsampled.hosts <- list()
      },

    
    ## ACCESSOR FUNCTIONS
    get.eventlog = function() { private$eventlog },
    set.eventlog = function(e) { private$eventlog <- e },
    
    get.initial.conds = function() { private$initial.conds },
    get.fixed.samplings = function() { private$fixed.samplings },
    
    get.types = function() { private$types },
    get.compartments = function() { private$compartments },
    get.lineages = function() { private$lineages },
    get.unsampled.hosts = function() { private$unsampled.hosts },
    
    get.extant.lineages = function(time) {
      # Retrieves a list of Lineage objects that are extant at the 
      # specified time.
      # Caches the result in member variable 'extant.lineages'
      # 
      # @param time = coalescent (cumulative time) of the simulation
      private$extant.lineages <- private$retrieve.extant.lineages(time)
      private$extant.lineages
    },
    
    
    get.pairs = function() {
      # function extracts and returns all the current pairs of pathogen lineages
      # that may coalesce
      private$choices
    },
    
    
    
    ## SIMULATION FUNCTIONS
    
    add.pair = function(L1, L2, host) {
      # adds a pair of pathogen lineages that may coalesce within a given 
      # compartment into list of `choices`
      # A. when a Lineage is moved from one compartment to another (transmission 
      #    or migration)
      # B. when a Lineage is sampled
      # C. can also be used to update the location of a pair
      # 
      # @param L1, L2 = Lineage objects
      # @param host = Compartment object
      pair <- sort(c(L1, L2))
      private$choices[[paste(pair[1], pair[2], sep=',')]] <- host
    },
    
    remove.pair = function(L1, L2) {
      # removes a pair of pathogen lineages that can no longer coalesce from 
      # the list of `choices`
      # A. when a coalescence occurs
      # B. when Lineages reach a transmission bottleneck, forcing coalescence
      # 
      # @param L1, L2 = <character> names of Lineage objects
      
      if (!is.character(L1) || !is.character(L2)) {
        stop("Error in remove.pair(): arguments must be character values.")
      }
      pair <- sort(c(L1, L2))
      private$choices[[paste(pair[1], pair[2], sep=',')]] <- NULL
    }, 
    
    add.lineage = function(lineage) {
      # @param lineage: an R6 object of class Lineage
      if ( !is.element('Lineage', class(lineage)) ) {
        stop("Error in add.lineage(): <lineage> must be R6 class Lineage object.")
      }
      private$lineages[[length(private$lineages)+1]] <- lineage
    },
    
    remove.lineage = function(lineage) {
      if ( !is.element('Lineage', class(lineage)) ) {
        stop("Error in remove.lineage: argument 'lineage' must be R6 class Lineage.")
      }
      temp <- sapply(private$lineages, function(x) x$get.name())
      target <- lineage$get.name()
      if (!is.element(target, temp)) {
        stop("Error in remove.lineage(): ", target, " not in Run lineages:\n",
             paste(temp, collapse=','))
      }
      
      private$lineages <- private$lineages[-which(temp==target)]
    },
    
    
    get.node.ident = function() {
      # returns unique identity for internal nodes (inner tree sim, ancestral lineages)
      result <- paste("Node", private$node.ident, sep='')
      private$node.ident <- private$node.ident + 1
      return(result)
    },
    
    
    generate.unsampled = function(type.name) {
      # function creates "blank" Compartment objects for Unsampled Hosts (US)
      # @param type.name = unique name of CompartmentType
      
      new.host <- Compartment$new(
        name=paste0('US_', type.name, '_', length(private$unsampled.hosts)+1),
        type=private$types[[type.name]], 
        unsampled=TRUE
        )
      
      private$unsampled.hosts[[new.host$get.name()]] <- new.host
      return(new.host)
    },
    
    clear.unsampled = function() {
      private$unsampled.hosts <- NULL
    }
    
  ),

  
  private = list(
    
    unsampled.hosts = NULL,
    extant.lineages = NULL,
    
    locations = NULL,         
    choices = NULL,
    
    node.ident = 1,  # used in simulation of inner tree for generating 
                     # unique idents for internal nodes of ancestral lineages
    
    
    # private functions
    
    retrieve.extant.lineages = function(time) {
      # intializes list of Lineages with sampling.time t=0
      unlist(sapply(private$lineages, function(b){
        if (b$get.sampling.time() <= time) {b}
      }))
    },
    
    
    init.locations = function() {
      # Helper function for private$init.pairs()
      # collect host locations of all extant lineages into a list 
      # of host1: [line1, line2, line3, ...]
      
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
      # Extract all pairs of extant lineages that may coalesce
      # 
      # @return:  list keyed by Lineage name tuples that are associated
      # to host names.
      
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
    }
  )
  
)


#' Wrapper around plot.EventLogger
#' @export
plot.Run <- function(run, transmissions=FALSE, migrations=FALSE, 
                     node.labels=FALSE) {
  # call S3 method of this Run's event log
  plot(run$get.eventlog())
}
