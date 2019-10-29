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
      # transfer objects from Model
      private$initial.conds <- model$get.initial.conds()
      private$types <- model$get.types()
      
      # deep copy of Compartments and Lineages
      private$compartments <- lapply(model$get.compartments(), function(comp) {
        comp$clone()
      })
      private$lineages <- lapply(model$get.lineages(), function(line) {
        line$clone()
      })
      
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
    get.initial.conds = function() { private$initial.conds },
    get.fixed.samplings = function() { private$fixed.samplings },
    get.types = function() { private$types },
    get.compartments = function() { private$compartments },
    get.lineages = function() { private$lineages },
    get.unsampled.hosts = function() { private$unsampled.hosts },
    
    get.extant.lineages = function(time) {
      # returns lineages extant at a given time
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
      # @param L1, L2 = Lineage objects
      pair <- sort(c(L1, L2))
      private$choices[[paste(pair[1], pair[2], sep=',')]] <- NULL
    }, 
    
    add.lineage = function(ancestral.lineage) {
      # at a coalescent event, an ancestral lineage must be created
      private$lineages[[length(private$lineages)+1]] <- ancestral.lineage
    },
    
    remove.lineage = function(lineage) {
      # at a coalescent event, lineages that coalesce must be removed
      lin.indices <- which(sapply(
        private$lineages, 
        function(x) { x$get.name() == lineage$get.name() }
        ))
      private$lineages <- private$lineages[-lin.indices]
    },
    
    
    get.node.ident = function() {
      # returns unique identity for internal nodes (inner tree sim, ancestral lineages)
      private$node.ident
    },
    
    update.node.ident = function() {
      # generates new unique identity for next `$get.node.ident()` call 
      # (inner tree sim, internal ancestral lineages)
      private$node.ident <- private$node.ident + 1
    }, 
    
    
    generate.unsampled = function(num.unsampled, t) {
      # function creates "blank" Compartment objects for Unsampled Hosts (US)
      # @param num.unsampled = number of unsampled
      # @param t = CompartmentType object
      new.hosts <- unlist(sapply(
        1:num.unsampled, 
        function(blank) {
          Compartment$new(
            name=paste0('US_', t$get.name(), '_', blank),
            type=t,
            unsampled=TRUE
          )
        }))
      names(new.hosts) <- sapply(new.hosts, function(x) x$get.name())
  
      private$unsampled.hosts <- c(private$unsampled.hosts, new.hosts)
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
    }
  )
  
)


