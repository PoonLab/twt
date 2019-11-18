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
      names(private$lineages) <- sapply(private$lineages, function(l) {
        l$get.name()
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
      #private$lineages[[length(private$lineages)+1]] <- lineage
      private$lineages[[lineage$get.name()]] <- lineage
      
      comp <- lineage$get.location()
      comp$add.lineage(lineage)
    },
    
    remove.lineage = function(lineage) {
      if ( !is.element('Lineage', class(lineage)) ) {
        stop("Error in remove.lineage: argument 'lineage' must be R6 class Lineage.")
      }
      
      if ( !is.element(lineage$get.name(), names(private$lineages)) ) {
        stop(paste("Error in remove.lineage: ", lineage$get.name(), " not in Run lineages:\n",
             paste(names(private$lineages), collapse=',')))
      }
      
      private$lineages[[lineage$get.name()]] <- NULL
      comp <- lineage$get.location()
      comp$remove.lineage(lineage)
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


#' plot.Run
#' 
#' S3 plot method for Run objects.  
#' 
#' If the Run object is the outcome of sim.outer.tree() or either of the 
#' related outer tree functions, then we call the internal function 
#' .plot.outer.tree().
#' 
#' If the Run object is the outcome of sim.inner.tree(), then the eventlog 
#' is converted into a Phylo object and plotted using the S3 method provided 
#' by the `ape` package.  
#' 
#' @param run:  R6 object of class `Run`
#' @param transmissions:  if TRUE, display transmission events on inner tree
#'        plot.
#' @param migrations: if TRUE, display migration events on inner tree plot.
#' @param node.labels: if TRUE, label internal nodes with names.
#' 
#' @examples 
#' # load model
#' path <- system.file('extdata', 'SI.yaml', package='twt')
#' 
#' # load file and parse to construct MODEL object
#' settings <- yaml.load_file(path)
#' mod <- MODEL$new(settings)
#' 
#' # simulate outer tree
#' outer <- sim.outer.tree(mod)
#' 
#' # simulate inner tree
#' tree <- sim.inner.tree(mod, outer)
#' plot(tree)
#' @export
plot.Run <- function(run, transmissions=FALSE, migrations=FALSE, 
                     node.labels=FALSE, 
                     pal=c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", 
                           "#66A61E", "#E6AB02", "#A6761D", "#666666")) {
  
  eventlog <- run$get.eventlog()
  evt <- eventlog$get.all.events()
  if (is.null(evt)) {
    cat("No events to display.")
  }
  else {
    if ( all(evt$event.type != 'coalescent') ) {
      # this is an outer tree
      .plot.outer.tree(run, pal)
    } 
    else {
      phy <- .eventlogger.to.phylo(
        eventlog=eventlog, 
        transmissions=transmissions, 
        migrations=migrations, 
        node.labels=node.labels)
      plot(phy)  # call plot.phylo S3 method
    }
  }
}

#' .reorder.events
#' 
#' A helper function that recursively geneates a list of node labels
#' by post-order traversal (outputting children before parents).
#' 
#' @keywords internal
.reorder.events <- function(events, node, result=c()) {
  children <- events$compartment1[events$compartment2==node]
  inf.times <- sapply(children, function(child) 
    events$time[events$compartment1==child])
  
  for (child in children[order(inf.times, decreasing=TRUE)]) {
    result <- .reorder.events(events, child, result)
  }
  result <- c(result, node)
  return(result)
}


#' .plot.outer.tree
#' 
#' A function that converts events stored in an EventLogger object into
#' an outer transmission tree.
#' 
#' @param e: EventLogger object
#' 
#' @examples
#' path <- system.file('extdata', 'structSI.yaml', package='twt')
#' settings <- yaml.load_file(path)
#' model <- Model$new(settings)
#' run <- sim.outer.tree(model)
#' plot(run)
#' 
#' @keywords internal
.plot.outer.tree <- function(run) {
  e <- run$get.eventlog()
  comps <- run$get.compartments()
  types <- run$get.types()
  
  events <- e$get.all.events()
  trans <- events[events$event.type=='transmission',]
  migrations <- events[events$event.type=='migration', ]
  
  # find root
  sources <- unique(trans$compartment2)
  root <- sources[!is.element(sources, trans$compartment1)]
  
  # sort nodes by post-order traversal (parents after children)
  nodes <- .reorder.events(trans, root)
  
  # prepare plot region
  par(mar=c(5,1,1,1))
  plot(NA, xlim=c(-max(trans$time*1.05), 0), ylim=c(0.5, length(nodes)+0.5),
       xlab='Time', yaxt='n', ylab=NA, bty='n')
  
  . <- sapply(nodes, function(node) {
    is.sampled <- FALSE
    if (is.element(node, names(comps))) {
      is.sampled <- TRUE
      comp <- comps[[node]]
    }
    
    if (node != root) {
      row <- trans[trans$compartment1 == node,]
      inf.time <- row[['time']]
      source <- row[['compartment2']]
    } else {
      inf.time <- max(trans$time)*1.05
      source <- NA
    }
    
    samp.time <- 0
    if (is.sampled) {
      samp.time <- -max(sapply(comp$get.lineages(), 
                              function(line) line$get.sampling.time()))
    } else {
      samp.time <- min(c(trans$time[trans$compartment2==node],
                         migrations$time[migrations$compartment2==node]))
    }
    
    segments(x0=-inf.time, x1=-samp.time, y0=which(nodes==node), 
             lwd=5, lend=2, col=ifelse(is.sampled, 'black', 'grey'))
    
    if (!is.na(source)) {
      arrows(x0=-inf.time, y1=which(nodes==node), y0=which(nodes==source), 
             length=0.08, lwd=2, col='orangered')  
    }
    
  })
  
  for (i in 1:nrow(migrations)) {
    m <- migrations[i,]
    time <- -m$time
    recipient <- which(nodes == m$compartment1)
    source <- which(nodes == m$compartment2)
    
    arrows(x0=time, y0=source, y1=recipient, col=rgb(0.37, 0.62, 0.63, 0.5), lwd=2, length=0.08)
  }
}
