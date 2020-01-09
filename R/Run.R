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
#' @field lineages:  list of Lineage objects
#' @field extant.lineages:  list of Lineage objects filtered by sampling time
#' @field locations:  a named list of Lineage objects keyed by Compartment
#' @field counts:  a data frame of population dynamics
#' 
#' @export
#' @note https://github.com/r-lib/roxygen2/issues/415
Run <- R6Class("Run",
  #lock_objects = FALSE,
  
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
      #private$extant.lineages <- private$retrieve.extant.lineages(0)
      private$locations <- private$locate.lineages()
      
      # placeholder to be populated by sim.outer.tree()
      # not used if user supplies transmission history
      private$unsampled.hosts <- list()
      },

    
    ## ACCESSOR FUNCTIONS
    get.eventlog = function() { private$eventlog },
    set.eventlog = function(e) { private$eventlog <- e },
    
    get.counts = function() { private$counts },
    set.counts = function(df) { private$counts <- df },
    
    get.initial.conds = function() { private$initial.conds },
    get.fixed.samplings = function() { private$fixed.samplings },
    
    get.types = function() { private$types },
    get.compartments = function() { private$compartments },
    get.lineages = function() { private$lineages },
    get.unsampled.hosts = function() { private$unsampled.hosts },
    
    get.extant.lineages = function(time, comp=NA) {
      #' Retrieves a list of Lineage objects that are extant at the 
      #' specified time.
      #' Note that because Lineages are being added and removed through
      #' coalescent and bottleneck events, the output will change
      #' over the course of inner tree simulation.
      #' 
      #' @param time: coalescent (cumulative time) of the simulation
      #' @param comp: R6 object of class Compartment
      #' 
      #' @return  a named list of Lineage objects.
      #' 
      if (is.environment(comp)) {
        extant.lineages <- unlist(sapply(comp$get.lineages(), function(b) {
          if (is.element(b$get.name(), names(private$lineages)) & 
              b$get.sampling.time() <= time) { b }
        }))
      }      
      else {
        extant.lineages <- unlist(sapply(private$lineages, function(b){
          if (b$get.sampling.time() <= time) { b }
        }))  
      }
      
      return(extant.lineages)
    },
    
    add.lineage = function(lineage) {
      # Add a Lineage object to private list.  Used to track ancestral 
      # lineages produced from a coalescent event.
      # 
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
      # 
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
    
    
    get.node.ident = function(prefix="Node") {
      # returns unique identity for internal nodes (inner tree sim, ancestral lineages)
      result <- paste(prefix, private$node.ident, sep='')
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
    eventlog = NULL,
    initial.conds = NULL,
    types = NULL,
    compartments = NULL,
    
    unsampled.hosts = NULL,
    lineages = NULL,
    fixed.samplings = NULL,
    
    locations = NULL,
    counts = NULL,
    
    node.ident = 1,  # used in simulation of inner tree for generating 
                     # unique idents for internal nodes of ancestral lineages
    
    
    # private functions
    locate.lineages = function() {
      # Report host locations of all extant lineages into a list,
      # given last time user called `get.extant.lineages`
      #
      # @return  List of form { host1 = c(line1, line2, line3), ... }
      
      private$locations <- list()      # reset list
      for (node in private$extant.lineages) {
        compName <- node$get.location()$get.name()
        
        if ( !is.element(compName, names(private$locations)) ) {
          private$locations[[compName]] <- list()
        }
        private$locations[[compName]] <- c(private$locations[[compName]], node$get.name())
      }
      private$locations
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
#' @param type:  what type of plot should be drawn.  Possible types are
#' \itemize{
#' \item "t" for transmission tree plot
#' \item "s" for stair step plot of population dynamics
#' }
#' 
#' @examples 
#' # load model
#' path <- system.file('extdata', 'SI.yaml', package='twt')
#' 
#' # load file and parse to construct MODEL object
#' settings <- yaml.load_file(path)
#' mod <- Model$new(settings)
#' 
#' # simulate outer tree
#' run <- sim.outer.tree(mod)
#' 
#' # generate plots side-by-side
#' par(mfrow=c(1,2))
#' plot(run, type='s')  # population dynamics
#' plot(run, type='t')  # transmission tree
#' 
#' @export
plot.Run <- function(run, type='t', ...) {
  eventlog <- run$get.eventlog()
  evt <- eventlog$get.all.events()
  if (is.null(evt)) {
    cat("No events to display.")
  }
  
  if (type=='t') {
    .plot.outer.tree(run, ...)
  }
  else if (type == 's') {
    .plot.dynamics(run, ...)
  }
  else {
    stop("Error: type '", type, "' is not recognized by plot.Run(). ",
         "Use 't' for outer tree plot or 's' for a stairstep plot ",
         "of the population dynamics.")
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
.plot.outer.tree <- function(run, ...) {
  e <- run$get.eventlog()
  comps <- c(run$get.compartments(), run$get.unsampled.hosts())
  types <- run$get.types()
  
  events <- e$get.all.events()
  trans <- events[events$event.type=='transmission',]
  migrations <- events[events$event.type=='migration', ]
  
  # find root
  sources <- unique(trans$compartment2)
  root <- sources[!is.element(sources, trans$compartment1)]
  if (length(root) > 1) {
    stop("Detected multiple roots in Run object: ", 
         paste(root, collapse=", "))
  }
  
  # sort nodes by post-order traversal (parents after children)
  nodes <- .reorder.events(trans, root)
  
  # prepare plot region
  par(mar=c(5,1,1,1))
  plot(NA, xlim=c(-max(trans$time*1.05), 0), ylim=c(0.5, length(nodes)+0.5),
       xlab='Time', yaxt='n', ylab=NA, bty='n', ...)
  
  #. <- sapply(nodes, function(node) {
  for (i in 1:length(nodes)) {
    node <- nodes[i]
    comp <- comps[[node]]
    is.sampled <- !comp$is.unsampled()
    
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
      # unsampled compartment - right limit at first transmission or migration event
      if (nrow(migrations) > 0) {
        samp.time <- min(c(trans$time[trans$compartment2==node],
                           migrations$time[migrations$compartment2==node]))  
      } else {
        if (is.element(node, trans$compartment2)) {
          samp.time <- trans$time[trans$compartment2 == node]  
        }
      }
    }
    
    segments(x0=-inf.time, x1=-samp.time, y0=which(nodes==node), 
             lwd=5, lend=2, col=ifelse(is.sampled, 'black', 'grey'))
    
    if (!is.na(source)) {
      arrows(x0=-inf.time, y1=which(nodes==node), y0=which(nodes==source), 
             length=0.08, lwd=2, col='orangered')  
    }
  }
  
  if (nrow(migrations) > 0) {
    # map migration events to transmission tree
    for (i in 1:nrow(migrations)) {
      m <- migrations[i,]
      time <- -m$time
      recipient <- which(nodes == m$compartment1)
      source <- which(nodes == m$compartment2)
      
      arrows(x0=time, y0=source, y1=recipient, col=rgb(0.37, 0.62, 0.63, 0.5), lwd=2, length=0.08)
    } 
  }
}


#' .plot.dynamics
#' Make a stairstep plot of the population dynamics.
#' @param run:  R6 object of class Run
#' @param ...:  additional arguments passed to generic plot()
#' 
#' @keywords internal
.plot.dynamics <- function(run, ...) {
  counts <- run$get.counts()
  if (is.null(counts)) {
    stop("Run object does not have counts field, may have been loaded from YAML")
  }
  
  # set up plot region
  par(mar=c(5,5,1,5), xpd=NA)
  plot(NA, xlim=rev(range(counts$time)), ylim=range(counts[,-1]), bty='n',
       xlab='Coalescent time', ylab='Population size', ...)
  for (i in 2:ncol(counts)) {
    label <- names(counts)[i]
    ctype <- gsub("[SI]\\.(.+)", "\\1", label)
    lines(counts$time, counts[,i], type='s', 
          col=ifelse(startsWith(label, 'S'), 'dodgerblue', 'firebrick2'),
          lwd=2)
    text(x=0, y=counts[nrow(counts), i], label=label, adj=0)
  }
}


