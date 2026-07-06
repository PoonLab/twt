#' Pathogen
#' 
#' \code{Pathogen} is an R6 class for objects that represent pathogen lineages
#' that are carried by Hosts and which comprise the "inner" tree of the 
#' simulation.
#' 
#' @param name:  character, unique identifier for Pathogen object
#' @param start.time:  numeric, start time of lineage; a coalescent event 
#'        terminates the parent and starts two child lineages.
#' @param end.time:   numeric, end time of lineage; associated with a 
#'        coalescent event (if parent) or sampling event.
#' @param parents:  list, parent Pathogen objects; length 1 for transmission,
#'        length 2 for recombination
#' @param children:  list, child Pathogen objects (NOTE: avoid cloning 
#'        Pathogen objects or there will be a circular reference problem!)
#' @param breakpoint:  integer, genomic position of recombination breakpoint;
#'        NA for non-recombinant lineages
#' @export
Pathogen <- R6Class(
  "Pathogen",
  public = list(
    initialize = function(name=NA, start.time=NA, end.time=NA, parents=list(),
                          children=list(), breakpoint=NA) {
      private$name <- name
      private$start.time <- start.time
      private$end.time <- end.time
      private$parents <- parents
      private$children <- children
      private$breakpoint <- breakpoint
    },
    
    print = function() {
      cat("twt Pathogen", self$get.name(), "\n")
      cat("  Start time:", self$get.start.time(), "\n")
      cat("  End time:", self$get.end.time(), "\n")
      cat("  Parents:", sapply(self$get.parents(), function(p) p$get.name()), "\n")
      cat("  Children:", sapply(self$get.children(), function(c) c$get.name()), "\n")
      cat("  Breakpoint:", self$get.breakpoint(), "\n")
    },
    
    # immutable attributes
    get.name = function() { private$name },
    get.end.time = function() { private$end.time },
    
    # mutables
    get.start.time = function() { private$start.time },
    set.start.time = function(time) { private$start.time <- time },
    get.parents = function() { private$parents },
    get.parent = function() { private$parents[[1]] },
    set.parent = function(parent) { private$parents <- list(parent) },
    add.parent = function(parent) {
      private$parents[[length(private$parents)+1]] <- parent
    },
    get.children = function() { private$children },
    add.child = function(child) { 
      private$children[[length(private$children)+1]] <- child 
    },
    get.breakpoint = function() { private$breakpoint },
    set.breakpoint = function(bp) { private$breakpoint <- bp }
    
  ),
  private = list(
    name = NULL,
    start.time = NULL,
    end.time = NULL,
    parents = NULL,
    children = NULL,
    breakpoint = NULL
  )
)
