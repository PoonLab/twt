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
#' @param parent:  character, name of parent Pathogen object
#' @param children:  list, names of child Pathogen objects
#' @export
Pathogen <- R6Class(
  "Pathogen",
  public = list(
    initialize = function(name=NA, start.time=NA, end.time=NA, parent=NA,
                          children=list()) {
      private$name <- name
      private$start.time <- start.time
      private$end.time <- end.time
      private$parent <- parent
      private$children <- children
    },
    
    # immutable attributes
    get.name = function() { private$name },
    get.end.time = function() { private$end.time },
    
    # mutables
    get.start.time = function() { private$start.time },
    set.start.time = function(time) { private$start.time <- time },
    get.parent = function() { private$parent },
    set.parent = function(parent) { private$parent <- parent },
    get.children = function() { private$children },
    add.child = function(child) { 
      private$children[[length(private$children)+1]] <- child 
    }
    
  ),
  private = list(
    name = NULL,
    start.time = NULL,
    end.time = NULL,
    parent = NULL,
    children = NULL
  )
)
