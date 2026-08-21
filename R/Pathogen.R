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
#' @param carrier.id:  identifier for the physical host individual this
#'        lineage's genome currently occupies. Distinct from the Pathogen's
#'        own name/lineage identity: recombination splits one lineage into
#'        two Pathogen objects that are still the same physical genome, so
#'        both inherit the same carrier.id. New (non-recombinant) Pathogens
#'        get a fresh carrier.id. Used to count physical individuals rather
#'        than tracked lineages (see InnerTree registry / .do.infection),
#'        so recombination-driven lineage growth doesn't inflate apparent
#'        population size.
#' @export
Pathogen <- R6Class(
  "Pathogen",
  public = list(
    initialize = function(name=NA, start.time=NA, end.time=NA, parents=list(),
                          children=list(), breakpoint=NA, carrier.id=NA) {
      private$name <- name
      private$start.time <- start.time
      private$end.time <- end.time
      private$parents <- parents
      private$children <- children
      private$breakpoint <- breakpoint
      private$carrier.id <- if (is.na(carrier.id)) name else carrier.id
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
    get.carrier.id = function() { private$carrier.id },
    set.carrier.id = function(id) { private$carrier.id <- id },
    
    # mutables
    get.start.time = function() { private$start.time },
    set.start.time = function(time) { private$start.time <- time },
    get.parents = function() { private$parents },
    get.parent = function() { private$parents[[1]] },
    set.parent = function(parent) {
      if (is.list(parent)) {
        if (!all(sapply(parent, function(p) is.R6(p) && is.element("Pathogen", class(p))))) {
          stop("set.parent: all elements of `parent` list must be Pathogen objects")
        }
        private$parents <- parent
      } else if (is.R6(parent) && is.element("Pathogen", class(parent))) {
        private$parents <- list(parent)
      } else {
        stop("set.parent: `parent` must be a Pathogen object or list of Pathogen objects, ",
             "not ", class(parent)[1], ". To set a parent by name, look up the Pathogen ",
             "object first (e.g. via InnerTree$get.pathogen()).")
      }
    },
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
    breakpoint = NULL,
    carrier.id = NULL
  )
)
