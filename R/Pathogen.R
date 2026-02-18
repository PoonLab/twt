#' Pathogen
#' 
#' \code{Pathogen} is an R6 class for objects that represent pathogen lineages
#' that are carried by Hosts and which comprise the "inner" tree of the 
#' simulation.
#' 
#' @param name: character, unique identifier for Pathogen object
#' @param sampling.time: numeric, time that the Pathogen was sampled; left to 
#'                       NA for unsampled Pathogens
#' @param host:  character, name of Host object associated with this Pathogen
#' 
#' @export
Pathogen <- R6Class(
  "Pathogen",
  public = list(
    initialize = function(name=NA, sampling.time=NA, host=NA) {
      private$name <- name
      private$sampling.time <- sampling.time
      private$host <- host
    },
    
    # immutable attributes
    get.name = function() { private$name },
    get.sampling.time = function() { private$sampling.time },
    
    get.host = function() { private$host },
    set.host = function(host) {
      private$host <- host
    }
  ),
  private = list(
    name = NULL,
    sampling.time = NULL,
    host = NULL
  )
)
