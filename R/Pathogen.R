#' Pathogen
#' 
#' \code{Pathogen} is an R6 class for objects that represent pathogen lineages
#' that are carried by Hosts and which comprise the "inner" tree of the 
#' simulation.
#' 
#' @param name: character, unique identifier for Pathogen object
#' @param sampling.time: numeric, time that the Pathogen was sampled; left to 
#'                       NA for unsampled Pathogens
#' @param host: a reference to a Host object
#' 
#' @export
Pathogen <- R6Class(
  "Pathogen",
  public = list(
    initialize = function(name=NA, type=NA, sampling.time=NA, location=NA) {
      private$name <- name
      private$sampling.time <- sampling.time
      private$host <- host
    },
    
    parent = NULL,
    
    copy = function(deep=FALSE) {
      # see https://github.com/r-lib/R6/issues/110
      if (deep) {
        parent <- private$host
        private$host <- NULL  # temporarily erase before cloning!
      }
      cloned <- self$clone(deep)
      if (deep) {
        private$location <- parent  # restore original reference
      }
      
      cloned
    },
    
    get.name = function() { private$name },
    get.sampling.time = function() { private$sampling.time },
    
    get.host = function() { private$host },
    set.host = function(host) {
      private$host <- host
    },
    
    set.host.by.name = function(hostList, new.hostName) {
      # FIXME: unused, deprecated?
      new.hostObj <- hostList[[ 
        which(sapply(hostList, function(x) {x$get.name()}) == new.hostName) 
      ]]
      private$host <- new.hostObj
    }
  ),
  private = list(
    name = NULL,
    sampling.time = NULL,
    host = NULL
  )
)
