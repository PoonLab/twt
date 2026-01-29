#' Host
#' 
#' \code{Host} is an R6 class representing a host individual carrying 
#' one or more sampled Pathogen lineages.  The transmission history among
#' hosts comprise the "outer" tree of the simulation.
#' 
#' @param name:  character, unique identifier for the host
#' @param compartment:  character, reference to Compartment the host currently
#'        belongs to.
#' @param source:  character, reference to another Host from which a Pathogen 
#'        lineage was transmitted to this Host.  May store more than one source
#'        in the case of superinfection.
#' @param transmission.time:  numeric, stores the time(s) that a Pathogen 
#'        lineage weas transmitted to this Host.  Uses a named vector to 
#'        make the association of time to source Host less ambiguous.
#' @param sampling.time:  numeric, the earliest time that a Pathogen lineage 
#'        was directly sampled from this Host.  If NA, then this is an unsampled
#'        Host, i.e., its lineages are ancestral to lineages in subsequent 
#'        hosts where they are sampled.
#' @param pathogens:  list, store objects of class Pathogen associated with 
#'        this Host
#' @export
Host <- R6Class(
  "Host",
  public = list(
    initialize = function(name=NA, compartment=NA, source=character(), 
                          transmission.time=numeric(), sampling.time=NA, 
                          unsampled=FALSE, pathogens=list()) {
      private$name <- name
      private$compartment <- compartment
      private$source <- source
      private$transmission.time <- transmission.time
      private$sampling.time <- sampling.time
      private$unsampled <- unsampled
      private$pathogens <- pathogens
    },
    
    # accessor functions
    get.name = function() { private$name },
    set.name = function(name) { private$name = name },
    
    get.compartment = function() { private$compartment },
    set.compartment = function(comp) { private$compartment=comp },
    
    get.source = function() { private$source },
    set.source = function(new.source) {
      if (length(private$source) > 0) {
        # superinfection
        private$source <- c(private$source, new.source)
      } else {
        private$source <- new.source
      }
    },
    
    get.transmission.time = function() { private$transmission.time },
    set.transmission.time = function(new.time) {
      if (length(private$transmission.time) > 0) {
        # superinfection
        private$transmission.time <- c(private$transmission.time, new.time)
      } else {
        private$transmission.time <- new.time  
      }
    },
    
    get.sampling.time = function() { private$sampling.time },
    set.sampling.time = function() {
      private$sampling.time <- min(sapply(private$pathogens, function(line) {
        line$get.sampling.time()
      }))
    },
    
    get.pathogens = function() { private$pathogens },
    add.pathogen = function(new.pathogen) {
      private$pathogens[[length(private$pathogens)+1]] <- new.pathogen
      self$set.sampling.time()
    },
    remove.pathogen = function(ex.pathogen) {
      private$pathogens[[ex.pathogen$get.name()]] <- NULL
      self$set.sampling.time()
    },
    
    is.sampled = function() {
      #!all(is.na(private$sampling.time))
      !private$unsampled
    }
  ),
  
  private = list(
    name = NULL,
    compartment = NULL,
    source = NULL,
    transmission.time = NULL,
    sampling.time = NULL,
    unsampled = NULL,
    pathogens = NULL
  )
)


#' HostSet
#' 
#' A mutable set of Host objects, used to make some methods more convenient
#' when simulating the outer tree.
#' @param name:  character, optional identifier
#' @param hosts:  list, used to build a HostSet from existing Host objects
#' @param index:  integer, internal counter
#' @export
HostSet <- R6Class(
  "HostSet",
  public = list(
    initialize = function(name=NA, hosts=list(), index=0) {
      private$name <- name
      private$hosts <- hosts
      private$index <- 0  # unique IDs for members
    },
    
    get.types = function() {
      # generate a vector of compartment memberships for current Hosts
      sapply(private$hosts, function(h) { h$get.compartment() })
    },
    
    count.type = function(type=NA) {
      host.types <- private$get.types()
      if (is.na(type)) {
        length(host.types)
      } else {
        sum(host.types == type)
      }
    },
    
    add.host = function(host) {
      # assign a unique name
      host$set.name(paste(
        ifelse(host$is.sampled(), "", "US"), 
        host$get.compartment(), private$index, 
        sep="_"))
      private$index <- private$index + 1
      private$hosts[[length(private$hosts)+1]] <- host
    },
    
    remove.host = function(idx) {
      if (idx < 1) { 
        stop("HostSet$remove.host cannot use negative index")
      } else if (idx > private$count.type()) {
        stop("HostSet$remove.host index ", idx, " exceeds number of hosts")
      }
      host <- private$hosts[idx]
      private$hosts <- private$hosts[-idx]
      return(host)
    },
    
    sample.host = function(type=NA, remove=FALSE) {
      if (is.na(type)) {
        # any host will do
        idx <- sample(1:length(private$count.type()), 1)
        if (remove) {
          private$remove.host(idx)
        } else {
          private$hosts[[idx]]
        }
      } else {
        if (private$count.type(type) > 0) {
          if (remove) {
            idx <- sample(which(private$get.types()==type), 1)
            private$remove.host(idx)
          } else {
            sample(unlist(private$hosts[private$get.types()==type]), 1)
          }
        } else {
          stop("HostSet$sample.host cannot return Host of type ", type)
        }
      }
    }
  ),
  
  private = list(
    name = NULL,
    hosts = NULL,
    index = NULL
  )
)


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

