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
#' @param sampling.comp:  character, name of sampling compartment.  Used for 
#'        checking whether target sample sizes are obtained.
#' @param pathogens:  list, store objects of class Pathogen associated with 
#'        this Host
#' @export
Host <- R6Class(
  "Host",
  public = list(
    initialize = function(
    name=NA, compartment=NA, source=character(), transmission.time=numeric(), 
    sampling.time=as.numeric(NA), sampling.comp=NA, unsampled=FALSE, 
    pathogens=list()) {
      private$name <- name
      private$compartment <- compartment
      private$source <- source
      private$transmission.time <- transmission.time
      private$sampling.time <- sampling.time
      private$sampling.comp <- sampling.comp
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
    
    get.sampling.comp = function() { private$sampling.comp },
    set.sampling.comp = function(comp) { private$sampling.comp <- comp },
    
    get.pathogens = function() { private$pathogens },
    count.pathogens = function () { length(private$pathogens) },
    add.pathogen = function(new.pathogen) {
      private$pathogens[[length(private$pathogens)+1]] <- new.pathogen
      #self$set.sampling.time()
    },
    remove.pathogen = function(ex.pathogen) {
      private$pathogens[[ex.pathogen$get.name()]] <- NULL
      #self$set.sampling.time()
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
    sampling.comp = NULL,
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
    initialize = function(name=NA, hosts=list(), index=1) {
      private$name <- name
      private$hosts <- hosts
      private$index <- index  # unique IDs for members
    },
    
    get.names = function() {
      sapply(private$hosts, function(h) { h$get.name() })
    },
    
    get.hosts = function() {
      private$hosts
    },
    
    get.host.by.name = function(host.name, remove=FALSE) {
      members <- self$get.names()
      if (host.name %in% members) {
        idx <- which(members == host.name)
        if (remove) {
          self$remove.host(idx)
        } else {
          private$hosts[[idx]]
        }
      } else {
        warning("HostSet does not contain Host ", host.name)
        invisible(NULL)
      }
    },
    
    get.types = function() {
      # generate a vector of compartment memberships for current Hosts
      sapply(private$hosts, function(h) { h$get.compartment() })
    },
    
    get.sampling.times = function() {
      sapply(private$hosts, function(h) { h$get.sampling.time() })
    },
    
    get.sampling.types = function() {
      sapply(private$hosts, function(h) { h$get.sampling.comp() })
    },
    
    get.transmission.times = function() {
      lapply(private$hosts, function(h) { h$get.transmission.time() })
    },
    
    is.sampled = function() {
      sapply(private$hosts, function(h) { h$is.sampled() })
    },
    
    get.sources = function() {
      lapply(private$hosts, function(h) { h$get.source() })
    },
    
    count.type = function(type=NA) {
      host.types <- self$get.types()
      if (is.na(type)) {
        length(host.types)
      } else {
        sum(host.types == type)
      }
    },
    
    count.sampling.type = function(type) {
      sampling.types <- self$get.sampling.types()
      sum(sampling.types == type)
    },
    
    add.host = function(host) {
      host.name <- host$get.name()
      if (is.na(host.name)) {
        # assign a unique name
        if (host$is.sampled()) {
          name <- paste(host$get.compartment(), private$index, sep="_")
        } else {
          name <- paste("US", host$get.compartment(), private$index, sep="_")
        }
        host$set.name(name)        
      } else {
        # check if named Host is already in this HostSet
        current <- self$get.names()
        if (host.name %in% current) {
          warning("Host ", host.name, " is already in this HostSet!")
          invisible(NULL)
        }
      }
      private$index <- private$index + 1
      private$hosts[[length(private$hosts)+1]] <- host
    },
    
    remove.host = function(idx) {
      if (idx < 1) { 
        stop("HostSet$remove.host cannot use negative index")
      } else if (idx > self$count.type()) {
        stop("HostSet$remove.host index ", idx, " exceeds number of hosts")
      }
      host <- private$hosts[[idx]]
      private$hosts <- private$hosts[-idx]
      return(host)
    },
    
    sample.host = function(type=NA, remove=FALSE) {
      if (is.na(type)) {
        # any host will do
        idx <- sample(1:self$count.type(), 1)
        if (remove) {
          self$remove.host(idx)  # without replacement
        } else {
          private$hosts[[idx]]  # with replacement
        }
      } else {
        # sample host of specific type
        if (self$count.type(type) > 0) {
          idx <- sample(which(self$get.types()==type), 1)
          if (remove) {
            self$remove.host(idx)  # without replacement
          } else {
            # with replacement
            private$hosts[[idx]]
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
