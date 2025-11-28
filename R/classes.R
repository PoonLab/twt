#' Deme
#' 
#' \code{Deme} is an R6 class that defines a population of host individuals
#' that share the same set of parameters, e.g., a transmission 
#' rate to other members of the deme.  Thus we can also think of this as 
#' "HostType".  A deme may represent a transmission risk group, or a regional 
#' subpopulation.  Members of a deme are either susceptible or infected.
#' 
#' @param name:  character, uniquely identifies the class
#' @param unsampled:  logical, if TRUE then no hosts in this deme can contain 
#'        sampled Pathogens (directly observed tips of the inner tree).
#' @param migration.rates:  named vector, rate that a host migrates from this 
#'        Deme to another Deme.
#' @param transmission.rates:  named vector, rates that an infected host in 
#'        this Deme transmits their infection to an *uninfected* member of 
#'        another Deme.
#' @param superinfection.rates:  named vector, rates that an infected host in 
#'        this Deme transmits their infection to an *infected* member of 
#'        another Deme.
#' @param bottleneck.size:  integer, the maximum number of Pathogen lineages 
#'        that can be transmitted to host in this Deme.
#' @param bottleneck.theta:  float, (optional) theta/size parameter of negative 
#'        binomial distribution of bottleneck size.  Converges to Poisson 
#'        distribution with larger values.  Defaults to 0 for constant 
#'        bottleneck size.
#' @param effective.size:  float, reciprocal of the rate at which Pathogen 
#'        lineages coalesce within a Host of this Deme.
#' @param popn.growth.dynamics:  character, an R expression for population 
#'        growth dynamics in forward time. If not NULL, can override 
#'        `effective.size`.
#' @param generation.time:  float, scales coalescent events to the time scale 
#'        of outer events (including transmission).
#' @param transmission.times:  numeric,  vector of transmission event times, 
#'        populated by `outer.tree.sim` from class parameters.
#' 
#' @examples 
#' 
#' # load Demes from a YAML object
#' path <- system.file('extdata', 'SI.yaml', package='twt')
#' settings <- yaml.load_file(path)
#' mod <- Model$new(settings)
#' mod$get.types()
#' 
#' # manually specify a Deme object (usually done by YAML)
#' deme <- Deme$new(name='host', transmission.rates=c(host=0.1), 
#'                  bottleneck.size=1, coalescent.rate=1)
#' 
#' @export
Deme  <- R6Class(
  "Deme",
  public = list(
    initialize = function(name=NA, unsampled = NA, migration.rates=NA, 
                          transmission.rates=NA, superinfection.rates=NA, 
                          bottleneck.size=NA, bottleneck.theta=NA,
                          effective.size=NA, popn.growth.dynamics=NA, 
                          generation.time=NA, transmission.times=NA) {
      private$name <- name
      private$unsampled <- unsampled
      
      # named vector of transmission rates corresponding to different Host 
      # objects
      private$migration.rates <- migration.rates
      private$transmission.rates <- transmission.rates
      private$superinfection.rates <- superinfection.rates
      private$bottleneck.size <- bottleneck.size
      private$bottleneck.theta <- bottleneck.theta
      
      # named vector of superinfection rates of different Hosts
      private$effective.size <- effective.size
      private$popn.growth.dynamics <- popn.growth.dynamics
      private$generation.time <- generation.time
      
      # populated after outer.tree.sim, tracked used and unused for Pathogen
      # lineage superinfection events in inner.tree.sim
      private$transmission.times <- transmission.times
    },
    
    # accessor functions
    get.bottleneck.size = function() {
      private$bottleneck.size
    },
    get.bottleneck.theta = function() {
      private$bottleneck.theta
    },
    
    get.name = function() {
      private$name
    },
    
    get.unsampled = function() {
      private$unsampled
    },
    
    get.migration.rates = function() {
      private$migration.rates
    },
    
    get.migration.rate = function(name.type) {
      private$migration.rates[[name.type]]
    },
    
    get.transmission.rates = function() {
      private$transmission.rates
    },
    
    get.transmission.rate = function(current.time, name.type) {
      # FIXME: current.time is in what time frame?
      if (length(private$transmission.rates) == 1) {
        # constant rate over time
        private$transmission.rates[[1]][[name.type]]
      }
      else {
        # rate variation over time
        rate.changes <- as.numeric(names(private$transmission.rates))
        index <- max(which(rate.changes >= current.time))
        private$transmission.rates[[index]][[name.type]]
      }
    },
    
    get.superinfection.rates = function() {
      private$superinfection.rates
    },
    
    get.superinfection.rate = function(name.type) {
      # FIXME: We might want to have time-varying superinfection rates too
      private$superinfection.rates[[name.type]]
    },
    
    set.superinfection.rate = function(recipient.type, new.superinfection.rate) {
      private$superinfection.rates[[recipient.type]] <- new.superinfection.rate
    }
    
    get.effective.size = function() {
      private$effective.size
    },
    
    get.popn.growth.dynamics = function() {
      private$popn.growth.dynamics
    },
    
    get.generation.time = function() {
      private$generation.time
    },
    
    get.transmission.times = function() {
      private$transmission.times
    },
    
    set.transmission.times = function(vector.transm.times) {
      private$transmission.times <- vector.transm.times
    },
    
  ),
  private = list(
    name = NULL,
    unsampled = NULL,
    
    transmission.rates = NULL,
    migration.rates = NULL,
    superinfection.rates = NULL,
    
    bottleneck.size = NULL,
    bottleneck.theta = NULL,
    
    effective.size = NULL,
    popn.growth.dynamics = NULL,
    generation.time = NULL,
    transmission.times = NULL
  )
)


print.Deme <- function(obj) {
  cat(paste(obj$name, ":"))
}


#' Host
#' 
#' \code{Host} is an R6 class for objects that represent the units of an
#' outer tree simulation, i.e., a host individual.
#' 
#' @param name:  character, unique identifier for Host
#' @param deme:  character, reference Deme by name
#' @param source:  character, reference to another Host from which a Pathogen 
#'                 lineage was first transmitted to this Host
#'                 TODO: we will need to track multiple sources for a 
#'                 transmission history (graph) that tracks superinfection
#' @param transmission.time:  numeric, stores the origin time of this Host, 
#'                            which corresponds to a transmission event in the 
#'                            "outer" tree.
#' @param unsampled:  logical, if TRUE, then any Pathogen carried by this Host 
#'                    is not directly observed, i.e., it does not represent a 
#'                    tip in the "inner" tree.
#' 
#' @examples 
#' # load Host objects from a YAML file
#' path <- system.file('extdata', 'SI.yaml', package='twt')
#' settings <- yaml.load_file(path)
#' mod <- MODEL$new(settings)
#' 
#' # display first Host object
#' host1 <- mod$get.hosts()[[1]]
#' host1
#' 
#' # manually initialize a new Host object
#' hostN <- Host$new(name='newHost', unsampled=TRUE)
#' hostN$set.deme(host1$get.deme())
#' 
#' @export 
Host <- R6Class("Host",
  public = list(
    initialize = function(name=NA, deme=NA, source=NA, transmission.time=NA, 
                          sampling.time=NA, unsampled=FALSE, pathogens=list()) {
      private$name <- name
      private$type <- type
      private$source <- source
      private$transmission.time <- transmission.time
      private$sampling.time <- sampling.time
      # attr req later when identifying new US Hosts to be promoted in superinfection events
      private$unsampled <- unsampled
      private$pathogens <- pathogens
    },
    
    copy = function(deep=FALSE) {
      # see https://github.com/r-lib/R6/issues/110
      cloned <- self$clone(deep)  # calls deep_clone method
      if (deep) {
        # attach new Pathogens to new Host
        for (pathogen in cloned$get.pathogens()) {
          pathogen$set.location(cloned)
        }
      }
      cloned  # return
    },
    
    # accessor functions
    get.name = function() {
      private$name
    },
    
    get.type = function() {
      private$type
    },
    
    get.source = function() {
      private$source
    },
    
    get.transmission.times = function() {
      private$transmission.time
    },
    
    get.sampling.time = function() {
      private$sampling.time
    },
    
    set.type = function(new.type) {
      private$type <- new.type
    },
    
    set.source = function(new.source) {
      private$source <- new.source
    },
    
    set.transmission.times = function(new.transmission.time) {
      private$transmission.time <- new.transmission.time
    },
    
    set.sampling.time = function() {
      private$sampling.time <- min(sapply(private$pathogens, function(line) {
        line$get.sampling.time()
      }))
    },
    
    set.unsampled = function(is.unsampled) {
      private$unsampled <- is.unsampled
    },
    
    is.unsampled = function() {
      private$unsampled
    },
    
    get.pathogens = function() {
      private$pathogens
    },
    
    add.pathogen = function(new.pathogen) {
      private$pathogens[[length(private$pathogens)+1]] <- new.pathogen
      self$set.sampling.time()
    },
    
    remove.pathogen = function(ex.pathogen) {
      private$pathogens[[ex.pathogen$get.name()]] <- NULL
      self$set.sampling.time()
    }
    
  ),
  private = list(
    name = NULL,
    type = NULL,  # reference to Deme object
    source = NULL,
    transmission.time = NULL,
    sampling.time = NULL,
    unsampled = NULL,
    pathogens = NULL,
    
    deep_clone = function(name, value) {
      if (name == 'pathogens') {
        # map deep clone to Pathogen copy() method
        lapply(value, function(pathogen) pathogen$copy(deep=TRUE))
      } 
      else {
        value
      }
    }
  )
)




#' Pathogen
#' 
#' \code{Pathogen} is an R6 class for objects that represent pathogen lineages
#' that are carried by Hosts and which comprise the "inner" tree of the 
#' simulation.
#' 
#' @param name: character, unique identifier for Pathogen object
#' @param type: character, reference to object of class PathogenType (not yet 
#'              implemented)
#' @param sampling.time: numeric, time that the Pathogen was sampled; left to 
#'                       NA for unsampled Pathogens
#' @param location: a reference to a Host object
#' 
#' @examples
#' # load model from a YAML file
#' path <- system.file('extdata', 'SI.yaml', package='twt')
#' settings <- yaml.load_file(path)
#' mod <- MODEL$new(settings)
#' 
#' # display Pathogens in first Host object
#' host <- mod$get.hosts()[[1]]
#' host$get.pathogens()  # display all 3 Pathogen objects
#' 
#' # manually add an unsampled Pathogen
#' pgn <- Pathogen$new(name="L0", location=comp)
#' host$add.pathogen(pgn)
#' 
#' 
#' @export
Pathogen <- R6Class("Pathogen",
  public = list(
    initialize = function(name=NA, type=NA, sampling.time=NA, location=NA) {
      private$name <- name
      private$type <- type
      private$sampling.time <- sampling.time
      private$location <- location
    },
    
    parent = NULL,
   
    copy = function(deep=FALSE) {
      # see https://github.com/r-lib/R6/issues/110
      if (deep) {
        parent <- private$location
        private$location <- NULL  # temporarily erase before cloning!
      }
     
      cloned <- self$clone(deep)
     
      if (deep) {
        private$location <- parent  # restore original reference
      }
     
      cloned
    },
   
    get.name = function() {
      private$name
    },
   
    get.type = function() {  # in the future, will be a pointer to a PathogenType object
      private$type
    },
   
    get.sampling.time = function() {
      private$sampling.time
    },
   
    get.location = function() {
      private$location
    },
     
    set.location = function(comp) {
      private$location <- comp
    },
   
    set.location.by.name = function(locationList, new.locationName) {
      new.locationObj <- locationList[[ 
        which(sapply(locationList, function(x) {x$get.name()}) == new.locationName) 
        ]]
      private$location <- new.locationObj
    }
    
  ),
  private = list(
    name = NULL,
    type = NULL,            # potential reference to PathogenType object
    sampling.time = NULL,
    location = NULL
  )
)


