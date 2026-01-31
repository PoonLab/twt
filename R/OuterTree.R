#' OuterTree
#' 
#' \code{OuterTree} is an R6 class defining a persistent object that records 
#' events while simulating the outer tree from population trajectories.
#' 
#' @param mod:  R6 object of class Model
#' @export
OuterTree <- R6Class(
  "OuterTree",
  public = list(
    initialize = function(mod) {
      private$outer.log <- data.frame(
        time=numeric(), event=character(), 
        host.anc=character(), host.des=character(), 
        comp.anc=character(), comp.des=character()
      )
      
      private$targets <- mod$sampling$targets
      
      private$sampled <- setNames(
        rep(0, length(private$targets)), 
        names(private$targets))
      
      private$infected <- mod$get.infected()
    },
    
    # accessor functions
    get.targets = function() { private$targets }, 
    get.log = function() { private$outer.log },
    get.nrow = function() { nrow(private$outer.log) },
    get.sampled = function() { private$sampled },
    get.infected = function() { private$infected },
    nsamples = function() { sum(private$sampled) },
    
    add.event = function(e) {
      if ( !all(names(e) == names(private$outer.log)) ) {
        stop("Event passed to OuterTree:add.event must be vector with ",
             "names: time, event, host.anc, host.des, comp.anc, comp.des")
      }
      private$outer.log[self$get.nrow()+1] <- e
    },
    
    add.sample = function(cn) {
      if ( !is.element(cn, names(private$sampled)) ) {
        stop("Error in OuterTree:add.sample() Unrecognized compartment ",
             "name `", cn, "`")
      }
      private$sampled[[cn]] <- private$sampled[[cn]] + 1
    }
  ),
  
  private = list(
    outer.log = NULL,
    targets = NULL,
    sampled = NULL,
    infected = NULL
  )
)