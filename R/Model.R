#' Model
#'
#' \code{Model} is an R6 class that defines an object that parses the model 
#' settings from a YAML (Parameters, Compartments, Sampling).
#' 
#' Four types of events affect Compartments:
#'   births         0->S      dS/dt = a [or] aS
#'   deaths         S->0      dS/dt = -a [or] aS
#'   migrations     I->S      dI/dt = -aI, dS/dt = aI
#'   transmissions  S+I->I+I  dS/dt = -bSI, dI/dt = bSI
#' In the case of superinfection, transmission occurs within a Compartment
#'   I + I -> I + I
#' This does not affect Compartment size and dynamics, but it is a branching
#' event.
#' 
#' It is important that individual model terms are separated!  For example,
#' there may be birth and death rates associated with a Compartment:
#'   dS/dt = lambda S - delta S [YES]
#'         = (lambda-delta)S    [NO!]
#' In other words, we want to model these stochastic events separately 
#' instead of the net rate of change.
#' 
#' @param settings: a named list returned by `yaml.load_file()` that contains 
#' user specifications of the simulation model.
#' 
#' @field parameters  model parameters
#' @field compartments  vector of Compartment objects
#' @field sampling  sampling conditions
#' 
#' @examples
#' require(R6)
#' require(yaml)
#' require(igraph)
#' settings <- yaml.load_file("examples/SIRS_serial.yaml")
#' mod <- Model$new(settings)
#' summary(mod)
#' plot(mod)
#' 
#' @export
Model <- R6Class(
  "Model",
  public = list(
    initialize = function(settings=NA) {
      private$load.parameters(settings)
      private$load.compartments(settings)
      private$load.sampling(settings)
    },
    
    # ACCESSOR FUNCTIONS (immutable object, no set methods)
    get.parameters = function() { private$parameters },
    get.compartments = function() { private$compartments },
    get.sampling = function() { private$sampling },
    
    # these return named vectors
    get.init.sizes = function() { private$init.sizes },
    get.infected = function(cn=NA) { 
      if (is.na(cn)) { 
        private$is.infected  # return entire list
      } else {
        if( is.element(cn, names(private$is.infected)) ) {
          private$is.infected[[cn]]  # retrieve one entry from list
        } else {
          NA  # failed to find key
        }
      }
    },
    get.birth.rates = function() { private$birth.rates },
    get.death.rates = function() { private$death.rates },
    
    # migration rates are stored in a K x K matrix
    get.migration.rates = function() { private$migration.rates },
    
    # transmission rates are in a K x K x K matrix (from comp, to comp, source)
    get.transmission.rates = function() { private$transmission.rates },
    
    # named vectors
    get.bottleneck.size = function(cname=NA) { 
      if (is.na(cname)) {
        return(private$bottleneck.sizes)  # return all sizes as named vector
        
      } else {
        if (cname %in% self$get.compartments()) {
          return (private$bottleneck.sizes[[cname]])
        }
        warning("Unrecognized compartment name ", cname)
        return (NULL)
      }
    },
    
    get.coalescent.rate = function(cname=NA) { 
      if (is.na(cname)) {
        return (private$coalescent.rates)  # return all rates as named vector
        
      } else {
        if (cname %in% self$get.compartments()) {
          return (private$coalescent.rates[[cname]])
        }
        warning("Unrecognized compartment name ", cname)
        return (NULL)
      }
    },
    
    get.graph = function() { private$graph }
  ),
  
  private = list(
    parameters = NULL,
    compartments = NULL,  # character, compartment names
    sampling = NULL,
    
    init.sizes = NULL,  # numeric, initial sizes
    is.infected = NULL,  # boolean, does compartment carry pathogens?

    birth.rates = NULL,  # list, birth rates per compartment
    death.rates = NULL,  # list, death rates per compartment
    migration.rates = NULL,  # 
    transmission.rates = NULL,
    
    bottleneck.sizes = NULL,  # note these are expressions, not numeric
    coalescent.rates = NULL,
    
    graph = NULL,
    
    load.parameters = function(settings) {
      if (is.null(settings$Parameters)) {
        stop("Error loading Model, Parameters key missing from settings")
      }
      
      # origin time is measured in reverse (prior to most recent sampled lineage)
      params <- settings$Parameters
      if (is.null(params$simTime)) {
        stop("Parameters: required key `simTime` is missing")
      }
      if (!is.numeric(params$simTime)) {
        stop("InitialConditions: simTime must be numeric")
      }
      if (params$simTime <= 0) {
        stop("InitialConditions: simTime must be positive")
      }
      private$parameters <- params
    },
    
    check.expression = function(s, env) {
      # is string `s` a valid R expression?
      tryCatch({
        x <- parse(text=s)
      }, error = function(e) {
        stop(paste("Invalid R expression `", s, "`"))
      })
      # are all parameters and variables declared?
      tryCatch({
        v <- eval(x, envir=env)
      }, error = function(e) {
        stop(paste(
          "Expression `", s, "` contains one or more undeclared variables\n",
          eval(parse(text="ls()"), envir = env)
        ))
      })
      if (!is.numeric(v)) {
        stop(paste(
          "Expression `", s, "` does not evaluate to a numeric value"
        ))
      }
    },
    
    # Parse Compartments from settings (YAML)
    load.compartments = function(settings) {
      if (is.null(settings$Compartments)) {
        stop("Missing 'Compartments' field in settings")
      }
      if (length(settings$Compartments) == 0) {
        stop("Empty 'Compartments' field in settings.")
      }
      
      # if parameters have not been validated, do it now
      if (is.null(private$parameters)) {
        load.parameters(settings)
      }
      
      # declare parameters and compartments in a new environment
      env <- new.env()
      for (key in names(private$parameters)) {
        eval(parse(text=paste(key, "<-", private$parameters[[key]])), envir=env)
      }
      private$compartments <- cnames <- names(settings$Compartments)
      k <- length(cnames)  # number of compartments
      for (cn in cnames) {
        eval(parse(text=paste(cn, "<-", 1)), envir=env)
      }
      # cat(eval(parse(text="ls()"), envir=env))  ## DEBUGGING
      
      # initialize containers
      private$init.sizes <- setNames(rep(0, k), cnames)
      private$is.infected <- setNames(rep(NA, k), cnames)
      private$bottleneck.sizes <- setNames(rep("1", k), cnames)
      private$coalescent.rates <- setNames(rep("Inf", k), cnames)
      
      private$birth.rates <- setNames(rep("0", k), cnames)
      private$death.rates <- setNames(rep("0", k), cnames)
      
      # 2D matrix (from compartment, to compartment)
      private$migration.rates <- matrix(
        "0", nrow=k, ncol=k, dimnames=list(cnames, cnames))
      
      # 3D matrix
      private$transmission.rates <- array( # from, to, infected by
        "0", dim=c(k, k, k), dimnames=list(cnames, cnames, cnames))
      
      for (src in cnames) {
        params <- settings$Compartments[[src]]
        
        # initial population size
        if (is.null(params$size)) {
          stop("Compartment `", src, "` must specify initial `size`")
        }
        private$init.sizes[[src]] <- params$size
        
        # does compartment carry pathogens?
        if (is.null(params$infected)) {
          stop("Compartment `", src, "` must specify `infected` (true/false)")
        }
        private$is.infected[[src]] <- params$infected
        
        # parse rate expressions
        if (!is.null(params$birth)) {
          if (is.list(params$birth)) {
            stop(src, "`birth` should be a number or R expression, not an",
                 "associative array.")
          }
          private$check.expression(params$birth, env)
          private$birth.rates[[src]] <- params$birth
        }
        
        if (!is.null(params$death)) {
          if (is.list(params$death)) {
            stop(src, "`death` should be a number or R expression, not an",
                 "associative array.")
          }
          private$check.expression(params$death, env)
          private$death.rates[[src]] <- params$death
        }
        
        if (!is.null(params$migration)) {
          if (!is.list(params$migration)) {
            stop(src, "`migration` should be an associative array.")
          }
          dests <- names(params$migration)
          if (any(!is.element(dests, cnames))) {
            stop(src, "`migration` contains undeclared Compartment(s): ",
                 dests[which(!is.element(dests, cnames))])
          }
          for (dest in names(params$migration)) {
            rate <- params$migration[[dest]]
            private$check.expression(rate, env)
            private$migration.rates[[src, dest]] <- rate
          }
        }
        
        if (!is.null(params$transmission)) {
          if (!is.list(params$transmission)) {
            stop(src, "`transmission` should be an associative array.")
          }
          dests <- names(params$transmission)
          if (any(!is.element(dests, cnames))) {
            stop(src, "`transmission` contains undeclared Compartment(s): ",
                 dests[which(!is.element(dests, cnames))])
          }
          for (dest in dests) {
            rates <- params$transmission[[dest]]
            if (!is.list(rates)) {
              stop("`transmission` from `", src, "` to `", dest, "` must be ",
                   "an associative array.")
            }
            
            for (source in names(rates)) {
              rate <- rates[[source]]
              private$check.expression(rate, env)
              private$transmission.rates[[src, dest, source]] <- rate
            }
          }
        }
        
        # check bottleneck setting
        if (!is.null(params$bottleneck.size)) {
          if (is.list(params$bottleneck.size)) {
            stop(src, "`bottleneck.size` should be a number or R expression,",
                 "not an associative array.")
          }
          private$check.expression(params$bottleneck.size, env)
          private$bottleneck.sizes[[src]] <- params$bottleneck.size
        }
        
        # check coalescent rate
        if (!is.null(params$coalescent.rate)) {
          if (is.list(params$coalescent.rate)) {
            stop(src, "`coalescent.rate` should be a number or R expression,",
                 "not an associative array.")
          }
          private$check.expression(params$coalescent.rate, env)
          private$coalescent.rates[[src]] <- params$coalescent.rate
        }
      }  # end loop
      
      # generate graph of compartments
      squash <- apply(private$transmission.rates, c(1,2), 
                      function(x) any(x!="0"))
      adj <- (squash | private$migration.rates != "0")
      private$graph <- igraph::graph_from_adjacency_matrix(
        adj, mode="directed", diag=F)
      is.orphan <- (igraph::degree(private$graph)==0)
      if (any(is.orphan)) {
        stop("Isolated compartments detected:", private$compartments[is.orphan])
      }
    },
    
    # Validate Sampling settings
    load.sampling = function(settings) {
      if (is.null(settings$Sampling)) {
        stop("'Sampling' is a required field in settings (YAML)")
      }
      private$sampling <- settings$Sampling
      if (is.null(private$sampling$mode)) {
        stop("'mode' is a required field within 'Sampling' in settings")
      }
      mode <- private$sampling$mode
    
      if (mode == "compartment") {
        # one or more Compartments represent sampled hosts, i.e., constant 
        # sampling rate.  Populate Host objects during trajectory simulation.
        
        if (is.null(private$sampling$targets)) {
          stop("Sampling field must specify `targets` (compartments)")
        }
        
        # check that targeted compartments have been defined by user
        if (is.null(private$compartments)) {
          load.compartments(settings)
        }
        for (cn in names(private$sampling$targets)) {
          if (!is.element(cn, private$compartments)) {
            stop(paste("In Sampling:targets compartment", cn, 
                       "has not been declared in Compartments"))
          }
          count <- private$sampling$targets[[cn]]
          if (!is.numeric(count) | count <= 0) {
            stop(paste("Sampling:targets:", cn, "must be a positive integer"))
          }
        }
      } else {
        stop(paste("Sampling mode `", mode, "` not recognized"))
      }
    }
  )
)


#' print.Model
#' S3 class function to display contents of a Model object
#' @export
print.Model <- function(obj) {
  cat("twt Model")
  cat("Parameters:\n")
  params <- obj$get.parameters()
  for (key in names(params)) {
    value <- params[[key]]
    cat("  ", key, ": ", value, "\n")
  }
  cat("Compartments: ")
  compartments <- obj$get.compartments()
  for (cn in compartments) {
    cat(cn, " ")
  }
  cat("\n")
}

#' summary.Model
#' S3 class function to summarize a Model object - simply a wrapper 
#' around print()
#' @export
summary.Model <- function(obj) {
  print(obj)
}


#' plot.Model
#' Plot a graph summarizing compartments and rates
#' TODO: label edges with rate expressions
#' @export
plot.Model <- function(obj) {
  igraph::plot.igraph(obj$get.graph())
}
