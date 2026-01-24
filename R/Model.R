#' Model
#'
#' \code{Model} is an R6 class that defines an object that generates all of the 
#' objects of the various classes (e.g., Compartment, Host, Pathogen) that 
#' define a simulation model.  A \code{Model} object is immutable - it should
#' not change over the course of simulation.  Instead, we derive a \code{Run}
#' object class that inherits from a \code{Model}.
#' 
#' @param settings: a named list returned by `yaml.load_file()` that contains 
#' user specifications of the simulation model.
#' 
#' @field initial.conds  Initial Conditions
#' @field compartments  vector of Compartment objects
#' @field hosts  vector of Host objects
#' @field lineages  vector of Lineage objects
#' @field fixed.samplings  list of Lineage names and sampling times for plotting
#' 
#' @export
Model <- R6Class("Model",
  #lock_objects = FALSE,  # FIXME: is this deprecated?
  
  public = list(
    initialize = function(settings=NA) {
      private$load.parameters(settings)
      private$load.compartments(settings)
      private$load.sampling(settings)
    },
    
    # ACCESSOR FUNCTIONS
    get.parameters = function() { private$parameters },
    get.compartments = function() { private$compartments },
    
    #' returns names of a given list of R6 objects
    #' @param listR6obj: list of R6 objects of class Compartment, Host or 
    #' Lineage
    get.names = function(listR6obj) {
      unname(sapply(listR6obj, function(x){ x$get.name() }))
    }
  ),
  
  private = list(
    parameters = NULL,
    compartments = NULL,
    sampling = NULL,
      
    load.parameters = function(settings) {
      if (is.null(settings$Parameters)) {
        stop("Error loading Model, Parameters key missing from settings")
      }
      
      # origin time is measured in reverse (prior to most recent sampled lineage)
      params <- settings$Parameters
      if (is.null(params$originTime)) {
        stop("Parameters: required key `originTime` is missing")
      }
      if (!is.numeric(params$originTime)) {
        stop("InitialConditions: originTime must be numeric")
      }
      if (params$originTime <= 0) {
        stop("InitialConditions: originTime must be positive")
      }
      private$parameters <- params
    },
    
    check.expression = function(s, env) {
      # is string `s` a valid R expression?
      tryCatch({
        x <- parse(text=s)
      }, error = function(e) {
        stop(paste("Invalid R expression `", x, "`"))
      })
      # are all parameters and variables declared?
      tryCatch({
        v <- eval(x, envir=env)
      }, error = function(e) {
        stop(paste(
          "Expression `", s, "` contains one or more undeclared variables"
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
      cnames <- names(settings$Compartments)
      for (cn in cnames) {
        eval(parse(text=paste(cn, "<-", 1)), envir=env)
      }
      cat(eval(parse(text="ls()"), envir=env))
      
      private$compartments <- sapply(cnames, function(cn) {
        params <- settings$Compartments[[cn]]
        if (is.null(params$size)) {
          stop("Compartment `", cn, "` must specify initial `size`")
        }
        
        if (!is.null(params$rates)) {
          # check rate specifications
          for (dest in names(params$rates)) {
            rate <- params$rates[[dest]]  # from `cn` to `dest`
            if (is.numeric(rate)) {
              pass  # constant rate
            } else if (is.character(rate)) {
              private$check.expression(rate, env)
            } else {
              stop("Unexpected type in load.compartments")
            }
          }          
        }
        
        # check bottleneck setting
        bottleneck.size <- NA
        if (!is.null(params$bottleneck.size)) {
          bottleneck.size <- params$bottleneck.size
          if (!is.numeric(bottleneck.size)) {
            private$check.expression(bottleneck.size, env)
          }
        }
        
        # check coalescent rate
        coalescent.rate <- NA
        if (!is.null(params$coalescent.rate)) {
          coalescent.rate <- params$coalescent.rate
          if (!is.numeric(coalescent.rate)) {
            private$check.expression(coalescent.rate, env)
          }
        }
        
        Compartment$new(
          name = cn,
          rates = params$rates,
          size = params$size,
          bottleneck.size = bottleneck.size,
          coalescent.rate = coalescent.rate
        )  # return value
      })
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
          if (!is.element(cn, names(private$compartments))) {
            stop(paste("In Sampling:targets compartment", cn, 
                       "has not been declared in Compartments"))
          }
          count <- private$sampling$targets[[cn]]
          if (!is.integer(count) | count <= 0) {
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
  for (cn in names(compartments)) {
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
plot.Model <- function(obj) {
  # TODO: work in progress!
}
