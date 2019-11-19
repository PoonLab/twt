#' sim.inner.tree
#' 
#' Simulate the coalescence of Lineages within Compartments, and resolve
#' migration events that may involve sampled Lineages
#' 
#' @param model:  R6 object of class Model or Run.  If user provides a Model 
#'        object, then an EventLogger object must also be provided.  The 
#'        expected use case is generating an eventlog from `eventlog.from.tree`
#'        and using Model to parameterize the inner tree simulation.
#' @param e: (optional)  R6 object of class EventLogger
#' @return R6 object of class EventLogger
#' 
#' @examples 
#' # load model
#' path <- system.file('extdata', 'SI.yaml', package='twt')
#' 
#' # load file and parse to construct MODEL object
#' settings <- yaml.load_file(path)
#' mod <- Model$new(settings)
#' 
#' # simulate outer tree - returns a Run object carrying EventLogger
#' run <- sim.outer.tree(mod)
#' 
#' # simulate inner tree
#' tree <- sim.inner.tree(run)
#' plot(tree)
#' 
#' @export
sim.inner.tree <- function(model, e=NA) {
  
  # handle arguments
  if ( is.element('Run', class(mod)) ) {
    run <- mod
    if ( nrow(run$get.eventlog()$get.all.events()) ==0 ) {
      stop("Error in sim.inner.tree(): empty EventLogger in Run object. ",
           "Did you run sim.outer.tree() before calling this function?")
    }
  }
  else if ( is.element('Model', class(mod)) ) {
    if (is.environment(e)) {
      run <- Run$new(mod)
      # don't modify the original EventLogger object
      run$set.eventlog(e$clone())
    } 
    else {
      stop("You must provide an EventLogger to sim.inner.tree() if `mod` is a Model object")
    }
  }
  else {
    stop("sim.inner.tree(): argument `mod` must be R6 class Model or Run.")
  }
  
  eventlog <- run$get.eventlog()
  events <- eventlog$get.all.events()
  events <- events[order(events$time), ]
  
  if (any(!is.element(run$get.compartments(), 
                      union(t.events$compartment1, t.events$compartment2)))) {
    stop(paste("Mismatch between compartment names in Run and EventLogger objects."))
  }
  
  # initialize simulation at most recent sampling time
  current.time <- 0.
  extant.lineages <- run$get.extant.lineages(current.time)
  num.extant <- length(extant.lineages)
  if (num.extant == 0) {
    stop ('There must be at least one lineage sampled at time t=0.')
  }
  
  for (i in 1:nrow(events)) {
    e <- events[i,]
    
    # sample waiting times to coalescent events
    c.times <- calc.coal.wait.times(run, current.time)
    
    if (length(c.times) == 0) {
      # no Compartments with two or more Lineages
      # resolve the next event
      
      if (e$event.type == 'transmission') {
        resolve.transmission(run, e)
      }
      else if (e$event.type == 'migration') {
        resolve.migration(run, e)
      }
    }
    else {
      # compete coalescent event against next outer event
      waiting.time <- e$time - current.time
      
    }
  }
  
}






