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
        time=numeric(),  # time of event
        event=character(),  # type of event, e.g., migration
        from.comp=character(),  # compartment before event
        to.comp=character(),  # compartment after event
        from.host=character(),  # source Host name (transmission only)
        to.host=character()  # recipient Host name (transmission only)
      )
      
      private$targets <- mod$get.sampling()$targets
      private$sampled <- HostSet$new()  # store sampled Hosts
      private$active <- HostSet$new()  # store Hosts carrying Pathogens
      private$retired <- HostSet$new()  # previously active Hosts
      private$is.infected <- mod$get.infected()  # copy whole named vector
    },
    
    # accessor functions
    get.targets = function() { private$targets }, 
    get.log = function() { private$outer.log },
    get.nrow = function() { nrow(private$outer.log) },
    get.sampled = function() { private$sampled },
    get.retired = function() { private$retired },
    get.active = function() { private$active },
    get.infected = function(cn=NA) { 
      if (is.na(cn)) { private$is.infected } else {
        if( is.element(cn, names(private$is.infected)) ) {
          private$is.infected[[cn]]
        } else { NA }
      }
    },
    nsamples = function() { private$sampled$count.type() },
    
    add.event = function(e) {
      if ( !all(names(e) == names(private$outer.log)) ) {
        stop("Event passed to OuterTree:add.event must be vector with ",
             "names: ", names(private$outer.log))
      }
      private$outer.log[nrow(private$outer.log)+1, ] <- e
    },
    
    add.sample = function(host) {
      private$sampled$add.host(host)
    },
    
    add.retired = function(host) {
      private$retired$add.host(host)
    }
  ),
  
  private = list(
    outer.log = NULL,
    targets = NULL,
    sampled = NULL,
    active = NULL,
    retired = NULL,
    is.infected = NULL
  )
)


#' .reorder.events
#' 
#' A helper function that recursively generates a list of node labels
#' by post-order traversal (outputting parents after children).
#' 
#' @keywords internal
.reorder.events <- function(events, node, result=c()) {
  children <- events$from.comp[events$to.comp==node]
  inf.times <- sapply(children, function(child) 
    events$time[events$from.comp==child])
  
  for (child in children[order(inf.times, decreasing=TRUE)]) {
    result <- .reorder.events(events, child, result)
  }
  result <- c(result, node)
  return(result)
}


#' plot.OuterTree
#' S3 generic plot method for objects of class `outerTree`.  This visualizes
#' the collection of events sampled in the outer simulation as a transmission 
#' tree.
plot.OuterTree <- function(obj) {
  # retrieve sampling times
  sampled <- obj$get.sampled()
  samp.times <- sampled$get.sampling.times()
  
  retired <- obj$get.retired()
  active <- obj$get.active()  # should only be one
  host.names <- c(active$get.names(), retired$get.names())
  host.types <- c(active$get.types(), retired$get.types())
  sources <- c(active$get.)
  
  events <- obj$get.log()
  trans <- events[events$event=='transmission', ]
  
  # find root
  sources <- unique(trans$from.host)
  root <- sources[!is.element(sources, trans$to.host)]
  if (length(root) > 1) {
    stop("Detected multiple roots in outer tree")
  }
  
  # sort nodes by post-order traversal
  nodes <- .reorder.events(trans, root)
  
  # prepare plot region
  par(mar=c(5,5,1,1), mfrow=c(1,1))
  plot(NA, xlim=c(-max(trans$time*1.05), 0), ylim=c(0.5, length(nodes)+0.5),
       xlab="Time", yaxt='n', ylab=NA, bty='n')
  
  for (node in nodes) {
    
  }
}