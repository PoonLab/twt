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
    },
    
    count.active = function() {
      private$active$count.type()
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
  children <- unique(events$to.host[events$from.host==node])
  inf.times <- sapply(children, function(child) 
    events$time[events$to.host==child])
  
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
plot.OuterTree <- function(obj, pad=1.05) {
  # retrieve sampling times
  sampled <- obj$get.sampled()
  stopifnot(sampled$count.type() > 0)
  samp.times <- setNames(
    as.numeric(sampled$get.sampling.times()), sampled$get.names())
  sampled <- sampled$get.names()
  
  # there should be only one Host in `active`
  active <- obj$get.active()
  stopifnot(active$count.type()==1)
  root <- active$sample.host()
  root.name <- root$get.name()
  
  hosts <- obj$get.retired()
  host.names <- hosts$get.names()
  
  # FIXME: what happens if there are multiple? i.e., superinfection
  #inf.times <- setNames(hosts$get.transmission.times(), host.names)
  sources <- setNames(
    lapply(hosts$get.sources(), function(h) h$get.name()), host.names)
  
  # starting from tips, build edge list
  events <- obj$get.log()
  events$time <- as.numeric(events$time)
  trans <- events[events$event=='transmission', ]
  
  # reduce to the earliest transmission to each host
  first.idx <- sapply(split(1:nrow(trans), trans$to.host), function(idx) {
    idx[which.min(trans$time[idx])]
  })
  first.trans <- trans[first.idx, ]
  
  # sort nodes by post-order traversal (children before parents)
  nodes <- .reorder.events(first.trans, root.name)
  
  # prepare plot region
  par(mar=c(5,1,1,5), mfrow=c(1,1))
  plot(NA, xlim=c(0, max(trans$time)*pad), ylim=c(0.5, length(nodes)+0.5),
       xlab="Time", yaxt='n', ylab=NA, bty='n')
  
  for (node in nodes) {
    is.sampled <- is.element(node, sampled)
    
    if (node == root.name) {
      inf.time <- 0
      source <- NA
    } else {
      inf.time <- first.trans$time[first.trans$to.host==node]
      source <- sources[[node]]
    }
    
    samp.time <- 0
    if (is.sampled) {
      samp.time <- samp.times[[node]]
    } else {
      # unsampled host - set right limit at first transmission event
      samp.time <- max(trans$time[trans$from.host==node])
    }
    
    segments(x0=inf.time, x1=samp.time, y0=which(nodes==node),
             lwd=2, lend=2, col=ifelse(is.sampled, 'black', 'grey'))
    text(x=max(samp.times)*pad, y=which(nodes==node), label=node, 
         xpd=NA, adj=0, cex=0.5)
    if (!is.na(source)) {
      arrows(x0=inf.time, y1=which(nodes==node), y0=which(nodes==source),
             length=0.08, lwd=2, col='orangered')      
    }
  }
}
