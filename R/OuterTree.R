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
.reorder.events <- function(events, node, order='postorder', result=c()) {
  if (order=='preorder') {
    result <- c(result, node)  # parent before children
  }
  
  children <- unique(events$to.host[events$from.host==node])
  
  inf.times <- sapply(children, function(child) 
    events$time[events$to.host==child])
  
  for (child in children[order(inf.times, decreasing=TRUE)]) {
    result <- .reorder.events(events, child, order, result)
  }
  
  if (order=='postorder') {
    result <- c(result, node)  # children before parent
  }
  
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

#' Generic S3 method for converting an OuterTree to a phylo object
#' 
#' `ape::phylo` is an S3 class object that represents a phylogenetic tree. We 
#' use this object class to represent the transmission tree.  Because this 
#' object is, by definition, an acyclic graph, it cannot represent cases of 
#' superinfection.  For this reason, only the earliest transmission event is 
#' used for each recipient host.  To represent the entire transmission 
#' history including superinfection, use `as.igraph.OuterTree`.
#' 
#' Note that interpreting the resulting transmission tree as a pathogen
#' phylogeny ignores all within-host evolution, e.g., lineage sorting.
#' 
#' @param obj:  R6 object of class `OuterTree`
#' @param singles:  bool, keep single nodes that represent transmission or 
#'                  migration events with only one descendant branch.
#' @return S3 object of class `phylo`
#' @export
as.phylo.OuterTree <- function(obj, singles=TRUE) {
  sampled <- obj$get.sampled()  # HostSet
  tips <- sampled$get.names()  # most recent sample first
  samp.times <- setNames(as.numeric(sampled$get.sampling.times()), tips)
  
  # get root node
  active <- obj$get.active()
  stopifnot(active$count.type()==1)
  root <- active$sample.host()
  root.name <- root$get.name()
  
  events <- obj$get.log()
  events$time <- as.numeric(events$time)
  
  # extract the earliest transmission to each host
  trans <- events[events$event=='transmission', ]
  first.idx <- sapply(split(1:nrow(trans), trans$to.host), function(idx) {
    idx[which.min(trans$time[idx])]
  })
  first.trans <- trans[first.idx, ]
  
  # sort nodes by preorder traversal
  preorder <- .reorder.events(first.trans, root.name, order='preorder')
  
  # reorder node list
  tips <- preorder[is.element(preorder, tips)]
  internals <- preorder[!is.element(preorder, tips)]
  nodes <- c(tips, internals)
  
  # generate edge list
  order.trans <- first.trans[na.omit(match(preorder, first.trans$to.host)), ]
  edges <- matrix(c(
    sapply(order.trans$from.host, function(n) which(nodes==n)),
    sapply(order.trans$to.host, function(n) which(nodes==n))
  ), byrow=FALSE, ncol=2)
  
  # locate first split
  root.time <- min(order.trans$time[
    order.trans$event=='transmission' & order.trans$from.host==root.name])
  
  edge.length <- sapply(1:nrow(order.trans), function(i) {
    e <- order.trans[i,]
    parent <- e$from.host
    child <- e$to.host
    
    end.time <- ifelse(
      child %in% tips,  # branch terminates
      samp.times[[child]],  
      min(order.trans$time[order.trans$from.host==child])  # next event
      )
    
    start.time <- e$time
    if( sum(order.trans$from.host==parent) > 1 ) {
      # go to previous event
      idx <- which(order.trans$from.host==parent)
      precedes <- (order.trans$time[idx] < e$time)
      if (any(precedes)) {
        start.time <- order.trans$time[idx][precedes]
      }
    }
    end.time - start.time  # return value
  })
  
  phy <- list(
    tip.label = tips,
    Nnode = length(internals),
    edge = edges,
    edge.length = edge.length
  )
  attr(phy, 'class') <- 'phylo'
  phy  # return object
}


