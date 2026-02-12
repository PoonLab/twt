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
#' @param events:  data.frame, transmission and migration events
#' @param node:  character, label a Host object named in at least one event
#' @param order:  character, 'preorder' or 'postorder' traversal
#' @param decreasing:  bool, sort events with respect to time (default: TRUE)
#' @param result:  character, node names are appended to this vector
#' @keywords internal
.reorder.events <- function(events, node, order='postorder', decreasing=TRUE,
                            result=c()) {
  if (order=='preorder') {
    result <- c(result, node)  # parent before children
  }
  children <- na.omit(unique(events$to.host[events$from.host==node]))
  if (length(children) == 0) { return(result) }
  
  inf.times <- sapply(children, function(child) 
    events$time[events$to.host==child])
  for (child in children[order(inf.times, decreasing=decreasing)]) {
    result <- .reorder.events(events, child, order, decreasing, result)
  }
  if (order=='postorder') {
    result <- c(result, node)  # children before parent
  }
  return(result)
}


#' plot.OuterTree
#' S3 generic plot method for objects of class `OuterTree`.  This visualizes
#' the collection of events sampled in the outer simulation as a transmission 
#' tree.
#' @param obj:  R6 object of class `OuterTree`
#' @param pad:  numeric, controls padding around tree in plot region
#' @export
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
  targets <- obj$get.targets()
  
  # get root node
  active <- obj$get.active()
  stopifnot(active$count.type()==1)
  root <- active$sample.host()
  root.name <- root$get.name()
  
  events <- obj$get.log()
  events$time <- as.numeric(events$time)
  
  events <- .filter.firsts(events)  # remove superinfections
  events <- .relabel.nodes(events, targets)
  
  # re-order events by preorder traversal (parents before children)
  preorder <- .reorder.events(events, root.name, order='preorder',
                              decreasing=FALSE)
  idx <- match(preorder, events$to.host)
  events <- events[idx[-1], ]
  
  root.time <- min(events$time[events$from.host==root.name])
  #events$branch.len <- ifelse(events$from.host==)
  ## FIXME: WORK IN PROGRESS - unfold transmission history into tree
  
  # each event creates a node - need to convert node list to edge list
  n.edges <- nrow(events)-1
  edge.list <- data.frame(parent=character(n.edges), 
                          child=character(n.edges),
                          length=numeric(n.edges),
                          label=character(n.edges))
  for (i in 2:nrow(events)) {
    e <- events[i,]
    
    if (events$to.host[i-1] == e$from.host) {
      e.prev <- events[i-1, ]  # next step in chain of events
      
    } else if (any(events$to.host[1:(i-1)] == e$from.host & 
                   events$time[1:(i-1)] < e$time)) {
      # preceding event from same host
      temp <- events[1:(i-1), ]
      temp <- temp[temp$to.host == e$from.host, ]
      e.prev <- temp[which.max(temp$time), ]
    } else if (e$from.host == root.name) {
      e.prev <- events[1,]  # root is a special case
    }
    edge.list[i-1, ] <- c(
      parent=paste(e.prev$from.host, e.prev$to.host, sep="__"),
      child=paste(e$from.host, e$to.host, sep="__"),
      length=e$time - e.prev$time,
      label=e$to.host
    )
  }
  
  
  node.label <- preorder[!is.element(preorder, tips)]
  tip.label <- preorder[is.element(preorder, tips)]
  nodes <- c(tip.label, node.label)
  
  edge <- matrix(0, nrow=nrow(edge.list), ncol=2)
  edge[,1] <- match(edge.list$parent, nodes)
  edge[,2] <- match(edge.list$child, nodes)
  
  phy <- list(
    tip.label = tips,
    Nnode = length(internals),
    edge = edges,
    edge.length = edge.length
  )
  attr(phy, 'class') <- 'phylo'
  phy  # return object
}


#' Remove any superinfections, keeping only the first transmission event 
#' to each host.
#' @param events:  data frame
#' @keywords internal
.filter.firsts <- function(events) {
  # extract the earliest transmission to each host
  idx <- which(events$event=='transmission')
  first.idx <- sapply(split(idx, events$to.host[idx]), function(i) {
    i[which.min(events$time[i])]
  })
  remove <- idx[!is.element(idx, first.idx)]
  if (length(remove) > 0) {
    events[-remove, ]
  } else {
    events
  }
}


#' Modify host names in the event log to make them unique.  This assumes 
#' that we have removed superinfection events from the log, i.e., with 
#' `.filter.firsts`, so there is a linear sequence of events for each host.
#' We assume the first event is the transmission of infection to the host.
#' Only migration events (not sampling) require a new label.
#' @keywords internal
.relabel.nodes <- function(events, targets) {
  # make sure events are ordered in time
  events <- events[order(events$time), ]
  
  # relabel host names for migration events
  nodes <- na.omit( unique(c(events$from.host, events$to.host)) )
  for (node in nodes) {
    idx <- which(events$from.host == node | events$to.host==node)
    new.node <- node
    label <- 1
    for (i in idx) {
      if (events$from.host[i] == node) {
        # overwrite node name in case it's been updated
        events$from.host[i] <- new.node        
      }
      
      if (events$event[i] == 'migration') {
        if (events$to.comp[i] %in% names(targets)) {
          # terminal node, can happen only once
          new.node <- paste(node, "sample", sep="_")
        } else {
          new.node <- paste(node, label, sep="_")  # internal node
          label <- label + 1
        }
        events$to.host[i] <- new.node
      }
    }
  }
  events
}
