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
      private$model <- mod
      
      private$outer.log <- data.frame(
        time=numeric(),  # time of event
        event=character(),  # type of event, e.g., migration
        from.comp=character(),  # compartment before event
        src.comp=character(),  # compartment of source Host (transmission only)
        to.comp=character(),  # compartment after event
        from.host=character(),  # Host name (source in case of transmission)
        to.host=character()  # recipient Host name (transmission only)
      )
      
      private$targets <- mod$get.sampling()$targets
      private$sampled <- HostSet$new()  # store sampled Hosts
      private$active <- HostSet$new()  # store Hosts carrying Pathogens
      private$retired <- HostSet$new()  # previously active Hosts
      private$is.infected <- mod$get.infected()  # copy whole named vector
    },
    
    # accessor functions
    get.model = function() { private$model },
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
    model = NULL,
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
#' by post-order (outputting children before parent), or preorder (parents
#' before children) traversal.  When there are multiple children from a node,
#' they are ordered in decreasing (increasing) order by event time.
#' 
#' @param events:  data.frame, transmission and migration events; must have 
#'                 `time`, `from.host` and `to.host` fields.
#' @param node:  character, label a Host object named in at least one event
#' @param order:  character, 'preorder' or 'postorder' traversal
#' @param decreasing:  bool, sort events with respect to time (default: TRUE)
#' @param result:  character, node names are appended to this vector
#' 
#' @keywords internal
#' @noRd
.reorder.events <- function(events, node, order='postorder', decreasing=TRUE,
                            result=c()) {
  if (order=='preorder') {
    result <- c(result, node)  # append parent before children
  }
  children <- na.omit(unique(events$to.host[events$from.host==node]))
  inf.times <- sapply(children, function(child) 
    events$time[events$to.host==child])
  for (child in children[order(inf.times, decreasing=decreasing)]) {
    result <- .reorder.events(events, child, order, decreasing, result)
  }
  if (order=='postorder') {
    result <- c(result, node)  # append children before parent
  }
  return(result)
}


#' .find.node.at
#'
#' Given an original (pre-relabeling) host name and a time point, return the
#' y-position (row index in \code{nodes}) of that host's segment in the
#' plotted outer tree.  The name of a host's segment changes whenever a
#' migration event occurs; this function traces those name changes to find
#' which relabeled segment is active at \code{time}.
#'
#' @param orig.name  character, original host name before relabeling
#' @param time       numeric, the time point of interest
#' @param rel.events data.frame, the relabeled filtered event log (output of
#'                   \code{.relabel.nodes})
#' @param nodes      character vector, ordered node names from
#'                   \code{.reorder.events} (defines y-positions)
#' @return integer y-position, or \code{NA} if not found
#' @keywords internal
#' @noRd
.find.node.at <- function(orig.name, time, rel.events, nodes) {
  current.name <- orig.name
  # trace through migration events in time order to follow name changes
  mig <- rel.events[rel.events$event == 'migration', , drop=FALSE]
  mig <- mig[order(mig$time), , drop=FALSE]
  for (i in seq_len(nrow(mig))) {
    if (mig$from.host[i] == current.name && mig$time[i] <= time) {
      current.name <- mig$to.host[i]
    }
  }
  y <- which(nodes == current.name)
  if (length(y) == 0L) NA_integer_ else y
}


#' plot.OuterTree
#' S3 generic plot method for objects of class `OuterTree`.  This visualizes
#' the collection of events sampled in the outer simulation as a transmission
#' tree.  Primary transmissions are drawn as solid orangered arrows;
#' superinfection events (secondary transmissions to an already-infected host)
#' are drawn as dashed steelblue arrows.
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
  
  # retrieve full event log before any filtering
  events <- obj$get.log()
  events <- events[!(events$event=='migration' & is.na(events$to.host)), ]
  events$time <- as.numeric(events$time)
  
  # identify superinfection events: any transmission to a host that already
  # received a (first) transmission, sorted by time so duplicated() works
  trans <- events[events$event == 'transmission', , drop=FALSE]
  trans <- trans[order(trans$time), , drop=FALSE]
  si.events <- trans[duplicated(trans$to.host), , drop=FALSE]
  
  # filter and relabel for primary-tree layout (unchanged logic)
  events <- .filter.firsts(events)
  events <- .relabel.nodes(events, obj$get.targets())
  
  # sort nodes by preorder traversal (parents before children)
  nodes <- .reorder.events(events, root.name, order="preorder", decreasing=FALSE)
  
  # prepare plot region
  par(mar=c(5,1,1,5), mfrow=c(1,1))
  plot(NA, xlim=c(0, max(events$time)*pad), ylim=c(0.5, length(nodes)+0.5),
       xlab="Time", yaxt='n', ylab=NA, bty='n')
  
  # primary transmission tree
  for (node in nodes) {
    is.sampled <- is.element(node, sampled)
    
    if (node == root.name) {
      inf.time <- 0
      source <- NA
    } else {
      inf.time <- events$time[events$to.host==node]
      source <- events$from.host[events$to.host==node]
    }
    
    samp.time <- 0
    if (is.sampled) {
      samp.time <- samp.times[[node]]
    } else {
      # unsampled host - right limit is its first outgoing transmission
      out.times <- events$time[events$from.host==node]
      samp.time <- if (length(out.times) > 0) max(out.times) else inf.time
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
  
  # superinfection arrows (dashed steelblue)
  if (nrow(si.events) > 0L) {
    for (i in seq_len(nrow(si.events))) {
      si <- si.events[i, ]
      donor.y <- .find.node.at(si$from.host, si$time, events, nodes)
      recip.y <- .find.node.at(si$to.host,   si$time, events, nodes)
      if (!is.na(donor.y) && !is.na(recip.y)) {
        arrows(x0=si$time, y0=donor.y, y1=recip.y,
               length=0.08, lwd=2, lty=2, col='steelblue')
      }
    }
    # legend only when superinfection events are present
    legend('bottomright', legend=c('transmission', 'superinfection'),
           col=c('orangered', 'steelblue'), lty=c(1, 2), lwd=2, bty='n', cex=0.8)
  }
}

#' as.phylo.OuterTree
#' 
#' A generic S3 method for converting an OuterTree to a phylo object.
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
#' @return S3 object of class `phylo`
#' @export
as.phylo.OuterTree <- function(obj) {
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
  
  # each event creates a node - need to convert node list to edge list
  n.edges <- nrow(events)-1
  edge.list <- data.frame(parent=character(n.edges), 
                          child=character(n.edges),
                          length=numeric(n.edges),
                          label1=character(n.edges),
                          label2=character(n.edges),
                          compartment=character(n.edges))
  for (i in 2:nrow(events)) {
    e <- events[i,]
    
    if (events$to.host[i-1] == e$from.host) {
      e.prev <- events[i-1, ]  # next step in chain of events
    } else {
      # preceding event from same host
      temp <- events[1:(i-1), ]
      temp <- temp[temp$from.host == e$from.host, ]
      e.prev <- temp[which.max(temp$time), ]
    }
    edge.list[i-1, ] <- c(
      parent=paste(e.prev$from.host, e.prev$to.host, sep="__"),
      child=paste(e$from.host, e$to.host, sep="__"),
      length=e$time - e.prev$time,
      label1=e.prev$to.host,
      label2=e$to.host,
      compartment=ifelse(e$event=='migration', e$from.comp, e$src.comp)
    )
  }
  
  preorder <- preorder[-1]  # discard first host label
  node.label <- preorder[!grepl("_sample$", preorder)]
  tip.label <- preorder[grepl("_sample$", preorder)]
  nodes <- c(tip.label, node.label)
  
  edge <- matrix(0, nrow=nrow(edge.list), ncol=2)
  edge[,1] <- match(edge.list$label1, nodes)
  edge[,2] <- match(edge.list$label2, nodes)
  
  phy <- list(
    tip.label = tip.label,
    node.label = node.label,
    Nnode = length(node.label),
    edge = edge,
    edge.length = as.numeric(edge.list$length),
    event = events$event[match(nodes, events$to.host)],
    compartment = edge.list$compartment,
    time = events$time[-1] - root.time
  )
  attr(phy, 'class') <- 'phylo'
  phy  # return object
}


#' .filter.firsts
#' Remove any superinfections, keeping only the first transmission event 
#' to each host.
#' @param events:  data frame with fields `time`, `event` and `to.host`
#' @keywords internal
#' @noRd
.filter.firsts <- function(events) {
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


#' .relabel.nodes
#' Modify host names in the event log to make them unique.  This assumes 
#' that we have removed superinfection events from the log, i.e., with 
#' `.filter.firsts`, so there is a linear sequence of events for each host.
#' We assume the first event is the transmission of infection to the host.
#' Only migration events require a new label.
#' 
#' @param events:  data frame with fields `time`, `event`, `from.comp`, 
#'                 `to.comp`, `from.host`, and `to.host`
#' @param targets:  list from Model$get.sampling()$targets
#' @keywords internal
#' @noRd
.relabel.nodes <- function(events, targets) {
  # make sure events have been filtered of superinfection
  n.inf <- table(events$to.host)
  if (any(n.inf > 1)) {
    stop("ERROR: .relabel.nodes assumes events have been filtered of ",
         "superinfection events.")
  }
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


#' Generic print function for R6 objects of class `OuterTree`
#' @export
#' @noRd
print.OuterTree <- function(obj) {
  cat("twt OuterTree\n")  # bold color!
  cat(" ", obj$get.sampled()$count.type(), "sampled Hosts\n")
  cat(" ", obj$get.active()$count.type(), "active Hosts\n")
  cat(" ", obj$get.retired()$count.type(), "retired Hosts\n")
  events <- obj$get.log()
  cat(" ", nrow(events), "events in outer log:")
  print(table(events$event))
}
