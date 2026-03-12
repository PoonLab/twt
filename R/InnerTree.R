#' InnerTree
#' 
#' \code{InnerTree} is an R6 class defining a persistent object that records 
#' events while simulating the inner tree within an outer tree.
#' 
#' The following types of events are recorded in the `inner.log`:
#'   sampling: `pathogen1` in `from.host`
#'   transmission: `from.host`->`to.host` of `pathogen1`
#'   coalescence: `pathogen1` and `pathogen2` within `from.host`
#' 
#' @param mod:  R6 object of class Model
#' @param prefix:  character, used for labeling Pathogen objects (default 'P')
#' @param p.index:  integer, counter for uniquely labeling Pathogens
#' @export
InnerTree <- R6Class(
  "InnerTree",
  public = list(
    initialize = function(outer, prefix='P', p.index=1) {
      private$inner.log <- data.frame(
        time=numeric(),  # time of event
        event=character(),  # type of event, e.g., coalescence
        from.comp=character(),
        to.comp=character(),
        from.host=character(),  # source Host name (transmission only)
        to.host=character(),  # recipient Host name (transmission only)
        pathogen1=character(),
        pathogen2=character()
      )
      
      private$mod <- outer$get.model()  # inherits model from OuterTree
      private$prefix <- prefix
      private$p.index <- p.index
      
      # extract and clone variables from OuterTree
      root <- outer$get.active()
      if (root$count.type() != 1) {
        stop("OuterTree does not converge to a single Host")
      }
      index.case <- root$sample.host()
      
      # inactive Hosts are part of the transmission history
      private$inactive <- HostSet$new()
      private$sampled <- HostSet$new()
      
      # note there is overlap between `retired` and `sampled` HostSets
      hosts <- outer$get.retired()
      sampled <- outer$get.sampled()
      for (host in hosts$get.hosts()) {
        if (host$get.name() %in% sampled$get.names()) { 
          private$sampled$add.host(host$clone()) 
        } else {
          private$inactive$add.host(host$clone())
        }
      }
      
      # handle index case
      if (index.case$get.name() %in% sampled$get.names()) {
        private$sampled$add.host(index.case$clone())
      } else {
        private$inactive$add.host(index.case$clone())
      }
      
      # track Hosts with active Pathogen lineages
      private$active <- HostSet$new()
      
    },
    
    # accessor functions
    get.log = function() { private$inner.log },
    get.model = function() { private$mod },
    get.prefix = function() { private$prefix },
    get.pindex = function() { private$p.index },
    
    get.inactive = function() { private$inactive },
    get.active = function() { private$active },
    n.active = function() { length(private$active) },
    get.sampled = function() { private$sampled },
    
    has.target = function(cname) { 
      is.element(cname, names(private$mod$get.sampling()$targets))
    },
    
    new.pathogen = function(time) {
      path <- Pathogen$new(
        name=paste(private$prefix, private$p.index, sep="_"),
        end.time=time,
      )
      private$p.index <- private$p.index + 1  # increment counter
      return(path)
    },
    
    add.event = function(e) {
      if ( !all(names(e) == names(private$inner.log)) ) {
        stop("Event passed to InnerTree:add.event must be vector with ",
             "names: ", names(private$inner.log))
      }
      private$inner.log[nrow(private$inner.log)+1, ] <- e
    }
  ),
  
  private = list(
    inner.log = NULL,
    mod = NULL,
    prefix = NULL,
    p.index = NULL,
    inactive = NULL,
    active = NULL,
    sampled = NULL
  )
)


#' as.phylo.InnerTree
#' 
#' A generic S3 method for converting an InnerTree object to an S3 object of 
#' class `ape::phylo` that represents a phylogenetic tree.  This inner 
#' phylogeny is nested within an "outer" transmission tree.  
#' 
#' @param obj:  R6 object of class `InnerTree`
#' @return object of class `phylo`
#' @export
as.phylo.InnerTree <- function(obj) {
  events <- obj$get.log()
  events$time <- as.numeric(events$time)
  
  active <- inner$get.active()
  if (active$count.type() != 1) {
    stop("Error, expected only one Host in active HostSet")
  }
  
  index.case <- active$get.hosts()[[1]]
  if (index.case$count.pathogens() != 1)  {
    stop("Error: expected only one Pathogen in index.case")
  }
  
  root <- index.case$get.pathogens()[[1]]$get.name()
  events <-  .relabel.inner.events(events)
  
  preorder <- .reorder.inner.events(events, root, order='preorder')
  idx <- match(preorder, events$pathogen2)
  events <- events[idx[-1], ]  # reorder events
  
  root.time <- min(events$time[events$pathogen1==root])
  
  n.edges <- nrow(events)-1
  edge.list <- data.frame(parent=character(n.edges), 
                          child=character(n.edges),
                          length=numeric(n.edges),
                          label1=character(n.edges),
                          label2=character(n.edges),
                          compartment=character(n.edges))
  for (i in 2:nrow(events)) {
    e <- events[i,]
    
    if (events$pathogen2[i-1] == e$pathogen1) {
      e.prev <- events[i-1, ]  # next step in chain of events
    } else {
      # preceding event involving same pathogen
      temp <- events[1:(i-1), ]
      temp <- temp[temp$pathogen1 == e$pathogen1, ]
      e.prev <- temp[which.max(temp$time), ]
    }
    
    edge.list[i-1, ] <- c(
      parent=paste(e.prev$pathogen1, e.prev$pathogen2, sep="__"),
      child=paste(e$pathogen1, e$pathogen2, sep="__"),
      length=e$time - e.prev$time,
      label1=e.prev$pathogen2,
      label2=e$pathogen2,
      compartment=e$from.comp
    )
  }
  
  preorder <- preorder[-1]  # discard first label
  node.label <- preorder[!grepl("_sample$", preorder)]
  tip.label <- preorder[grepl("_sample$", preorder)]
  nodes <- c(tip.label, node.label)
  
  edge  <- matrix(0, nrow=nrow(edge.list), ncol=2)
  edge[,1] <- match(edge.list$label1, nodes)
  edge[,2] <- match(edge.list$label2, nodes)
  
  phy <- list(
    tip.label=tip.label, node.label=node.label,
    Nnode=length(node.label), edge=edge, 
    edge.length=as.numeric(edge.list$length),
    event=events$event[match(nodes, events$pathogen2)],
    compartment=edge.list$compartment
  )
  attr(phy, 'class') <- 'phylo'
  phy
}


#' Modify Pathogen names in the inner tree event log to make them unique.
#' This makes it easier to label single nodes (transmission, migration events)
#' in the tree. 
#' 
#' @param events:  data frame including fields `time`, `event`, `pathogen1` and 
#'                 `pathogen2`
#' @return data frame
#' 
#' @keywords internal
#' @noRd
.relabel.inner.events <- function(events) {
  # re-order events in forward time
  events$event <- factor(
    events$event, 
    levels=c("transmission", "coalescent", "migration", "sampling"))
  events <- events[order(events$time, events$event), ]
  
  nodes <- na.omit(unique(c(events$pathogen1, events$pathogen2)))
  for (node in nodes) {
    idx <- which(events$pathogen1==node | events$pathogen2==node)
    new.node <- node  # for storing the revised label
    label <- 1
    for (i in idx) {
      if (events$pathogen1[i] == node) {
        # overwrite original node name in case it's been updated
        events$pathogen1[i] <- new.node
      }
      
      if (is.na(events$pathogen2[i])) {
        if (events$event[i] == 'sampling') {
          new.node <- paste(node, "sample", sep="_")
        } else {
          new.node <- paste(node, label, sep="_")
          label <- label + 1
        }
        events$pathogen2[i] <- new.node
      }
    }
  }
  return(events)
}


#' .reorder.inner.events
#' Recursive function for traversing inner tree event log.
#' 
#' @param events:  data.frame
#' @param node:  character, Pathogen name
#' @param order:  character, 'preorder' or 'postorder'
#' @param result:  character, vector to append results by recursive calls
#' @return character, Pathogen names ordered by tree traversal
#' 
#' @keywords internal
#' @noRd
.reorder.inner.events <- function(events, node, order='postorder', result=c()) {
  if (order == 'preorder') {
    result <- c(result, node)
  }
  children <- na.omit(unique(events$pathogen2[events$pathogen1==node]))
  for (child in children) {
    result <- .reorder.inner.events(events, child, order, result)
  }
  if (order=='postorder') {
    result <- c(result, node)
  }
  return(result)
}


#' Generic print function for R6 objects of class `InnerTree`
#' @export
#' @noRd
print.InnerTree <- function(obj) {
  cat("\033[93m\033[1mtwt InnerTree\033[22m\033[37m\n")  # bold color!
  cat(" ", obj$get.sampled()$count.type(), "sampled Pathogens\n")
  cat(" ", obj$get.active()$count.type(), "active Pathogens\n")
  cat(" ", obj$get.inactive()$count.type(), "inactive Pathogens\n")
  events <- obj$get.log()
  cat(" ", nrow(events), "events in inner log: ")
  print(table(events$event))
} 
