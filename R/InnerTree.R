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
    initialize = function(outer, mod, prefix='P', p.index=1) {
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
      
      private$mod <- mod
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
      private$inactive$add.host(index.case$clone())
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
  
  root <- index.case$get.pathogens()[[1]]
  
  tips <- events$pathogen1[events$event=='sampling']
}


#' @keywords
.relabel.inner.events <- function(events, tips) {
  # re-order events in forward time
  events <- events[order(events$time), ]
  
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
}