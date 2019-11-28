#' EventLogger
#' 
#' \code{EventLogger} is an R6 class for an object that tracks migration, 
#' transmission, and coalescent events.  Note that bottleneck events are logged as 
#' coalescent events.
#' 
#' @param events: a data frame where each row represents an event with a time 
#' stamp in forward time.
#' 
#' @examples
#' # manually initialize an EventLog object
#' e <- EventLogger$new()
#' # note this log entry is not linked to existing Lineage or Compartment objects
#' e$add.event("transmission", time=1, line1="NA", comp1="hist1", comp2="host2")
#' e$get.all.events()
#' 
#' @export
EventLogger <- R6Class("EventLogger", 
  public = list(
    initialize = function(events = NA) {
      if (is.na(events)) {
        self$clear.events()
      } else {
        private$events <- events  
      }
    },
   
   
    get.all.events = function() {
      if (nrow(private$events) == 0) {
        #cat('No events to display.')
        NULL
      } else {
        private$events  # default eventlog shows cumulative time b/c more user friendly
      }
    },
    
    
    get.events = function(event.type) {
      eventList <- private$events[ which(private$events$event.type == event.type), ]
      
      if (nrow(eventList) != 0) {
        eventList
      } else {
        # cat('No events of type "', event.type, '".\n')
        NULL
      }
    },
    
    
    #' @param type: event type, one of 'transmission', 'transition', 'migration', 
    #'        'coalescence' or 'bottleneck'.
    #' @param time: CUMULATIVE time that event has occurred between two compartments 
    #'        in a transmission/migration/coalescent event
    add.event = function(type, time, line1=NA, line2=NA, comp1=NA, comp2=NA, 
                         type1=NA, type2=NA) {
      
      if ( !is.element(type, c(
        'transmission', 'migration', 'transition', 'coalescent', 'bottleneck'
        ))) {
        stop("Error in EventLogger:add.event(), unrecognized event type ", type)
      }
      
      e <- list(event.type=type, time=time, lineage1=line1, lineage2=line2, 
                compartment1=comp1, compartment2=comp2, 
                type1=type1, type2=type2)
      private$events <- rbind(private$events, e, stringsAsFactors=F)
    }, 
    
    clear.events = function() {
      private$events <- self$blank.events()
    },
    
    blank.events = function() {
      data.frame(
        event.type=character(),
        time=numeric(),
        lineage1=character(),
        lineage2=character(),
        compartment1=character(),
        compartment2=character(),
        type1=character(),
        type2=character(),
        stringsAsFactors = FALSE
        )
    },
    
    
    #' Record which Lineages are transmitted from source to 
    #' recipient (should only be one entry with Compartment as recipient)
    #' 
    #' @param recipient:  Compartment object
    #' @param lineages:  a vector of names of Lineages to transfer out
    #'                   of recipient Compartment
    record.transmission = function(recipient, lineages) {
      # locate transmission event
      events <- self$get.all.events()
      idx <- which( events$compartment1 == recipient$get.name() &
                      events$event.type == 'transmission' )
                      #is.na(as.logical(events$lineage1)) )
      
      if (length(idx) != 1) {
        stop("Error in eventLogger:record.transmission - ",
             ifelse(length(idx)>1, "multiple", "no"), 
             " transmission events match recipient ",
             recipient$get.name())
      }
      
      # replace event with new rows
      e <- events[idx, ]
      cache <- events[(idx+1):nrow(events), ]
      
      if (idx > 1) {
        events <- events[1:(idx-1), ]  
      } else {
        # blank data frame
        events <- self$blank.events()
      }
      
      for (l in lineages) {
        e$lineage1 <- l$get.name()
        events <- rbind(events, as.list(e), stringsAsFactors=F)
      }
      private$events <- rbind(events, cache)
    },
    
    
    #' Record which Lineages are transmitted from source to recipient
    #' through a migration event
    #'
    #' @param recipient:  Compartment object
    #' @param source:  Compartment object
    #' @param time:  double, time of migration event
    #' @param lineages:  a list of Lineage objects
    record.migration = function(recipient, source, time, lineages) {
      events <- self$get.all.events()
      idx <- which(events$compartment1 == recipient$get.name() & 
                     events$compartment2 == source$get.name() & 
                     events$time == time &
                     is.na(events$lineage1))
      if (length(idx) != 1) {
        stop("Error in eventLogger:record.migration - ",
             ifelse(length(idx)>1, "multiple", "no"), 
             " migration events match recipient ",
             recipient$get.name(), " and source ", source$get.name())
      }
      
      e <- events[idx,]
      cache <- events[(idx+1):nrow(events), ]
      events <- events[1:(idx-1), ]
      for (l in lineages) {
        e$lineage1 <- l$get.name()
        events <- rbind(events, as.list(e), stringsAsFactors=F)
      }
      private$events <- rbind(events, cache)
    },
    
    
    get.fixed.samplings = function() {
      # retrieves the fixed sampling times of the tips of a MODEL object
      private$fixed.samplings.storage
    },
   
    store.fixed.samplings = function(model.fixed.samplings) {
      # stores the fixed sampling times of the tips of a MODEL object
      private$fixed.samplings.storage <- model.fixed.samplings
    }
    
  ),  # end public
  
  
  private = list(
    events = NULL,
    fixed.samplings.storage = NULL,
    
    generate.events = function(events, root, tips) {
      
      # inner recursive helper function
      generate.indiv.event <- function(node, parent_time) {
        # recursive function to generate cumulative times for each individual event
        # returns a childEvent or NULL to be added to the eventlog data frame
        if (node %in% tips) {
          return (NULL)
        } else {
          nodeEvents <- events[ which(events$compartment2 == node), ]
          if (length(row.names(nodeEvents)) == 0) {
            return(NULL)
          } else {
            for (x in 1:nrow(nodeEvents)) {
              childEvent <- nodeEvents[x,]
              # traverse descendants
              generate.indiv.event(
                as.character(childEvent['compartment1']), as.numeric(childEvent['time'])
              )
              childEvent['time'] <- parent_time - as.numeric(childEvent['time'])
              
              # append to data frame
              private$events <- rbind(private$events, childEvent, stringsAsFactors=F)
            }
            return(private$events)
          }
         
        }
      }
      
      # beginning of function generate.events()
      rootEvents <- events[ which(events$compartment2 == root), ]
      maxRootTime <- max(rootEvents$time)
      for (x in 1:nrow(rootEvents)) {
        parentEvent <- rootEvents[x,]
        # traverse descendants
        generate.indiv.event(as.character(parentEvent['compartment1']), as.numeric(parentEvent['time']))
        
        # root's individualt delta t from when it was infected to when it made its 
        # first transmission is 'undefined'
        # 0 or 1 by convention (see treeswithintrees closed issue #29)
        parentEvent['time'] <- maxRootTime - parentEvent['time']
        private$events.noncumul <- rbind(private$events, parentEvent, stringsAsFactors=F)
      }
      
      indices <- grep('NA', row.names(private$events), ignore.case=T, invert=T)
      match.cumul.ordering <- order(as.numeric(row.names(private$events[indices,])))
      
      private$events[indices,][match.cumul.ordering,]
    }
    
  )  # end private
)



#' print.EventLogger
#' 
#' `print.EventLogger` is simply a wrapper function on EventLogger's 
#' `get.all.events()`.
#' 
#' @param eventlog: object of class 'EventLog' 
#' 
#' @examples
#' e <- EventLogger$new()
#' e  # no events to display
#' 
#' e$add.event("transmission", time=1, line1="NA", comp1="hist1", comp2="host2")
#' e
#' 
#' @export
print.EventLogger <- function(eventlog) {
  events <- eventlog$get.all.events()
  if (is.null(events)) {
    print("No events to display.")
  } else {
    print(events)
  }
}


#' plot.EventLogger
#' Generic plot method for EventLogger class as an S3 object.
#' Converts the EventLogger object into an ape:phylo S3 object.
#' 
#' @param eventlog:  R6 object of class EventLogger
#' @param transmissions:  if TRUE, display transmission events on inner tree
#'        plot.
#' @param migrations: if TRUE, display migration events on inner tree plot.
#' @param node.labels: if TRUE, label internal nodes with names.
#' 
#' @examples 
#' path <- system.file('extdata', 'structSI.yaml', package='twt')
#' settings <- yaml.load_file(path)
#' structSI <- Model$new(settings)
#' run <- sim.outer.tree(structSI)
#' eventlog <- sim.inner.tree(run)
#' tr <- as.phylo(eventlog)
#' plot(tr)
#' 
#' @export
plot.EventLogger <- function(eventlog, transmissions, migrations, node.labels) {
  phy <- as.phylo.EventLogger(
    eventlog=eventlog, 
    transmissions=transmissions, 
    migrations=migrations, 
    node.labels=node.labels)
  
  # call plot.phylo S3 method
  plot(phy)
}



#' as.phylo
#' 
#' An S3 method converts events stored in the EventLogger object into an inner 
#' coalescent tree w/ option to include/exclude transmission events.
#' 
#' @param eventlog: EventLogger object
#' @param transmissions: logical; if TRUE, transmission events included, 
#' else excluded otherwise.
#' @return ape::phylo object
#' 
#' @examples 
#' path <- system.file('extdata', 'structSI.yaml', package='twt')
#' settings <- yaml.load_file(path)
#' structSI <- Model$new(settings)
#' run <- sim.outer.tree(structSI)
#' eventlog <- sim.inner.tree(run)
#' tr <- as.phylo(eventlog) 
#' tr
#' 
#' @export
as.phylo.EventLogger <- function(eventlog, transmissions=FALSE, migrations=FALSE, 
                                 node.labels=FALSE) {
  # retrieve fixed sampling times of tips
  fixed.sampl <- eventlog$get.fixed.samplings()
  
  # retrieve events from log
  events <- eventlog$get.all.events()
  core.events <- events[is.element(events$event.type, c('bottleneck', 'coalescent')), ]
  
  # separate nodes into root, tips, and internals
  # note lineage1 is child, lineage2 is parent
  root <- unique(core.events$lineage2[
    which(!is.element(core.events$lineage2, core.events$lineage1))
    ])
  
  if (length(root) > 1) stop('ERROR! Found multiple root nodes: ', root)
  if (length(root) == 0) stop("Error in as.phylo.EventLogger(): failed to locate root.")
  
  
  if (transmissions | migrations) {
    # retrieve transmission/migration events
    targets <- c()
    if (migrations) targets <- c(targets, 'migration')
    if (transmissions) targets <- c(targets, 'transmission')
    
    tm.events <- events[!is.na(events$lineage1) & 
                          is.element(events$event.type, targets), ]
    
    for (node in unique(tm.events$lineage1)) {
      # locate origin of this Lineage among core events
      idx <- which(core.events$lineage1==node)
      e <- core.events[idx, ]
      if (nrow(e) > 1) {
        stop("Error in as.phylo.EventLogger: found multiple origins for ",
             "node ", node)
      }
      parent <- e$lineage2
      
      # modify Lineage names in sequence of events
      temp <- tm.events[tm.events$lineage1==node, ]
      temp <- temp[order(temp$time ), ]  # start with most recent
      child <- node
      for (i in 1:nrow(temp)) {
        tm.event <- temp[i, ]
        new.node <- paste(node, i, sep=tm.event$event.type)
        
        tm.event$lineage1 <- child
        tm.event$lineage2 <- new.node
        child <- new.node
        core.events <- rbind(core.events, tm.event)
      }
      
      e$lineage1 <- new.node
      core.events[idx, ] <- e
    }
  }
  
  
  # gather lists of node names
  nodes <- union(core.events$lineage1, core.events$lineage2)
  tips <- unique(core.events$lineage1[
    which(!is.element(core.events$lineage1, core.events$lineage2))
    ])
  internals <- setdiff(nodes, tips)
  

  # populate edge list by postorder traversal (children before parent)
  preorder <- .reorder.inner.events(core.events, root, order='preorder')
  
  # reorder node lists
  tips <- preorder[is.element(preorder, tips)]
  internals <- preorder[is.element(preorder, internals)]
  # numbered 1:n for tips, (n+1):(n+m) for internal nodes, starting with root
  nodes <- c(tips, internals)
  
  # generate edgelist
  edges <- core.events[na.omit(match(preorder, core.events$lineage1)), ]
  edgelist <- as.matrix(t(sapply(1:nrow(edges), function(i) {
    c(which(nodes==edges$lineage2[i]), which(nodes==edges$lineage1[i])) 
  })))  

  
  # calculate edge.lengths
  edge.length <- c()
  for (i in 1:nrow(edges)) {
    e <- edges[i,]
    child <- e$lineage1
    if (is.element(child, tips)) {
      bl <- e$time - as.numeric(fixed.sampl[child])
    }
    else {
      next.t <- unique(edges$time[edges$lineage2 == child])
      if (length(next.t) > 1) {
        stop("Error in as.phylo.EventLogger: more than one time associated with node ", 
             child)
      }
      bl <- e$time - next.t
    }
    edge.length <- c(edge.length, bl)
  }
  
  # create phylo object
  phy <- list(
    tip.label = tips,
    node.label = internals,
    Nnode = length(internals),
    edge = edgelist,
    edge.length = edge.length,
    event.type = edges$event.type,
    compartment1 = edges$compartment1,
    compartment2 = edges$compartment2,
    type1 = edges$type1,
    type2 = edges$type2
    )
  attr(phy, 'class') <- 'phylo'
  phy
}


#' .reorder.inner.events
#' Recursive helper function to generate ordering of events by preorder
#' traversal of inner tree.  Similar to Run::.reorder.events().
#' @param events:  data frame of core (coalescent, bottleneck) events
#' @param node:  character, unique identifier of a node
#' @param result:  vector to pass between recursive calls
#' @keywords internal
.reorder.inner.events <- function(events, node, order, result=c()) {
  if (order=='preorder') {
    result <- c(result, node)  # parent before children
  }
  
  children <- events$lineage1[events$lineage2 == node]
  for (child in children) {
    result <- .reorder.inner.events(events, child, order, result)
  }
  if (order=='postorder') {
    result <- c(result, node)  # children before parent
  }
  return(result)
}







