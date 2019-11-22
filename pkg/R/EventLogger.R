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
      idx <- which(events$compartment1 == recipient$get.name() && 
                     events$compartment2 == source$get.name() &&
                     events$time == time &&
                     is.na(events$lineage1))
      if (length(idx) != 1) {
        stop(paste("Error in eventLogger:record.migration - ",
                   ifelse(length(idx)>1, "multiple", "no"), 
                   " transmission events match recipient ",
                   recipient$get.name()))
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
#' @export
as.phylo.EventLogger <- function(eventlog, transmissions=FALSE, migrations=FALSE, 
                                  node.labels=FALSE) {
  
  # fixed sampling times of tips and heights
  fixed.sampl <- eventlog$get.fixed.samplings()
  
  # transmission events, can be optionally included
  if (transmissions) {
    t_events <- eventlog$get.events('transmission')
  } else {
    t_events <- NULL
  }
  
  # migration events, can be optionally included
  if (migrations) {
    m_events <- eventlog$get.events('migration')
  } else {
    m_events <- NULL
  }
  
  # bottleneck events, may or may not be present
  b_events <- eventlog$get.events('bottleneck')
  
  # coalescent events
  c_events <- eventlog$get.events('coalescent')
  
  # minimum set of events to be included
  core_events <- rbind(c_events, b_events)
  
  # initialize attributes of an ape::phylo object
  tip.label <- vector()
  edge.length <- vector()
  Nnode <- 0
  record.node.label <- vector()   # will keep track of node labels one-to-one with edge matrix, will pare down in function recursive.populate.node.labels()
  node.label <- vector()
  edge <- data.frame()      # eventually will convert into a static edge matrix

    
  # separate nodes into root, tips, and internals
  # note lineage1 is child, lineage2 is parent
  root <- unique(core_events$lineage2[
    which(!is.element(core_events$lineage2, core_events$lineage1))
    ])
  
  if (length(root) > 1) stop('ERROR! Found multiple root nodes: ', root)
  if (length(root) == 0) stop("Error in as.phylo.EventLogger(): failed to locate root.")
  
  nodes <- union(core_events$lineage1, core_events$lineage2)
  
  tips <- unique(core_events$lineage1[
    which(!is.element(core_events$lineage1, core_events$lineage2))
  ])
  
  internals <- setdiff(nodes, tips)
  
  
  # initialize ape::phylo indices to be assigned to root, tips, and internals
  tip.no <- 1
  node.no <- length(tips) + 1
  

  recursive.populate.branchlength <- function(node) {
    # postorder traversal to populate ape::phylo object
    # visit and populate branch lengths of all children before visiting/populating itself
    # @param node = name of a node of type character()
    # @return branch.length = branch length of node of type numeric()
    
    if (node %in% tips) {
      
      # add node.sampling.time as branch length
      individual.branch.length <- fixed.sampl$tip.height[ which(fixed.sampl$tip.label == node) ]
      tip.label[tip.no] <- node
      
      tip.no <- tip.no + 1
      
      # return tip's branch length (initial sampling time)
      return (c(individual.branch.length, (tip.no -1)))
      
    } else {
      
      eventRow <- core_events[ which(core_events$compartment1 == node), ]
      
      if (eventRow$event.type == 'coalescent') {
        # this is a coalescent event
        children <- c(eventRow$lineage1, eventRow$lineage2)
      } else {
        # this is a bottleneck event
        children <- split.bottleneck.lineages(eventRow$lineage1)
      }

      # record current node.no before recursive call
      indiv.node.no <- node.no
      
      node.no <<- node.no + 1
      Nnode <<- Nnode + 1
      
      final.node.no <- indiv.node.no
      
      for (child in children) {
          
        if (migrations) {
          # check for migration events involving `node`; must incorporate singleton node
          if (child %in% m_events$lineage1) {
            
            mig.event <- m_events[ which(m_events$lineage1 == child), ]
            
            mig.node.no <- node.no
            node.no <<- node.no + 1
            Nnode <<- Nnode + 1
            
            # create singleton node here
            migration.branch.length <- eventRow$time - mig.event$time
            
            edge.length[nrow(edge)+1] <<- migration.branch.length
            edge <<- rbind(edge, c(indiv.node.no, mig.node.no))

            # continue with recursion
            add.branch(mig.event, child, child, mig.node.no)
            
          } else {
            
            if (transmissions) {
              if (child %in% t_events$lineage1) {
                
                t.event <- t_events[ which(t_events$lineage1 == child), ]
                t_events <- t_events[-which(t_events$lineage1 == child), ]
                
                t.node.no <- node.no
                node.no <<- node.no + 1
                Nnode <<- Nnode + 1
                
                # create a singleton node here
                transmission.branch.length <- eventRow$time - t.event$time
                
                edge.length[nrow(edge)+1] <<- transmission.branch.length
                edge <<- rbind(edge, c(indiv.node.no, t.node.no))
                
                # continue with recursion
                add.branch(t.event, child, child, t.node.no)
                
              } else {
                add.branch(eventRow, node, child, indiv.node.no)
              }
              
            } else {
              add.branch(eventRow, node, child, indiv.node.no)
            }
          
          }
          
        } else {
          
          if (transmissions) {
            if (child %in% t_events$lineage1) {
              
              t.event <- t_events[ which(t_events$lineage1 == child), ]
              t_events <- t_events[-which(t_events$lineage1 == child), ]
              
              t.node.no <- node.no
              node.no <<- node.no + 1
              Nnode <<- Nnode + 1
              
              # create a singleton node here
              transmission.branch.length <- eventRow$time - t.event$time
              
              edge.length[nrow(edge)+1] <<- transmission.branch.length
              edge <<- rbind(edge, c(indiv.node.no, t.node.no))
              
              # continue with recursion
              add.branch(t.event, child, child, t.node.no)
              
            } else {
              add.branch(eventRow, node, child, indiv.node.no)
            }
            
          } else {
            add.branch(eventRow, node, child, indiv.node.no)
          }

        }
        
      }
      
      # return node's branch length
      return (c(eventRow$time, final.node.no))
      
      
    } 
      
  }
  
  
  add.branch <- function(event, parent.node, child.node, parent.node.no) {
    
    recursive.call <- recursive.populate.branchlength(child.node)
    child.branch.length <- recursive.call[1]
    child.node.no <- recursive.call[2]
    
    individual.branch.length <- event$time - child.branch.length
    
    record.node.label[nrow(edge)+1] <<- parent.node
    
    edge.length[nrow(edge)+1] <<- individual.branch.length
    edge <<- rbind(edge, c(parent.node.no, child.node.no))
    
  }
  

  recursive.populate.node.labels <- function(node) {
    # ape::phylo populates node labels in plot function by preorder traversal
    # visit node, then go right, then go left
    
    if (node %in% 1:length(tips) == FALSE) {
      label <- record.node.label[which(edge[,1] == node)[1]]
      if (label %in% node.label == FALSE) {
        # record name in node.label vector
        node.label[length(node.label)+1] <<- label
        
        # access children, right first, then left
        children <- edge[which(edge[,1] == node), 2]
        iterating.list <- children[order(children)]
        for (child in iterating.list) {
          recursive.populate.node.labels(child)
        }
      }
    }
  
  }
  
  
  if (transmissions) {
    
    if (root %in% t_events$lineage1) {
      t.event <- t_events[ which(t_events$lineage1 == root), ]
      
      t.node.no <- node.no
      node.no <- node.no + 1
      Nnode <- Nnode + 1
      
      # now begin recursion
      add.branch(t.event, root, root, t.node.no)
    } else {
      recursive.populate.branchlength(root)
    }
    
  } else {
    recursive.populate.branchlength(root)
  }

  
  # convert edge dataframe into matrix
  edge <- as.matrix(edge)
  
  recursive.populate.node.labels(length(tips)+1)
  
  if (node.labels) {
    phy <- list(tip.label=tip.label, Nnode=Nnode, edge.length=as.numeric(edge.length), edge=edge, node.label=node.label)
  } else {
    # default behaviour: no node labels bc INDELible doesn't like them
    phy <- list(tip.label=tip.label, Nnode=Nnode, edge.length=as.numeric(edge.length), edge=edge)
  }
  
  if (length(phy$edge.length) != nrow(phy$edge)) {
    cat('Nrow does not equal edge length.')
  }
  
  attr(phy, 'class') <- 'phylo'
  attr(phy, 'order') <- 'cladewise'
  phy
}


