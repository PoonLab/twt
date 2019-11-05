#' EventLogger
#' 
#' \code{EventLogger} is an R6 class for an object that tracks migration, 
#' transmission, and coalescent events.  Note that bottleneck events are logged as 
#' coalescent events.
#' 
#' @param events: a data frame where each row represents an event with a time 
#' stamp in forward time.
#' @param migration.events.storage: a data frame of migration events, returned by
#' `simMigrations.R:.calc.migration.events()`.
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
    initialize = function(events = NA, migration.events.storage = NA) {
      if (is.na(events)) {
        self$clear.events()
      } else {
        private$events <- events  
      }
      
      if (is.na(migration.events.storage)) {
        private$migration.events.storage <- data.frame(stringsAsFactors = FALSE)
      } else {
        private$migration.events.storage <- migration.events.storage
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
    
    add.event = function(type, time, line1, line2, comp1, comp2) {
      # @param type: event type, one of 'transmission', 'migration', 'coalescence',
      # or 'bottleneck'.
      # @param time: CUMULATIVE time that event has occurred between two compartments 
      # in a transmission/migration/coalescent event
      
      if (is.element(type, c('transmission', 'migration'))) {
        e <- list(event.type=type, time=time, lineage1=NA, lineage2=NA,
                  compartment1=comp1, compartment2=comp2)
      } 
      else if (is.element(type, c('coalescent', 'bottleneck'))) {
        e <- list(event.type=type, time=time, lineage1=line1, lineage2=line2,
                  compartment1=comp1, compartment2=NA)
      } 
      else {
        stop("Error, unrecognized type argument in add.event()")
      }
      
      private$events <- rbind(private$events, e, stringsAsFactors=F)
    }, 
    
    
    clear.events = function() {
      private$events <- data.frame(
        event.type=character(),
        time=numeric(),
        lineage1=character(),
        lineage2=character(),
        compartment1=character(),
        compartment2=character(),
        stringsAsFactors = FALSE
        )
    },
    
    
    modify.event = function(transmission.time, lineages) {
      # when inner tree simulation has reached a transmission event, 
      # need to fill in the lineage column w/ the lineages that are present 
      # at transmission time
      #
      # @param transmission.time:  time of event
      # @param lineages:  
      
      # in the case of bottleneck events, will have same time so have to isolate 
      # transmission event times
      transmission.events <- self$get.events('transmission')
      
      index <- which(transmission.events$time == transmission.time)
      rowname <- rownames(transmission.events)[index]
      eventlog.index <- which(rownames(self$get.all.events()) == rowname)
      private$events[eventlog.index, 'lineage1'] <- lineages
    },
    
    get.migration.events = function() {
      private$migration.events.storage
    },
    
    store.migration.events = function(migration.events) {
      private$migration.events.storage <- migration.events
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
    #events.noncumul = NULL,  ## DEPRECATED (issue #58)
    migration.events.storage = NULL,
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



# Functions to convert the events stored in the EventLogger into an 
# ape::phylo object that can then be plotted or printed.
# These are S3 methods defined outside of the R6 class definition.

#' plot.EventLogger
#' 
#' Generate a plot summarizing the content of an EventLogger object
#' as a phylogenetic tree.
#' 
#' @param eventlog:  Object of class 'EventLogger'
#' @param transmissions:  If TRUE, display transmission events in plot.
#' 
#' @examples 
#' # load model
#' path <- system.file('extdata', 'SI.yaml', package='twt')
#' 
#' # load file and parse to construct MODEL object
#' settings <- yaml.load_file(path)
#' mod <- MODEL$new(settings)
#' 
#' # simulate outer tree
#' outer <- sim.outer.tree(mod)
#' 
#' # simulate inner tree
#' tree <- sim.inner.tree(mod, outer)
#' plot(tree)
#' 
#' @export
plot.EventLogger <- function(eventlog, transmissions=FALSE, migrations=FALSE, 
                             node.labels=FALSE) {
  evt <- eventlog$get.all.events()
  if (is.null(evt)) {
    cat("No events to display.")
  }
  else {
    if (all(is.element(evt$event.type, c('transmission', 'migration')))) {
      # this is an outer tree
      .plot.outer.tree(eventlog)
    } 
    else {
      phy <- .eventlogger.to.phylo(
        eventlog=eventlog, 
        transmissions=transmissions, 
        migrations=migrations, 
        node.labels=node.labels)
      plot(phy)
    }
  }
}


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



#' .plot.outer.tree
#' 
#' A function that converts events stored in an EventLogger object into
#' an outer transmission tree.
#' 
#' @param e: EventLogger object
#' 
#' @return ape::phylo object
#' @keywords internal
.plot.outer.tree <- function(e) {
  events <- e$get.all.events()
  stopifnot(all(events$event.type == 'transmission'))
  
  recip <- events$compartment1
  source <- events$compartment2
  comps <- unique(c(recip, source))
  
  # source that is never a recipient is the root
  root <- unique(source[!is.element(source, recip)])
  if (length(root) != 1) {
    stop("Eventlog contains more than one index compartment (root)")
  }
  
  # TODO work in progress!
  
  # 1. sort compartments in order of infection time
  
  # 2. draw horizontal line segments for compartments
  
  # 3. draw vertical line segments for transmissions
  
  
}


#' .eventlogger.to.phylo
#' 
#' function converts events stored in the EventLogger object into an inner 
#' coalescent tree w/ option to include/exclude transmission events.
#' 
#' @param eventlog: EventLogger object
#' @param transmissions: logical; if TRUE, transmission events included, 
#' else excluded otherwise.
#' @return ape::phylo object
#' @keywords internal
.eventlogger.to.phylo <- function(eventlog, transmissions=FALSE, migrations=FALSE, 
                                  node.labels=FALSE) {
  
  # helper function: if the event being examined is a bottleneck event, 
  # must split column named "$lineage1" into the bottlenecking lineages
  split.bottleneck.lineages <- function(lineages.names) {
    as.vector(strsplit(lineages.names, ',', fixed=T)[[1]])
  }
  
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
  if (is.null(b_events)) {
    b_events_lineages <- NULL
  } else {
    b_events_lineages <- sapply(1:nrow(b_events), function(x) {
      split.bottleneck.lineages(b_events[x,]$lineage1)
    })
  }
  
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
  root <- unlist(
            setdiff(
              c(c_events$compartment1, b_events$compartment1), 
              c(c_events$lineage1, c_events$lineage2, unlist(b_events_lineages))
            )
          )
  if (length(root) > 1) stop('ERROR! Root > 1: ', root)   
  
  tips <- unlist(
            setdiff(
              c(c_events$lineage1, c_events$lineage2, unlist(b_events_lineages)), 
              c(c_events$compartment1, b_events$compartment1)
            )
          )
  
  internals <-  unlist(
    intersect(
      c(c_events$compartment1, b_events$compartment1), 
      c(c_events$lineage1, c_events$lineage2, unlist(b_events_lineages))
      )
    )
  
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
      tip.label[tip.no] <<- node
      
      tip.no <<- tip.no + 1
      
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



