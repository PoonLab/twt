# functions to convert the events stored in the EventLogger into an ape::phylo object that can then be plotted or printed

plot.EventLogger <- function(eventlog, transmissions=FALSE, migrations=FALSE, node.labels=FALSE) {
  phy <- .eventlogger.to.phylo(eventlog=eventlog, transmissions=transmissions, migrations=migrations, node.labels=node.labels)
  plot(phy)
}



print.EventLogger <- function(eventlog) {
  eventlog$get.all.events()
}



.eventlogger.to.phylo <- function(eventlog, transmissions=FALSE, migrations=FALSE, node.labels=FALSE) {
  # function converts events stored in the EventLogger object into an inner coalescent tree w/ option to include/exclude transmission events
  # @param eventlog = EventLogger object
  # @param transmissions = logical; if TRUE, transmission events included, else excluded otherwise
  # @return phy = ape::phylo object
  
  # helper function: if the event being examined is a bottleneck event, must split column named "$lineage1" into the bottlenecking lineages
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
                
                t.node.no <- node.no
                node.no <<- node.no + 1
                Nnode <<- Nnode + 1
                
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
              
              t.node.no <- node.no
              node.no <<- node.no + 1
              Nnode <<- Nnode + 1
              
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
  
  
  add.branch <- function(event, node, child, node.no) {
    
    recursive.call <- recursive.populate.branchlength(child)
    child.branch.length <- recursive.call[1]
    child.node.no <- recursive.call[2]
    
    individual.branch.length <- event$time - child.branch.length
    
    record.node.label[nrow(edge)+1] <<- node
    
    edge.length[nrow(edge)+1] <<- individual.branch.length
    edge <<- rbind(edge, c(node.no, child.node.no))
    
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




