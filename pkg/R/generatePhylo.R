## these functions were extracted from MODEL class definition
## TODO: modify function .outer.tree.to.phylo to accommodate this change
get.leaves.names = function(e) {
  # returns a vector of Compartment object names that are terminal nodes (only a recipient)
  # @param e = EventLogger object
  t_events <- e$get.events('transmission', cumulative=F)
  setdiff(unlist(t_events$compartment1), unlist(t_events$compartment2))
}

get.nonterminals = function(e) {
  # return an iterator over all names of internal nodes of the transmission tree
  # @param e = EventLogger object
  t_events <- e$get.events('transmission', cumulative=F)
  intersect(unlist(t_events$compartment1), unlist(t_events$compartment2))
}

get.node.heights = function() {
  # calculate node heights for all nodes of the tree
  # annotate nodes with heights in place
}



.outer.tree.to.phylo <- function(eventlog) {
  # function converts the transmission events stored in an EventLogger object into a transmission tree
  # @param eventlog = EventLogger object
  # @return phy = ape::phylo object
  
  t_events <- eventlog$get.events('transmission')
  tips <- unlist(setdiff(t_events$compartment1, t_events$compartment2))
  internals <- unlist(intersect(t_events$compartment1, t_events$compartment2))
  root <- unlist(setdiff(t_events$compartment2, t_events$compartment1))
  
  # initialize attributes of an ape::phylo object
  tip.label <- vector()
  edge.length <- vector()
  Nnode <- nrow(t_events)
  node.label <- vector()
  # edge matrix can't be determined as of yet, so we're going to record into a data frame
  edge <- data.frame()
  
  # intiialize indices to be assigned to root, tips, and internals
  tip.no <- 1
  root.no <- length(c(tips, internals, root)) + 1
  node.no <- root.no + 1
  
  
  # helper function (recursive) for STEP 2
  .assign.parent.branch.lengths <- function(node, node.ind, child.length) {
    if (node == root) {
      NULL
    } else {
      node.length <- t_events[which(t_events$compartment1 == node), 'time'] - child.length
      ancestor <- t_events[which(t_events$compartment1 == node), 'compartment2']
      ancestor.ind <- node.no 
      
      node.label[ancestor.ind] <- paste0(ancestor, '->', node)
      edge <- rbind(edge, c(ancestor.ind, node.ind), stringsAsFactors=F)
      
      node.no <- node.no + 1
      .assign.parent.branch.lengths(ancestor, ancestor.ind, node.length)
    }
  }
  
  
  # STEP 1: start at tips and assign branch lengths
  for (t in tips) {
    branch.length <- t_events[which(t_events$compartment1 == t), 'time']
    parent <- t_events[which(t_events$compartment1 == t), 'compartment2']
    parent.ind <- node.no
    
    tip.label[tip.no] <- t
    edge.length[tip.no] <- branch.length
    node.label[parent.ind] <- paste0(parent, '->', t)
    edge <- rbind(edge, c(parent.ind, tip.no), stringsAsFactors=F)
    
    tip.no <- tip.no + 1
    node.no <- node.no + 1
    
    # STEP 2: follow each tip's parents up to root and assign branch lengths
    .assign.parent.branch.lengths(parent, parent.ind, branch.length)
  }
  
  # STEP 3: populate inter-node branch lengths
  for (i in c(root, internals)) {
    specific.events <- t_events[which(t_events$compartment2 == i), ]
    reordered.events <- specific.events[order(specific.events$time),]
    
    if (nrow(reordered.events) > 1) {
      for (pair in 1:nrow(reordered.events)-1) {
        node.below <- reordered.events[pair,]
        node.above <- reordered.events[pair+1,]
        
        node.below.name <- paste0(node.below$compartment2, '->', node.below$compartment1)
        if (node.below.name %in% node.label) {
          node.below.ind <- which(node.label == node.below.name)
        } else {
          node.below.ind <- node.no
          node.label[node.below.ind] <- node.below.name
          node.no <- node.no + 1
        }
        
        node.above.name <- paste0(node.above$compartment2, '->', node.above$compartment1)
        if (node.above.name %in% node.label) {
          node.above.ind <- which(node.label == node.above.name)
        } else {
          node.above.ind <- node.no
          node.label[node.above.ind] <- node.above.name
          node.no <- node.no + 1
        }
        
        edge.length[node.above.ind] <- node.above$time - node.below$time
        edge <- rbind(edge, c(node.above.ind, node.below.ind), stringsAsFactors=F)
          
      }
    } else {
      next
    }
  }
  
  # STEP 4: finally, populate root and internals final branch lengths up to extant time (t=0) as singleton nodes
  for (n in c(root, internals)) {
    singleton.length <- min(t_events[which(t_events$compartment2 == n), 'time'])
    r.name <- t_events[which(t_events$time == singleton.length), 'compartment1']
    node.ind <- which(node.label == paste0(n, '->', r.name))
    
    tip.label[tip.no] <- n
    edge.length[tip.no] <- singleton.length
    edge <- rbind(edge, c(node.ind, tip.no), stringsAsFactors=F)
    
    tip.no <- tip.no + 1
  }
  
  # edge matrix complete, convert data.frame to matrix
  edge.mat <- as.matrix(edge)
  
  phy <- list(tip.label=tip.label, Nnode=Nnode, edge.length=as.numeric(edge.length), edge=edge.mat, node.label=node.label[!is.na(node.label)])
  attr(phy, 'class') <- 'phylo'
  attr(phy, 'order') <- 'cladewise'
  phy
}



plot.EventLogger <- function(eventlog, fixed.samplings=fixed.samplings) {
  phy <- .inner.tree.to.phylo(eventlog=eventlog, fixed.samplings=fixed.samplings)
  plot(phy)
}



.inner.tree.to.phylo <- function(eventlog, fixed.samplings, transmissions=FALSE, migrations=FALSE) {
  # function converts coalescent & migration events stored in the EventLogger object into an inner coalescent tree w/ option to include/exclude transmission events
  # @param eventlog = EventLogger object
  # @param transmissions = logical; if TRUE, transmission events included, else excluded otherwise
  # @return phy = ape::phylo object
  
  if (transmissions) { t_events <- eventlog$get.events('transmission')
  } else {t_events <- NULL}
  
  # if the event being examined is a bottleneck event, must split column named "$lineage1" into the bottlenecking lineages
  split.bottleneck.lineages <- function(lineages.names) {
    as.vector(strsplit(lineages.names, ',', fixed=T)[[1]])
  }
  
  m_events <- eventlog$get.events('migration')
  c_events <- eventlog$get.events('coalescent')
  b_events <- eventlog$get.events('bottleneck')
  b_events_lineages <- sapply(1:nrow(b_events), function(x) {
    split.bottleneck.lineages(b_events[x,]$lineage1)
  })
  
  total_events <- rbind(m_events, c_events, b_events)
  
  # initialize attributes of an ape::phylo object
  tip.label <- vector()
  edge.length <- vector()
  Nnode <- 0
  record.node.label <- vector()   # will keep track of node labels one-to-one with edge matrix, will pare down in function recursive.populate.node.labels()
  node.label <- vector()
  edge <- data.frame()      # eventually will convert into a static edge matrix
  
  # separate nodes into root, tips, and internals
  root <- unlist(setdiff(c_events$compartment1, c(c_events$lineage1, c_events$lineage2, b_events_lineages)))
  tips <- unlist(setdiff(c(c_events$lineage1, c_events$lineage2, b_events_lineages), c(c_events$compartment1, b_events$compartment1)))
  internals <- unlist(intersect(c_events$compartment1, c(c_events$lineage1, c_events$lineage2, b_events_lineages)))
  
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
      individual.branch.length <- fixed.samplings$tip.height[ which(fixed.samplings$tip.label == node) ]
      tip.label[tip.no] <<- node
      
      tip.no <<- tip.no + 1
      
      # return tip's branch length (initial sampling time)
      return (c(individual.branch.length, (tip.no -1)))
      
    } else {
      
      eventRow <- total_events[ which(total_events$compartment1 == node), ]
      
      if (eventRow$event.type == 'bottleneck') {
        children <- split.bottleneck.lineages(eventRow$lineage1)
      } else {
        # this is a coalescent event
        children <- c(eventRow$lineage1, eventRow$lineage2)
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
            
            #node.label[nrow(edge)+1] <<- node
            edge.length[nrow(edge)+1] <<- migration.branch.length
            cat(indiv.node.no, ' ', mig.node.no, ' ', migration.branch.length, '\n')
            edge <<- rbind(edge, c(indiv.node.no, mig.node.no))

            # continue with recursion
            recursive.call <- recursive.populate.branchlength(child)
            child.branch.length <- recursive.call[1]
            child.node.no <- recursive.call[2]
            
            individual.branch.length <- mig.event$time - child.branch.length
            
            record.node.label[nrow(edge)+1] <<- child
            
            edge.length[nrow(edge)+1] <<- individual.branch.length ###
            edge <<- rbind(edge, c(mig.node.no, child.node.no))
            
            #final.node.no <- mig.node.no
            #return (c(eventRow$time, mig.node.no))
            
          } else {
            
            recursive.call <- recursive.populate.branchlength(child)
            child.branch.length <- recursive.call[1]
            child.node.no <- recursive.call[2]
            
            individual.branch.length <- eventRow$time - child.branch.length
            
            record.node.label[nrow(edge)+1] <<- node
            
            edge.length[nrow(edge)+1] <<- individual.branch.length ###
            edge <<- rbind(edge, c(indiv.node.no, child.node.no))
            
          }
          
        } else {
          
          recursive.call <- recursive.populate.branchlength(child)
          child.branch.length <- recursive.call[1]
          child.node.no <- recursive.call[2]
          
          individual.branch.length <- eventRow$time - child.branch.length
          
          record.node.label[nrow(edge)+1] <<- node
          
          edge.length[nrow(edge)+1] <<- individual.branch.length ###
          edge <<- rbind(edge, c(indiv.node.no, child.node.no))

        }
        
      }
      
      # return node's branch length
      return (c(eventRow$time, final.node.no))
      
      
    } 
    
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
    
    #return(node.label)
  }
  
  
  recursive.populate.branchlength(root)

  # convert edge dataframe into matrix
  edge <- as.matrix(edge)
  
  recursive.populate.node.labels(length(tips)+1)
  
  phy <- list(tip.label=tip.label, Nnode=Nnode, edge.length=as.numeric(edge.length), edge=edge.mat, node.label=node.label)
  attr(phy, 'class') <- 'phylo'
  attr(phy, 'order') <- 'cladewise'
  phy
  
}




