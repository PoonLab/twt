## these functions were extracted from MODEL class definition
## TODO: modify function .outer.tree.to.phylo to accomodate this change
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




plot.EventLogger <- function(eventlog) {
  # function plots the population trajectories of susceptibles and infected over time
  t_events <- eventlog$get.events('transmission')
  
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





.inner.tree.to.phylo <- function(eventlog, transmissions=FALSE) {
  # function converts the coalescent events stored in the EventLogger object into an inner coalescent tree with an option to include or exclude transmission events
  # TODO: make this functional to include migration events as well
  # @param eventlog = EventLogger object
  # @param transmissions = logical; if TRUE, transmission events are included, with transmission events being excluded otherwise
  # @return phy = ape::phylo object
  
  if (transmissions) {
    t_events <- eventlog$get.events('transmission')
  } else {
    t_events <- NULL
  }
  c_events <- eventlog$get.events('coalescent')         # TODO: will need to include migration events here later
  
  tips <- unlist(setdiff(union(c_events$lineage1, c_events$lineage2), c_events$compartment1))      # determined by which lineages are not in the ancestral (`compartment1`) column
  root <- unlist(setdiff(c_events$compartment1, union(c_events$lineage1, c_events$lineage2)))      # determined by which lineage is not present in a coalescent event (`lineage1` and `lineage2`)
  internals <- unlist(intersect(c_events$compartment1, union(c_events$lineage1, c_events$lineage2)))
  
  # initialize attributes of an ape::phylo object
  tip.label <- vector()
  edge.length <- vector()
  # Nnode --> can't be determined as of yet
  node.label <- vector()
  # edge matrix can't be determined as of yet, so recording in a data frame
  edge <- data.frame()
  
  # intialize indices to be assigned to root, tips, and internals
  tip.no <- 1 
  root.no <- length(tips) + 1
  node.no <- root.no + 1
  
  
  # helper function (recursive) for STEP 2
  
  
  # STEP 1: start at tips and assign branch lengths
  
  # STEP 2a: check if parents of tips/nodes already included in `tip.label` or `edge.label` to minimize redundancy
  
  # STEP 2b: follow parents of tips up to root and assign branch lengths
  
  # STEP 3: transmission events are treated like any other internal node, save for the edge label should reflect the transmission event (ie. 'A->B')
  
  
}





