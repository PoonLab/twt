.to.transmission.tree <- function(eventlog) {
  # function converts transmission events stored in an EventLogger object into a transmission tree
  # @param eventlog = EventLogger object
  # @return phy = ape::phylo object
  
  t_events <- eventlog$get.events('transmission')
  root <- unlist(setdiff(t_events$compartment2, t_events$compartment1))
  tips <- unlist(setdiff(t_events$compartment1, t_events$compartment2))
  
  storage <- data.frame()
  
  # initializing attributes of an ape::phylo object
  tip.label <- vector()
  edge.length <- vector()
  Nnode <- nrow(t_events)
  edge <- matrix(nrow=nrow(t_events)*2, ncol=2)
  node.label <- vector()
  
  # initialize index assignment numbers
  tip.no <- 1
  root.no <- length(tips) + 1
  node.no <- root.no + 1
  edge.count <- 1
  
  root.source.events <- t_events[ which(t_events$compartment2 == root), ]
  root.event.ind <- which(t_events$time == max(root.source.events$time))
  
  for (x in 1:nrow(t_events)) {
    source <- t_events[x,]$compartment2
    recipient <- t_events[x,]$compartment1
    
    # arbitrarily assign index if recipient not a tip
    if (recipient %in% tips) {
      recipient.ind <- tip.no
      tip.label[recipient.ind] <- recipient
      tip.no <- tip.no + 1
    } else {
      recipient.ind <- node.no
      node.label[recipient.ind] <- recipient
      node.no <- node.no + 1
    }
    
    # arbitrarily assign index if not specifically the root event
    if (x == root.event.ind) {
      source.ind <- root.no
      node.label[source.ind] <- source
    } else {
      source.ind <- node.no
      node.label[source.ind] <- source
      node.no <- node.no + 1
    }
    
    if (source %in% storage$compartment1) {        # if a source in an earlier recorded event as a recipient
      storage.indices <- storage[which(storage$compartment1 == source), ]
      s.ind <- storage.indices$s.ind
      edge[edge.count,] <- as.numeric(c(s.ind, source.ind))
      edge.length[edge.count] <- storage.ind$time - t_events[x,]$time
      edge.count <- edge.count + 1
    }
    
    if (source %in% storage$compartment2) {        # if a source in earlier recorded event(s) as a source
      storage.indices <- storage[which(storage$compartment2 == source), ]
      earlier.times <- which(sapply(1:nrow(storage.indices), function(y) storage[y,]$time > t_events[x,]$time ) == T)
      later.times <- storage.indices[-earlier.times, ]
      
      closest.earlier.time <- min(sapply(earlier.times, function(z) storage[z,]$time ))
      cet.ind <- which(sapply(earlier.times, function(z) storage[z,]$time) == closest.earlier.time)
      in.between.time <- closest.earlier.time - t_events[x,]$time
      edge[edge.count,] <- as.numeric(c(storage.indices[cet.ind,]$r.ind, source.ind))
      edge.length[edge.count] <- in.between.time
      edge.count <- edge.count + 1
      
      closest.later.time <- max(sapply(later.times, function(z) storage[z,]$time ))
      clt.ind <- which(sapply(later.times, function(z) storage[z,]$time) == closest.later.time)
      modify.time <- t_events[x,]$time - closest.later.time
      edge.to.modify <- which(edge[,2] == clt.ind)
      edge[edge.to.modify,] <- as.numeric(c(source.ind, clt.ind))
      edge.length[edge.to.modify] <- modify.time
    }
    
    if (recipient %in% storage$compartment2) {        # if a recipient in an earlier recorded event as source(s)
      storage.ind <- storage[which(storage$compartment2 == recipient), 'r.ind']
      for (s in storage.ind) {
        edge[edge.count,] <- as.numeric(c(s, recipient.ind))
        edge.length[edge.count] <- t_events[x,]$time
        edge.count <- edge.count + 1
      }
    }
    
    store.event <- c(t_events[x,], r.ind=recipient.ind, s.ind=source.ind)
    storage <- rbind(storage, store.event, stringsAsFactors=F) 
  }
  
  phy <- list(tip.label=tip.label, Nnode=Nnode, edge.length=as.numeric(edge.length), edge=edge, node.label=node.label[!is.na(node.label)])
  attr(phy, 'class') <- 'phylo'
  attr(phy, 'order') <- 'cladewise'
  phy
  
}