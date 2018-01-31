
# generate coalescent events

time1 <- function(rate, n_pathogens) {
  # @param rate = rate of coalescence for the particular Compartment
  # @param n_pathogens = number of extant lineages at current draw time
  # @return a randomly generated waiting time (in continuous time) for 2 pathogens to coalesce out of a sample of n_pathogens
  npairs <- (0.5 * n_pathogens * (n_pathogens - 1))  # number of pairings of extant lineages
  rexp(1, rate = (npairs*rate))
}


simulate.coalescent.events <- function(inputs) {
  # Draw a coalescence time, t, from an expoential distributions with a rate of choose(k,2), where k is the number of nodes.
  # Two nodes are selected at random from nodes, and coalesced: a new ancesotr node is created, and the two nodes are added
  # as children to this node with an edge length such that the total length from tip to the ancestral node is equal to the 
  # depth of the deepest child + t. The two nodes are removed from the list of nodes and the new node is added to it.
  # The process repeats until k = 1.
  
  extant <- inputs$get.lineages()
  
  while (length(extant) > 1) {
    # pick two nodes to coalesce at random (sampled without replacement)
    to_coalesce <- sample(extant, 2)
    
    # draw relevant coalescence rate
    lineages_types <- unique(sapply(to_coalesce, function(x){x$get.location()$get.type()$get.name()}))
    if (length(lineages_types) == 1) rate <- inputs$get.types()[[lineages_types]]$get.coalescent.rate()
    
    # draw waiting time to coalescent event
    tau <- time1(rate, length(extant))
    
    # create the new ancestor of these nodes
    # add the nodes as child nodes of the new node
    # set the ancestor's edge length
    
    # add waiting time to extant lineages
    sapply(extant, function(x) {
      update.time <- x$get.sampling.time() + tau
      x$set.sampling.time(update.time)
    })
    
  }
  
}





##########################################################################################################################################

simulate <- function(inputs, eventlog) {
  # Function that performs a coalescent simulation of lineages w/in the host tree, w/ transmission events provided as an EventLogger object
  # udpates event log with coalescent events tracking the lineages tree(s)
  
  # all Lineage sampling times must be stated, either with a vector, or a `numeric`, assumed to be applied to all w/in that Lineage pop'n
  # there must be at least 1 Lineage that starts with a sampling time of 0
  # is it possible that a Lineage's initial `location` is not one of the available transmission tips, and is an internal node of the host tree?
      # if so, then this is a problem. B/c having a tip compartment further down the line implies that there is an associated Lineage
      # earlier than sampling time of 0, which isn't possible
  
  # potential solution: when generating transmission events, can look at all Lineage sampling times and group earliest sampling times
  # but then we would have to check the source compartment's lineages' sampling times every time before assigning a recipient's source compartment
  
  # would this happen in a practical sense anyways?
  # ...actually I don't think so 
}