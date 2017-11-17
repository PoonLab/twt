# simulate a tree under the coalescent model
# end goal --> ape::phylo object with the following attributes:
  # edge = 1x2 matrix int [1:ntips*2-1, 1:2] 
  # Nnode = ntips - 1
  # tip.label = chr [1:ntips] 'A' 'B' etc
  # edge.length = branch lengths (indices correspond to tip.label or node.label)
  # attr(*, 'class') = chr 'phylo'
  # attr(*, 'order') = chr 'cladewise'
harvest.tree <- function(ntips, transmission, migration, compartments, simulationTime) {
  extant <- 1:ntips                   # list of nodes that have not yet coalesced
  nodes <- ntips+1 : ntips*2-1
  
}


# need to calculate node height for a given node of the tree and then annotate heights for the node in place
.generate.nodeheight <- function(ntips, node, children, simulationTime) {
   
}


# look into .coalesce.lineages() function in Kaphi
.coalesce.paths <- function(n1, n2) {
  
}


# modify ancestors of given child
.set.parents <- function(child, parent) {
  
}


# modify children of ancestors
.set.children <- function(parent, children) {
  
}



# returns a single waiting time for the coalescence of two pathogen lineages w/ an exponential distribution
# randomly generated waiting time (in continuous time) for 2 pathogens to coalesce out of a sample of n_pathogens
time1 <- function(rate, ntips) {
  k <- choose(ntips, 2)
  return( rexp(n=1, rate=(k*rate)) )
}


# coalescent interval = time interval btwn two adjacent coalescent events
# the length of each interval is an exponential random variable
k <- choose(ntips, 2)          # number of pairings of extant lineages
rate <- 20                     # where is this defined?
# returns a generator that yields (n - 1) waiting times; each pass through the loop deals with one coalescent interval
#accessible through nextElem(time1)
timer <- iter( sapply(1:(ntips-1), function(x) {rexp(n=1, rate = (k*rate))}) )   


while(length(extant) > 1) {
  # draw waiting time to coalescent event
  tau <- time1(rate, length(extant))
  
  # pick two nodes to coalesce at random (sampled without replacement)
  to.coalesce <- sample(extant, 2)
  
  # create the new ancestor of these nodes
  new.ancestor <- sample(nodes, 1)
  
  
  
  # add waiting time to extant lineages
  
}
