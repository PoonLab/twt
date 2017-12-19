## test script
require(R6)
require(yaml)
setwd('~/git/treeswithintrees')
settings <- yaml.load_file('example1.yaml')
test <- NestedCoalescent$new(settings)



# Load all of the different objects into one larger class
NestedCoalescent <- R6Class("NestedCoalescent",
  public = list(
    settings = NULL,
    types = NULL,
    compartments = NULL,
    lineages = NULL,

    choices = NULL,
    
    initialize = function(settings=NA) {
      private$load.types(settings)
      private$load.compartments(settings)
      private$set.sources()
      private$load.lineages(settings)
    },
    
    get.types = function() {self$types},
    get.compartments = function() {self$compartments},
    get.lineages = function() {self$lineages},
    
    
    ## function to extract all pairs of lineages that may coalesce
    get.pairs = function() {
      locations <- private$get.locations()
      self$choices <- list()
      sapply(locations, function(x) {
        if (length(x) > 1) {
          for (pair in combn(x, 2)) {
            self$choices[[pair]] <- c(self$choices[[pair]], x)    # update list of pathogen pairs in same host
          }
        }
      })
    }
    
  ),
  private = list(
    locations = NULL,
    
    
    ## collect host locations of all extant lineages into named list of host1:[pathogen1, pathogen2, ...] name-value pairs
    get.locations = function() {
      private$locations <- list()  # reset the list
      sapply(self$lineages, function(node) {    # should be lineages which are extant during this cycle  <- self$lineages might not cut it
        # TODO: check that lineage is extant
        my.comp.type <- node$get.location()$get.type()
        # append this lineage to the vector
        private$locations[[my.comp.type]] <- c(private$locations[[my.comp.type]], node)
      })
      private$locations
    },
    
    
    load.types = function(settings) {
      types <- sapply(names(settings$CompartmentType), function(x) {
        params <- settings$CompartmentType[[x]]
        x <- CompartmentType$new(name = x,
                                 transmission.rates = params$transmission.rates,
                                 migration.rates = params$migration.rates,
                                 coalescent.rate = params$coalescent.rate,
                                 bottleneck.size = params$bottleneck.size
        )
      })
      self$types <- types
    },
    
    
    load.compartments = function(settings) {
      compartments <- sapply(names(settings$Compartments), function(x) {
        compartX <- list()
        params <- settings$Compartments[[x]]
        # set 'pointer' to CompartmentType object for type
        if (params$type %in% names(self$types)) {
          typeObj <- self$types[[ which(names(self$types) == params$type) ]]
        } else {
          stop(params$type, ' of Compartment ', x, ' is not a specified Compartment Type object')
        }
        nIndiv <- params$pop.size
        for(obj in 1:nIndiv) {
          x <- Compartment$new(type = typeObj,
                               source = params$source,        
                               inf.time = params$inf.time,
                               sampling.time = params$sampling.time
          )
          compartX[[obj]] <- x
        }
        compartX
      })
      self$compartments <- compartments
    },
    
    # set 'pointers' to other Compartment objects, after all compartments have been generated with private$load.compartments()
    set.sources = function() {
      compartments <- sapply(self$compartments, function(x) {
        if (x$get.source() %in% names(self$compartments)) {
          sourceObj <- self$compartments[[ which(names(self$compartments) == x$get.source()) ]]
          x$set.source(sourceObj)
        }
        x
      })
      self$compartments <- compartments
    },
    
    
    load.lineages = function(settings) {
      lineages <- sapply(names(settings$Lineages), function(x) {
        lineageX <- list()
        params <- settings$Lineages[[x]]
        # set 'pointer' to Compartment object for location
        if (params$location %in% names(self$compartments)) {
          locationObj <- self$compartments[[ which(names(self$compartments) == params$location) ]]
        } else {
          stop(params$location, ' of Lineage ', x, ' is not a specified Compartment object')
        }
        nIndiv <- params$pop.size
        for (obj in 1:nIndiv) {
          x <- Lineage$new(type = params$type,
                           sampling.time = params$sampling.time,
                           location = locationObj
          )
          lineageX[[obj]] <- x
        }
        lineageX
      })
      self$lineages <- lineages
    }
    
  )
)






# simulate a (bifurcating) tree under the coalescent model
# end goal --> ape::phylo object with the following attributes:
  # edge = 1x2 matrix int [1:ntips*2-2, 1:2] 
  # Nnode = ntips - 1
  # tip.label = chr [1:ntips] 'A' 'B' etc
  # edge.length = branch lengths (indices correspond to tip.label or node.label)
  # attr(*, 'class') = chr 'phylo'
  # attr(*, 'order') = chr 'cladewise'
simulate.tree <- function(test, simulationTime) {
  require(ape)
  
  ntips <- length(test$get.lineages())
  
  edge <- matrix(nrow=(ntips*2-2), ncol=2)
  Nnode <- ntips - 1
  tip.label <- names(lineages)
  edge.length <- numeric(ntips*2-2)
  
  extant <- 1:ntips                   # list of nodes that have not yet coalesced
  nodes <- ntips+1 : ntips*2-1        # list of ancestors to sample from
  
  while(length(extant) > 1) {
    # create a list of all possible combinations, and then eliminate the ones that cannot coalesce
    #  * lineages cannot coalesce before the sources coalesce -- is this "looking forward" in time?
    all.pairs <- t(combn(extant, 2))
    choices <- list()
    ind <- 1
    for (row in 1:nrow(all.pairs)) {
      pair <- all.pairs[row,]
      sources <- sapply(pair, function(x) {get.source(compartments, lineages, x)})
      # if the infection times of any of the remaining compartments precedes the infection times of the pair's sources, 
        # this means the tips haven't coalesced yet, so the pair isn't possible
      pair.inf.times <- sapply(sources, function(x) {x$inf.time})
      other.inf.times <- sapply(setdiff(compartments, sources), function(x) {min(pair.inf.times) < x$inf.time && max(pair.inf.times) > x$inf.time})
      # if not 0, then at least one of the remaining compartments have an infection time earlier than those of the pair chosen
      if (sum(other.inf.times) == 0) {    
        choices[[ind]] <- pair
        ind <- ind + 1
      }
    }
    
    # pick two nodes to coalesce at random (sampled without replacement)
    to.coalesce <- sample(choices, 1)
    
    # draw relevant transmission rate
    lineageName <- lineages[[to.coalesce[1]]]
    compartmentName <- lineageName$location
    compartmentType <- compartments[[compartmentName]]$type
    rate <- compartmentTypes[[compartmentType]]$coalescent.rate
      
    # draw waiting time to coalescent event
    tau <- time1(rate, length(extant))
    
    # create the new ancestor of these nodes
    new.ancestor <- sample(nodes, 1)
    
    # set new ancestor to children, and add children to ancestor
    edge.ind <- 1
    edge <- sapply(to.coalesce, function(x) {
      edge[edge.ind, 1] <- new.ancestor
      edge[edge.ind, 2] <- x
      edge.ind <- edge.ind + 1
    })
    
    # set branch length of children and for the rest of the extant lineages, add waiting time (will be cumulative as a node remains in extant list)
    # branch length of ancestor will be set when it is picked later as a node to coalesce
    # "a new ancestor node is created, and the two nodes are added as childrent to this node with an edge length such the total length from tip to the
    # ancestral node is equal to the depth of the deepest child + tau"
    edge.length <- sapply()
    
    
    
    # remove children and add ancestor to list of extant lineages
    extant <- union( setdiff(extant, to.coalesce), new.ancestor)
    # remove chosen ancestor node from list of available nodes
    nodes <- setdiff(nodes, new.ancestor)
    
  }
  
  tree <- list(edge=edge, Nnode=Nnode, tip.label=tip.label, edge.length=edge.length)
  class(tree) <- 'phylo'
  order(tree) <- 'cladewise'
  tree
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



# need to calculate node height for a given node of the tree and then annotate heights for the node in place
.generate.nodeheight <- function(ntips, node, children, simulationTime) {
   
}


# look into .coalesce.lineages() function in Kaphi
.coalesce.paths <- function(n1, n2) {
  
}




