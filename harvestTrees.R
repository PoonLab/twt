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
    unsampled.hosts = NULL,
    compartments = NULL,
    lineages = NULL,

    extant = NULL,
    choices = NULL,
    
    initialize = function(settings=NA) {
      private$load.types(settings)
      private$load.unsampled.hosts()
      private$load.compartments(settings)
      private$set.sources()
      private$load.lineages(settings)
      # TODO: populate extant with extant lineages
    },
    
    get.types = function() {self$types},
    get.unsampled.hosts = function() {self$unsampled.hosts},
    get.compartments = function() {self$compartments},
    get.lineages = function() {self$lineages},
    
    
    ## collect host locations of all extant lineages into named list of host1:[pathogen1, pathogen2, ...] name-value pairs
    get.locations = function() {
      locations <- list()  # reset the list
      for (node in self$extant) {    # sapply does not work here... will create list of `$1.host`, `$2.host`, etc
        # TODO: check that lineage is extant
        my.comp.type <- node$get.location()$get.type()$get.name()
        # append this lineage to the vector
        locations[[my.comp.type]] <- c(locations[[my.comp.type]], node)
      }
      locations
    },
    
    
    ## function to extract all pairs of lineages that may coalesce
    get.pairs = function() {
      locations <- self$get.locations()
      choices <- list()
      sapply(locations, function(x) {
        if (length(x) > 1) {
          pairs <- t(combn(1:length(x), 2))
          for (row in 1:nrow(pairs)) {
            pair <- pairs[row,]
            choices[[pair]] <- c(choices[[pair]], x)    # update list of pathogen pairs in same host
          } # TODO: store not in tuples, but in some other data structure
        }
      })
      self$choices <- choices
    }
    
    
  ),
  
  
  private = list(
    
    ## function creates CompartmentType objects
    ## within each CompartmentType, there are distinct compartments with individual transmission & migration rates, and unsampled host & susceptible populations 
    load.types = function(settings) {
      types <- sapply(names(settings$CompartmentType), function(x) {
        params <- settings$CompartmentType[[x]]
        x <- CompartmentType$new(name = x,
                                 no.unsampled = eval(parse(text=paste('list', params$no.unsampled))),
                                 no.susceptible = eval(parse(text=paste('list', params$no.susceptible))),
                                 transmission.rates = eval(parse(text=paste('list', params$transmission.rates))),
                                 migration.rates = eval(parse(text=paste('list', params$migration.rates))),
                                 coalescent.rate = params$coalescent.rate,
                                 bottleneck.size = params$bottleneck.size
        )
      })
      self$types <- types
    },
    
    
    ## function creates "blank" Compartment objects for Unsampled Hosts (US) 
    ## stored in lists for each section within a CompartmentType object
    load.unsampled.hosts = function() {
      us.hosts <- sapply(self$types, function(x) {
        types.unsampled <- names(x$no.unsampled)           # accessing a private variable here; maybe add another public method in CompartmentType instead
        indiv <- lapply(types.unsampled, function(y) {
          compartY <- list()
          nBlanks <- x$get.no.unsampled(y)
          for(blank in 1:nBlanks) {
            x <- Compartment$new()
            compartY[[blank]] <- x
          }
          compartY
        })
        indiv
      })
      self$unsampled.hosts <- us.hosts
    },
    
    
    ## function creates Compartment objects
    ## `type` attr points directly back to a CompartmentType object, and `name` attr is a unique identifier
    load.compartments = function(settings) {
      compartments <- sapply(names(settings$Compartments), function(comp) {
        compartX <- list()
        params <- settings$Compartments[[comp]]
        if (params$type %in% names(self$types)) {
          typeObj <- self$types[[ which(names(self$types) == params$type) ]]    # pointer to CompartmentType object
        } else {
          stop(params$type, ' of Compartment ', comp, ' is not a specified Compartment Type object')
        }
        nIndiv <- params$pop.size
        for(obj in 1:nIndiv) {
          x <- Compartment$new(name = paste0(comp,'_', obj),                     # unique identifier
                               type = typeObj,
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
    
    ## re-iterates over generated Compartment objects and populates `source` attr with R6 objects
    ## sets 'pointers' to other Compartment objects after all have been generated with function call private$load.compartments()
    set.sources = function() {
      compNames <- sapply(self$compartments, function(n){n$get.name()})
      compartments <- sapply(self$compartments, function(x) {
        if (paste0(x$get.source(),'_1') %in% compNames) {                                        # FIXME: arbitrary assignment of source
          sourceObj <- self$compartments[[ which(compNames == paste0(x$get.source(),'_1')) ]]    # FIXME: arbitrarily assigning source to first object in Compartment with $name == x$get.source()
          x$set.source(sourceObj)
        } # TODO: else statement { if source is 'undefined' or not in the list, must be assigned to an unsampled host (US) }
        x
      })
      self$compartments <- compartments
    },
    
    
    ## function creates Lineage objects
    ## `location` attr points directly to a Compartment object, and `name` attr is unique identifier
    load.lineages = function(settings) {
      lineages <- sapply(names(settings$Lineages), function(label) {
        lineageX <- list()
        params <- settings$Lineages[[label]]
        
        if (params$location %in% names(self$compartments)) {                            # set 'pointer' to Compartment object for location
          locationObj <- self$compartments[[ which( sapply(self$compartments, function(y){which(y$get.name() == paste0(params$location,'_1'))}) ==1) ]]  # FIXME: arbitrarily assigns location to first object in Compartment with $name == params$location
        } else {
          stop(params$location, ' of Lineage ', label, ' is not a specified Compartment object')
        }
        nIndiv <- params$pop.size
        if (is.character(params$sampling.time)) {
          vec <- unlist(strsplit(params$sampling.time, split='[`(`|,|`)`]'))
          sampleTimes <- as.double(vec[nzchar(x=vec)])
          if (length(sampleTimes) != nIndiv) { stop(paste('Lineage', label, 'sampling.time does not match pop.size specified for respective Lineage.'))}
        } else {
          sampleTimes <- params$sampling.time
        }
        for (obj in 1:nIndiv) {
          x <- Lineage$new(name = paste0(label,'_',obj),                                 # unique identifier
                           type = params$type,
                           sampling.time = sampleTimes[obj],
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




