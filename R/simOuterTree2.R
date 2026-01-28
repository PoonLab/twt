#' sim.outer.tree
#'
#' Starting from the most recently sampled host, we proceed backwards in 
#' simulation time and follow 
#' 
#' @param mod:  R6 object of class Model
#' @param eventlog:  character or data frame of event log from forward-time 
#'                   simulation of population trajectories
#' @export
sim.outer.tree <- function(mod, eventlog, chunk.size=100) {
  # container for recording outer events
  outer.log <- data.frame(
    time=numeric(),
    event=character(),
    host.anc=character(),
    host.des=character(),
    comp.anc=character(),
    comp.des=character()
  )
    
  if (is.character(eventlog) & file.exists(eventlog)) {
    # load logfile into memory
    eventlog <- read.csv(eventlog)
  }
  counts <- get.counts(eventlog, mod)
  stopifnot(all(counts$time == c(0, eventlog$time)))
  
  cnames <- mod$get.compartments()
  k <- length(cnames)
  stopifnot(all(names(counts)[2:ncol(counts)] == cnames))
  
  # merge event log and counts
  for (cn in cnames) {
    eventlog[[cn]] <- counts[[cn]][2:nrow(counts)]
  }
  
  sampling <- mod$get.sampling()
  if (sampling$mode != "compartment") {
    stop("Unsupported sampling mode", sampling$mode)
  }
  
  # confirm sufficient numbers in sampled compartments
  for (cn in names(sampling$targets)) {
    size <- sum(eventlog$event=="migration" & eventlog$dest==cn)
    target <- sampling$targets[[cn]]
    if (target < size) {
      stop("Insufficient number in sampling compartment ", 
           cn, " (", size , "<", target, ")")
    }
  }
  nsamples <- sum(unlist(targets))
  
  active <- HostSet$new()  # store active Host lineages
  
  # iterate through events in reverse (start with most recent)
  for (row in seq(nrow(eventlog), 1, -1)) {
    # check stopping criterion
    if (nrow(sampled)==nsamples &  # sampled all hosts
        sum(active) == 1) {  # reached root of transmission tree
      break  
    }
    
    e <- eventlog[row, ]
    
    if (e$event == "migration") {
      
      if (e$dest %in% names(targets)) {
        # this was a sampling event
        
        if (nrow(sampled) < nsamples) {
          host <- Host$new(compartment=e$src, sampling.time=e$time)
          active$add.host(host)          
        }

      } else {  
        # this was NOT a sampling event
        
        # check if an active lineage underwent a migration
        active.in.comp <- active$count_type(e$dest)
        if (active.in.comp > 0) {
          tot <- e[[e$dest]]  # size of compartment AFTER migration
          prob <- active.in.comp / tot
          if (runif(1) < prob) {
            
            # choose an eligible Host at random to migrate
            host <- active$sample.host(e$dest)
            host$set.compartment(e$src)
            
            # record this event
            outer.log[nrow(outer.log)+1, ] <- 
              c(time=e$time, event=e$event, host.dec=host$get.name(), 
                comp.anc=e$src, comp.dec=e$dest)
          }
          # otherwise ignore this migration
        }
      }
    }
    if (e$event == "transmission") {
      
      # does this involve an active lineage?
      
    }
  }
}