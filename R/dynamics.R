#' Simulate population trajectories for all compartments forward in time
#' @param mod:  R6 object of class Model
simulate.dynamics <- function(mod) {
  # unpack the Model object
  params <- mod$get.parameters()
  cnames <- mod$get.compartments()
  k <- length(cnames)
  sampling <- mod$get.sampling()
  
  # instantiate model parameters
  e <- new.env()
  for (key in names(params)) {
    eval(parse(text=paste(key, "<-", params[[key]])), envir=e)
  }
  
  # instantiate model variables (compartments)
  counts <- mod$get.init.sizes()
  for (cn in cnames) {
    eval(parse(text=paste(cn, "<-", counts[[cn]])), envir=e)
  }

  # instantiate rate matrices
  birth.rates <- sapply(mod$get.birth.rates(), function(x) {
    eval(parse(text=x), envir=e)
  })
  death.rates <- sapply(mod$get.death.rates(), function(x) {
    eval(parse(text=x), envir=e)
  })
  migration.rates <- apply(mod$get.migration.rates(), MARGIN=c(1,2), 
                           function(x) { eval(parse(text=x), envir=e) })
  transmission.rates <- apply(mod$get.transmission.rates(), MARGIN=c(1,2),
                              function(x) { eval(parse(text=x), envir=e) })
  
  # are there other stopping criteria?
  targets <- NULL
  if (sampling$mode == "compartment") {
    targets <- sampling$targets
  }
  
  # main simulation loop
  cur.time <- params$originTime
  while (cur.time >= 0) {
    # update total rate
    rates <- c(sum(birth.rates), sum(death.rates), sum(migration.rates), 
               sum(transmission.rates))
    total.rate <- sum(rates)
    if (total.rate == 0) {
      break  # nothing can occur, end simulation
    }
    
    # sample waiting time to next event
    wait <- rexp(1, total.rate)
    current.time <- current.time - wait
    if (current.time <= 0) {
      break  # end simulation
    }
    
    # which event?
    eidx <- sample(1:4, 1, prob=rates/total.rate)
    
    if (eidx == 1) {  
      # BIRTH
      cn <- sample(cnames, 1, prob=birth.rates)
      eval(parse(text=paste(cn, "<-", cn, "+1")), envir=e)

    } else if (eidx == 2) {  
      # DEATH
      cn <- sample(cnames, 1, prob=death.rates)
      eval(parse(text=paste(cn, "<-", cn, "-1")), envir=e)
      
    } else if (eidx == 3) { 
      # MIGRATION
      src <- sample(cnames, 1, prob=apply(migration.rates, 1, sum))
      dest <- sample(cnames, 1, prob=migration.rates[src,])
      
      eval(parse(text=paste(src, "<-", src, "-1")), envir=e)
      eval(parse(text=paste(dest, "<-", dest, "+1")), envir=e)
      
    } else {
      # TRANSMISSION
      src <- sample(cnames, 1, prob=apply(transmission.rates, 1, sum))
      dest <- sample(cnames, 1, prob=transmission.rates[src,])
      
      eval(parse(text=paste(src, "<-", src, "-1")), envir=e)
      eval(parse(text=paste(dest, "<-", dest, "+1")), envir=e)
      
    }
    
    # update rates
    birth.rates <- sapply(mod$get.birth.rates(), function(x) {
      eval(parse(text=x), envir=e)
    })
    death.rates <- sapply(mod$get.death.rates(), function(x) {
      eval(parse(text=x), envir=e)
    })
    migration.rates <- apply(mod$get.migration.rates(), MARGIN=c(1,2), 
                             function(x) { eval(parse(text=x), envir=e) })
    transmission.rates <- apply(mod$get.transmission.rates(), MARGIN=c(1,2),
                                function(x) { eval(parse(text=x), envir=e) })
  }
  
}