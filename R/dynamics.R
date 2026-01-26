
simulate.dynamics <- function(mod) {
  # unpack the Model object
  params <- mod$get.parameters()
  comps <- mod$get.compartments()
  cnames <- names(comps)
  sampling <- mod$get.sampling()
  
  # instantiate model parameters
  e <- new.env()
  for (key in names(params)) {
    eval(parse(text=paste(key, "<-", params[[key]])), envir=e)
  }
  
  # instantiate model variables (compartments)
  for (cn in cnames) {
    init.size <- comps[[cn]]$get.size()
    eval(parse(text=paste(cn, "<-", init.size)), envir=e)
  }

  # generate rate matrix
  k <- length(comps)
  rates <- matrix(0, nrow=k, ncol=k, dimnames=list(cnames, cnames))
  for (i in 1:k) {
    src <- comps[[cnames[i]]]
    for (dest in names(src$get.rates())) {
      j <- which(cnames==dest)
      rates[i,j] <- eval(parse(text=src$get.rates()[[dest]]), envir=e)
    }
  }
  
  # are there other stopping criteria?
  targets <- NULL
  if (sampling$mode == "compartment") {
    targets <- sampling$targets
  }
  
  # main simulation loop
  cur.time <- params$originTime
  while (cur.time >= 0) {
    # determine total rates
    counts <- sapply(colnames(rates), function(cn) 
      eval(parse(text=cn), envir=e))
    pop.rates <- counts %*% rates
    total.rate <- sum(pop.rates)
    
    # sample compartment to increment
    
    
    # update rate matrix
    
  }
  
}