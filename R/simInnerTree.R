#' sim.inner.tree
#' 
#' Simulate the coalescence of Pathogen lineages within Hosts, and resolve 
#' superinfection events that may involve sampled lineages.
#' 
#' @param outer:  R6 object of class `OuterTree`
#' @param mod:  R6 object of class `Model`
#' 
#' @return R6 object of class `InnerTree`
sim.inner.tree <- function(outer, mod) {
  inner <- InnerTree$new(outer, mod)
  active <- inner$get.active()
  inactive <- inner$get.inactive()
  
  env <- new.env()  # to instantiate bottleneck/coalescence expressions
  
  # iterate through events in outer tree in reverse time (most recent first)
  events <- outer$get.log()
  events$time <- as.numeric(events$time)
  events <- events[order(events$time, decreasing=TRUE), ]
  
  row <- 1
  time.delta <- -diff(events$time)  # time to next event
  p.index <- 1  # uniquely label Pathogens
  
  while (row <= nrow(events)) {
    e <- events[row, ]  # retrieve this event
    e.prev <- ifelse(row < nrow(events), events[row+1, ], NA)  # FIXME
    
    # check if coalescence occurs before next event
    if (inner$n.active() > 0) {
      wait.time <- .rcoal(active, mod, env)
      if (!is.na(wait.time) & wait.time < time.delta[row]) {
        .do.coalescent(e, inner)
        time.delta[row] <- time.delta[row] - wait.time
      }
    }
    
    if (e$event == "migration" & inner$has.target(e$to.comp)) {
      # this is a sampling event - any other migration is ignored
      .do.sampling(e, inner)
      row <- row + 1
      next  # go to next event
    } 
    
    
  }
}


#' Draw exponential waiting times to coalescence in all active Hosts and 
#' return the shortest time.
#' @param active:  R6 object of class `HostSet`, containing Hosts carrying 
#'                 Pathogens that may coalesce
#' @param mod:  R6 object of class `Model`
#' @return list containing:
#'         - dt: numeric, minimum waiting time to coalescence
#'         - host name
#' @keywords internal
#' @noRd
.rcoal <- function(active, mod, envir=baseenv()) {
  wait.times <- sapply(active$get.hosts(), function(h) {
    k <- h$count.pathogens()
    if (k < 2) { return(NA) } else {
      # only Hosts with multiple Pathogen lineages can have coalescent events
      comp <- h$get.compartment()
      expr <- mod$get.coalescent.rate(comp)
      
      # TODO: compute time since infection for time-varying rates
      rate <- eval(parse(text=expr), envir=envir)
      if (rate == 0) { return(NA) } else {
        return(rexp(1, choose(k, 2)*rate))
      }
    }
  })
  wait.times <- setNames(wait.times, active$get.names())
  
  if ( all(is.na(wait.times)) ) { return(NA) } else {
    return(list(
      dt=min(wait.times, na.rm=TRUE),
      host=names(wait.times)[which.min(wait.times)]
      ))
  }
}


#' @param e:  row from event log
#' @param inner:  R6 object of class `InnerTree`
#' @keywords internal
#' @noRd
.do.sampling <- function(e, inner) {
  # transfer Host from sampled to active HostSets
  sampled <- inner$get.sampled()
  host <- sampled$get.host.by.name(e$from.host, remove=TRUE)
  active <- inner$get.active()
  active$add.host(host)
  
  # create a Pathogen and link to sampled Host
  path <- inner$new.pathogen(e$time)
  host$add.pathogen(path)
  
  event <- list(
    time=e$time, event='sampling', from.host=host$get.name(), to.host=NA,
    pathogen1=path$get.name(), pathogen2=NA
  )
  inner$add.event(event)
}


#' .do.infection
#' Transfer Pathogen from 
#' @param e:  row from outer event log
#' @param inner:  R6 object of class `InnerTree`
#' @keywords internal
#' @noRd
.do.infection <- function(e, inner) {
  active <- inner$get.active()
  
  # remove Pathogen from recipient Host
  recipient <- active$get.host.by.name(e$to.host)
  path <- recipient$sample.pathogen(remove=TRUE)
  if (recipient$count.pathogens() == 0) {
    active$remove.host(recipient)  # deactivate if Host is empty
  }
  
  # transfer Pathogen to source Host
  source <- active$get.host.by.name(e$from.host)
  source$add.pathogen(path)
}

#' @param e:  row from outer event log
#' @param inner:  R6 object of class `InnerTree`
#' @keywords internal
#' @noRd
.do.coalescent <- function(e, inner) {
  
}