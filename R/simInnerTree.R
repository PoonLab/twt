#' sim.inner.tree
#' 
#' Simulate the coalescence of Pathogen lineages within Hosts, and resolve 
#' superinfection events that may involve sampled lineages.
#' 
#' @param outer:  R6 object of class `OuterTree`
#' 
#' @return R6 object of class `InnerTree`
#' @export
sim.inner.tree <- function(outer) {
  inner <- InnerTree$new(outer)
  mod <- inner$get.model()
  active <- inner$get.active()
  inactive <- inner$get.inactive()
  
  env <- new.env()  # to instantiate bottleneck/coalescence expressions
  
  # iterate through events in outer tree in reverse time (most recent first)
  events <- outer$get.log()
  events$time <- as.numeric(events$time)
  events <- events[order(events$time, decreasing=TRUE), ]
  
  time.delta <- -diff(events$time)  # time to next event
  row <- 1
  p.index <- 1  # uniquely label Pathogens
  
  while (row <= nrow(events)) {
    e <- events[row, ]
    t.delta <-  ifelse(row==1, NA, time.delta[row-1])
    
    # check if coalescence occurs before next event
    if (row > 1 & inner$n.active() > 0) {
      wait.time <- .rcoal(active, mod, env)
      if ( !all(is.na(wait.time)) ) {
        if ( wait.time$dt < t.delta ) {
          t.delta <- t.delta - wait.time$dt 
          .do.coalescent(wait.time$host, inner, e$time + t.delta)
          next
        }
      } 
    }
    
    # otherwise no coalescence, handle the event
    if (e$event == "migration") {
      if ( inner$has.target(e$to.comp) ) {
        # this is a sampling event, add a new Pathogen
        .do.sampling(e, inner)
      } else {
        # record the migration event for any Pathogens carried by the Host
        .migrate.pathogens(e, inner)
      } 
    } else if (e$event == "transmission") {
      .do.infection(e, inner)
    }
    
    row <- row + 1  # go to next event
  }
  
  count <- active$count.type()
  if (count > 1) {
    warning("Multiple (", count, ") hosts remain in active HostSet at end of ",
            "simulation.")
  } else if (count == 0) {
    warning("Empty active HostSet at end of simulation!")
    
  } else {
    # finish coalescence in last host
    root <- active$get.hosts()[[1]]
    while (root$count.pathogens() > 1) {
      wait.time <- .rcoal(active, mod, envir=env)
      .do.coalescent(root$get.name(), inner, e$time - wait.time$dt)
    }
  }
  
  return(inner)
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
    time=e$time, event='sampling', 
    from.comp=e$from.comp, to.comp=e$to.comp,
    from.host=host$get.name(), to.host=NA,
    pathogen1=path$get.name(), pathogen2=NA
  )
  inner$add.event(event)
}


#' .do.infection
#' 
#' Transfer Pathogen from one Host to another.  If this is the first infection 
#' of the recipient Host, and it carries multiple Pathogens, then do a 
#' bottleneck.  If it is a subsequent infection (superinfection), one of the 
#' Pathogens is transferred with probability equal to the coalescence rate.
#' 
#' We assume that superinfection samples [b.size] Pathogens from 
#' a population of [p.size] at random without replacement.  [count] of 
#' these are actively tracked lineages, so 0, 1, ..., [bsize] of them 
#' may be transferred to the super-infecting source Host.  It is possible 
#' that all active lineages are transferred at this event, leaving none 
#' for the first transmission to this Host, in which case the Host is 
#' deactivated.
#' 
#' @param e:  row from outer event log
#' @param inner:  R6 object of class `InnerTree`
#' @return R6 object of class `Host` if it is empty (no more Pathogens)
#' 
#' @keywords internal
#' @noRd
.do.infection <- function(e, inner, envir=baseenv()) {
  active <- inner$get.active()
  inactive <- inner$get.inactive()
  
  # we are going back in time, so a recipient must already be active
  recipient <- active$get.host.by.name(e$to.host)
  
  # source may or may not be an active Host
  source.is.active <- TRUE
  suppressWarnings({
    source <- active$get.host.by.name(e$from.host)
    if (is.null(source)) {
      source <- inactive$get.host.by.name(e$from.host, remove=TRUE)
      source.is.active <- FALSE
    }
  })
  
  # prepare log entry
  event <- list(
    time=e$time, event='transmission', 
    from.comp=e$src.comp, to.comp=e$to.comp,
    from.host=e$from.host, to.host=e$to.host,
    pathogen1=NA, pathogen2=NA
  )
  
  is.infected <- inner$get.model()$get.infected()
  
  if (is.infected[[e$from.comp]]) {
    # superinfection
    count <- recipient$count.pathogens()
    
    # determine bottleneck and population sizes
    expr <- inner$get.model()$get.bottleneck.size(e$to.comp)
    b.size <- eval(parse(text=expr), envir=envir)
    expr <- inner$get.model()$get.pop.size(e$to.comp)
    p.size <- eval(parse(text=expr), envir=envir)
    
    n.transfer <- rhyper(1, count, p.size-count, b.size)
    if (n.transfer > 0) {
      for (i in 1:n.transfer) {
        path <- recipient$remove.pathogen(1)
        source$add.pathogen(path)
        event$pathogen1 <- path$get.name()
        inner$add.event(event)
      }
      if (!source.is.active) { active$add.host(source) }
      
    } else {
      # no Pathogens transferred to source, restore to inactive HostSet
      inactive$add.host(source)
    }
    
  } else {
    # first infection
    
    if (recipient$count.pathogens() > 1) {
      # bottleneck remaining Pathogen lineages
      .do.coalescent(recipient$get.name(), inner, e$time, bottleneck=TRUE)
    }
    
    # move remaining Pathogens from recipient to source
    count <- recipient$count.pathogens()
    for (i in 1:count) {
      path <- recipient$remove.pathogen(1)
      source$add.pathogen(path)
      event$pathogen1 <- path$get.name()  # copy-on-modify
      inner$add.event(event)
    }
    
    if (!source.is.active) { active$add.host(source) }
  }
  
  # deactivate if Host is empty
  if (recipient$count.pathogens() == 0) {
    host <- active$remove.host(recipient)
    inactive$add.host(host)
  }
}


#' Propagate an outer tree migration event to every Pathogen carried by 
#' the affected Host, updating the inner event log.
#' 
#' @param e:  row from outer event log
#' @param inner:  R6 object of class `InnerTree`
#' @keywords internal
#' @noRd
.migrate.pathogens <- function(e, inner) {
  active <- inner$get.active()
  host <- active$get.host.by.name(e$from.host)
  for (path in host$get.pathogens()) {
    event <- list(
      time=e$time, event='migration', 
      from.comp=e$from.comp, to.comp=e$to.comp,
      from.host=host$get.name(), to.host=NA,
      pathogen1=path$get.name(), pathogen2=NA
    )
    inner$add.event(event)
  }
}


#' .do.coalescent
#' 
#' A coalescent event is the merging of two Pathogen lineages to their 
#' common ancestor.  To facilitate subsequent steps, both Pathogens are 
#' replaced by a new Pathogen representing their ancestor.  This function is 
#' also used to handle transmission bottlenecks, which are associated with 
#' the first event for the recipient Host.
#' 
#' @param host.name:  character, name of Host for coalescence
#' @param inner:  R6 object of class `InnerTree`
#' @param time:  numeric, time of event in simulation time, to update 
#'               Pathogen start.time or end.time
#' @param bottleneck:  bool, if TRUE then determine bottleneck size and 
#'                     coalesce all Pathogens to that size
#' @param envir:  environment, for evaluating bottleneck size expression
#'                (defaults to base env)
#' 
#' @keywords internal
#' @noRd
.do.coalescent <- function(host.name, inner, time, bottleneck=FALSE,
                           envir=baseenv()) {
  active <- inner$get.active()
  host <- active$get.host.by.name(host.name)
  comp <- host$get.compartment()
  
  count <- host$count.pathogens()
  if (count < 2) {
    warning("Cannot coalesce fewer than two Pathogens in Host ", 
            host$get.name())
    return(NULL)
  }
  
  size <- count - 1  # default to single coalescent event
  if (bottleneck) {
    # determine bottleneck size
    expr <- inner$get.model()$get.bottleneck.size(comp)
    size <- eval(parse(text=expr), envir=envir)
  }
  
  while (count > size) {
    p1 <- host$sample.pathogen(remove=TRUE)
    p1$set.start.time(time)
    
    p2 <- host$sample.pathogen(remove=TRUE)
    p2$set.start.time(time)
    
    anc <- inner$new.pathogen(time)  # sets end.time
    host$add.pathogen(anc)
    
    # assign ancestral/descendant relations
    anc$add.child(p1)
    anc$add.child(p2)
    p1$set.parent(anc)
    p2$set.parent(anc)
    
    # update inner log
    event <- list(
      time=time, event='coalescent', from.comp=comp, to.comp=NA,
      from.host=host$get.name(), to.host=NA,
      pathogen1=anc$get.name(), pathogen2=p1$get.name()
    )
    inner$add.event(event)
    event$pathogen2 <- p2$get.name()
    inner$add.event(event)
    
    count <- host$count.pathogens()  # update count
  }
}