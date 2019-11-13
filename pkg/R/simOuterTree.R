
#' eventlog.from.tree
#' 
#' \code{eventlog.from.tree} converts a Newick tree string or ape::phylo object
#' into a sequence of transmission events that are stored as an eventlog object
#' 
#' @param tree: either a Newick tree string or an object of class 'phylo' (ape)
#' that represents the transmission tree (history).  Internal node labels must 
#' be present to indicate transmission sources.
#' 
#' @return
#'   An object of class EventLogger initialized with a list of fixed transmission
#'   events that were extracted from 'tree'.
#' 
#' @examples
#' e <- eventlog.from.tree('(((A:1,B:1)B:1,C:1)C:1,D:1)D:1;')
#' e  # print contents
#' 
#' 
#' @seealso sim.outer.tree, init.branching.events
#' @export
eventlog.from.tree <- function(tree) {
  # intiialize return value
  e <- EventLogger$new()
  
  # check input
  if (is.character(tree)) {
    phy <- read.tree(text=tree)
  } else if (any(class(tree) == 'phylo')) {
    phy <- tree
  } else {
    stop("Error: arg 'tree' must be Newick tree string or 'phylo' object.")
  }
  
  # tree should be rooted and binary
  if (!is.rooted(phy)) {
    stop("Error in eventlog.from.tree(), tree must be rooted")
  }
  if (!is.binary(phy)) {
    stop("Error in eventlog.from.tree(), tree must be binary")
  }
  if (is.null(phy$node.label)) {
    stop('Node labels must be present in host transmission tree to indicate sources.')
  }
  phy$labels <- c(phy$tip.label, phy$node.label)
  phy$depths <- node.depth.edgelength(phy)
  phy$heights <- max(phy$depths) - phy$depths
  
  # iterate over edges in tree
  . <- sapply(1:nrow(phy$edge), function(x) {
    s_ind <- phy$edge[x, 1]  # source
    r_ind <- phy$edge[x, 2]  # recipient
    
    # nodes are numbered terminal first (1,2,...,n) and internal nodes
    # are numbered root first.  Sources are always internal
    source_label <- phy$labels[s_ind]
    recipient_label <- phy$labels[r_ind]

    if (source_label != recipient_label) {
      branching_time <- phy$edge.length[x]
      e$add.event('transmission', phy$heights[s_ind], NA, NA, 
                  recipient_label, source_label)
    }  # otherwise no transmission on branch
  })
  
  e
}


#' init.branching.events
#' 
#' `init.branching.events`` initializes an EventLogger object with 
#' a list of fixed transmission events that are specified by the user 
#' in a YAML file that has been parsed into a MODEL object.
#' 
#' @param model: an object of class MODEL from which we will extract the 
#' transmission tree information.
#' @param eventlog: an object of class EventLog to update in-place
#' 
#' @examples
#' # load model
#' path <- system.file('extdata', 'chain.yaml', package='twt')
#' settings <- yaml.load_file(path)
#' mod <- Model$new(settings)
#' 
#' run <- init.branching.events(mod)
#' plot(run$get.eventlog())
#' 
#' @seealso eventlog.from.tree, sim.outer.tree
#' @export
init.branching.events <- function(model, eventlog=NA) {
  # initialize Run object from Model
  run <- Run$new(model)
  
  if (is.environment(eventlog)) {
    # bind EventLogger to Run object
    run$set.eventlog(eventlog)
  } else {
    eventlog <- run$get.eventlog()
  }
  
  # store fixed sampling times of the tips for plotting functions
  eventlog$store.fixed.samplings(run$get.fixed.samplings())
  
  # if the user input includes a tree (host tree) then add transmission events
  comps <- run$get.compartments()
  lineages <- run$get.lineages()
  locations <- sapply(lineages, function(l) l$get.location()$get.name())

  . <- sapply(comps, function(comp) {
    branching.time <- comp$get.branching.time()
    
    if (is.numeric(branching.time)) {
      
      if ( is.element('R6', class(comp$get.source())) ) {
        source <- comp$get.source()$get.name()
      } 
      else {
        # FIXME: why is this necessary?
        source <- comp$get.source()
      }
      
      # add transmission event to EventLogger object
      eventlog$add.event('transmission', branching.time, NA, NA, 
                         comp$get.name(), source)
    }
    else {
      
      if (branching.time == 'NA' && comp$get.source() == 'NA') {
        # the index case will have no branching time specified
      } 
      else {
        stop ("Error in init.branching.events(): cannot parse Compartment ",
              "with branching.time ", branching.time, 
              " and source ", comp$get.source())
      }
    }
    
  })
  
  return(run)
}


#' sim.outer.tree
#' 
#' After the objects are generated from user inputs ('loadInputs.R'), 
#' we need to initialize the list of fixed events.  We anticipate 
#' three use cases:
#' 1. User provides a string or object representing the transmission tree
#'    that should be converted into an EventLogger object to use for
#'    inner tree simulation (`eventlog.from.tree`).
#' 2. User manually inputs a host transmission tree into YAML format under 
#'    Compartments header (`init.branching.events`).
#' 3. No host tree provided, transmission events need to be simulated under
#'    user-specified model (`sim.outer.tree`).
#'    
#' \code{sim.outer.tree} simulates transmission events and fixes them to the 
#' timeline of lineage sampled events.
#' 
#' @param model: object of class 'MODEL'
#' @return object of class 'EventLogger'
#' 
#' @examples
#' 
#' #' # load susceptible-infected (SI) compartmental model
#' path <- system.file('extdata', 'SI.yaml', package='twt')
#' settings <- yaml.load_file(path)
#' model <- Model$new(settings)
#' 
#' run <- sim.outer.tree(model)
#' run$get.eventlog()
#' 
#' @seealso init.branching.events, eventlog.from.tree
#' @export
sim.outer.tree <- function(model) {
  run <- Run$new(model)
  
  # initialize a new object as return value
  eventlog <- run$get.eventlog()
  
  # ptm <- proc.time()   # benchmark start time
  
  # store fixed sampling times of Lineages as specified by model
  eventlog$store.fixed.samplings(model$get.fixed.samplings())
  
  # extract objects from Run object
  comps <- run$get.compartments()
  types <- run$get.types()  # CompartmentType objects
  
  # maps Compartments to CompartmentTypes by index
  indiv.types <- sapply(unlist(comps), function(a){a$get.type()$get.name()})
  
  
  # record population totals and transmission rates for all Types
  init.conds <- run$get.initial.conds()
  popn.totals <- init.conds$size
  popn.rates <- .calc.popn.rates(types)  # matrix of transmission rates
  
  
  # record max sampling times (furthest back in time) of lineages for each Compartment
  time.bands <- .store.initial.samplings(comps, model$get.lineages(), popn.rates)
  
  
  # generate transmission events based on population dynamics and Compartments' 
  # (time, recipient Type, source Type)
  t_events <- .calc.transmission.events(init.conds, popn.rates, time.bands)
  
  
  # generate unsampled hosts
  last.indiv <- T  # FIXME: what is the role of this flag?
  for (t in types) {
    # difference between number of recipient Compartments and 
    #  number of sampled Compartments of this Type
    num.unsampled <- length(which(t_events$r_type == t$get.name())) - 
      length(which(indiv.types == t$get.name()))
    
    if (num.unsampled > 0) {
      if (last.indiv) {
        # this is the reason why unsampled are not generated at the 
        # initialization of the MODEL object
        # FIXME: why do we have one extra unsampled host of the first Type?
        run$generate.unsampled(num.unsampled+1, t)      
        last.indiv <- F
      } else {
        run$generate.unsampled(num.unsampled, t)
      }
    }
  }
  
  
  # expand `time.bands` to include unsampled hosts
  time.bands <- c(time.bands, rep(NA, length(run$get.unsampled.hosts())))
  names(time.bands) <- c(
    names(time.bands)[nzchar(x=names(time.bands))], 
    model$get.names(run$get.unsampled.hosts())
    )
  
  
  # unsampled infected now calculated and generated, set source population
  infected <- c(comps, run$get.unsampled.hosts())
  
  # assign transmission times
  recipient.types <- names(apply(popn.rates, 2, sum) > 0)
  
  . <- sapply(recipient.types, function(rtype) {
    r.indices <- which( sapply(infected, function(comp) {
      comp$get.type()$get.name() == rtype 
      }) )
    r.comps <- infected[r.indices]
    r.init.samplings <- time.bands[ which(names(time.bands) %in% names(infected)[r.indices]) ]
    r.events <- t_events[ which(t_events$r_type == rtype), ]
    
    if (nrow(r.events) > 0) {
      type.comp <- types[[ which(sapply(types, function(t) t$get.name() == rtype)) ]]
      assigned.times <- .assign.transmission.times(
        r.comps, r.events, r.init.samplings, type.comp
        )
    }
  })
  
  # After transmission times are matched with infected Compartments as recipients, 
  # now have to assign source Compartments
  # Order infection times from most recent to furthest back in time
  btimes <- sapply(comps, function(x) x$get.branching.time())
  comps[sapply(btimes, is.null)] <- NULL
  new.order <- order(sapply(comps, function(x) x$get.branching.time()))
  comps <- comps[ new.order ]
  
  
  # each assigned transmission time is associated with an event, which determines 
  # what TYPE its source is, just not which Compartment in particular
  source.popns = 
    source.popns.names <- setNames(vector(length(types), mode="list"), 
                                   names(types))

  # separate source populations into lists that are as many as the number of distinct 
  # Types in the model
  for (x in 1:length(names(source.popns))) {
    specific.type <- names(source.popns)[x]
    
    # get all source Compartments of this Type
    s.pop.by.type <- sapply(infected, function(y) {
      if (y$get.type()$get.name() == specific.type) {
        
        # include unsampled Compartments as potential sources
        # root no branching time
        if (length(y$get.branching.time()) == 0 
            || is.null(y$get.branching.time())) {y}
        
        # exclude Compartments sampled *after* the recipient
        else if (y$get.branching.time() > comps[[1]]$get.branching.time()) {y}
      } 
    })
    
    s.pop.by.type[sapply(s.pop.by.type, is.null)] <- NULL
    source.popns[[specific.type]] <- s.pop.by.type
  }
  
  eventlog <- .assign.source.compartments(eventlog, comps, source.popns, t_events)

  return(run)
}



#' .assign.source.compartments
#' 
#' Assign source Compartments for all sampled and promoted unsampled infected
#' Compartments.
#' 
#' @keywords internal
.assign.source.compartments <- function(eventlog, comps, source.popns, t_events) {
  
  numActive <- length(comps)
  
  while (numActive >= 1) {
    # first recipient == most recent infection time == largest list of 
    # `sources` --> efficiency
    r_ind_comps <- 1                            
    recipient <- comps[[r_ind_comps]]
    r_name <- recipient$get.name()
    r_type <- recipient$get.type()$get.name()
    
    # remove chosen recipient from relevant lists (`comps` and possibly `source.popns`)
    comps[[ r_ind_comps ]] <- NULL
    ind_source_popn <- which(names(source.popns[[r_type]]) == r_name)
    if (length(ind_source_popn) != 0) {
      # TRUE for all iterations except initial iteration case (bc starting 
      # recipient already removed)
      source.popns[[r_type]][[ind_source_popn]] <- NULL
    }
    
    # list of possible sources is based off of the `s_type` recorded in the 
    # event associated with the transmission time in master copy `t_events`
    if (length(recipient$get.branching.time()) == 0 
        || is.null(recipient$get.branching.time())) {
      # root case, "final" resolved transmission event (aka first recorded 
      # transmission, furthest back in time)
      break
      
    } else {
      
      s_type <- t_events[ which(t_events$time == recipient$get.branching.time()), 's_type']
      
      # refine list of sources previously separated by Type also by earlier branching times than recipient's branching time
      refined.list.sources <- sapply(source.popns[[s_type]], function(s) {
        if (length(s$get.branching.time()) == 0 || is.null(s$get.branching.time())) s
        else if (s$get.branching.time() > recipient$get.branching.time()) s
        else NULL
      })
      
      if (length(refined.list.sources) == 0) {
        print(numActive)
        print(source.popns)
        browser()
      }
      refined.list.sources[sapply(refined.list.sources, is.null)] <- NULL   # cleanup
      
      # select source from a list of sources previously separated by Type
      s_ind_s_popn <- sample.int(length(refined.list.sources), 1)
      source <- refined.list.sources[[s_ind_s_popn]]
      s_name <- source$get.name()
      
      eventlog$add.event('transmission', recipient$get.branching.time(), NA, NA, 
                         r_name, s_name)
      #print(eventlog)  ### DEBUGGING ###
    }
    
    # if source is an unsampled infected Compartment, now holds a sampled lineage 
    # we care about (promote us_comp)
    if (source$is.unsampled() 
        && !is.element(source$get.name(), names(comps))
        && length(source$get.branching.time()) != 0) {
      # add us_comp to list `comps` (once first a source, can now be a recipient)
      # AND ONLY if source is not "final" transmission event (b/c it has no source, 
      # therefore can't be recipient)
      comps[[s_name]] <- source
    }
    
    # update recipient object `source` attr
    recipient$set.source(source)
    
    # update number of active compartments (excludes recipients, includes promoted us_comps)
    new.order <- order(sapply(comps, function(x) x$get.branching.time()))
    comps <- comps[ new.order ]
    numActive <- length(comps)
  }
  
  return(eventlog)
}



#' .calc.popn.rates
#' 
#' Stores population rates for each CompartmentType specified by the user.
#
#' @param types:  list of CompartmentType objects
#' @return matrix of population transmission rates for each type to type comparison
#' @keywords internal
.calc.popn.rates <- function(types) {
  popn.rates <- matrix(nrow=length(types),  # source types
                       ncol=length(types),  # recipient types
                       dimnames=list(names(types), names(types)))
  
  for (source in types) {                                                      
    for (recipient in names(types)) {
      # store instrinsic transmission rates for all x->y pairs
      rate <- source$get.branching.rate(recipient)
      popn.rates[source$get.name(), recipient] <- rate
    }
  }
  
  popn.rates
}



#' .store.initial.samplings
#' 
#' Support function for sim.outer.tree.  Stores first sampling time of a Lineage 
#' for each sampled infected Compartment.
#'
#' @param infected:  list of sampled infected Compartment objects
#' @param lineages:  list of Lineage objects
#' @param popn.rates:  named matrix of transmission rates (row = source)
#'
#' @return named vector of first Lineage sampling times for each Compartment
#' @keywords internal
.store.initial.samplings <- function(infected, lineages, popn.rates) {
  # Compartment name for each Lineage
  lineage.locations <- sapply(lineages, function(x) {x$get.location()$get.name()})
  
  # iterate over every sampled infected compartment
  sapply(infected, function(comp) {
    my.type <- comp$get.type()$get.name()
    
    # all rates *to* this recipient's CompartmentType
    my.rates <- popn.rates[,which(colnames(popn.rates) == my.type)]
    if (all(my.rates==0)) {
      # not possible for this compartment to be a recipient
      # sampling Lineage must be from migration event
      NULL
    }
    else{
      my.lineages <- lineages[which(lineage.locations == comp$get.name())]
      sampling.times <- sapply(my.lineages, function(x) x$get.sampling.time())
      max(sampling.times)  # furthest back in time
    }
  })
}



#' .calc.transmission.events
#' 
#' Support function for sim.outer.tree.  Generates transmission events only, 
#' based on population dynamics of the Model.
#' 
#' @param init.conds: initial conditions, comprising (1) time scale of simulation
#'                    (2) initial population size per CompartmentType and (3) Type
#'                    of index case
#' @param popn.rates: rates of transmission between different CompartmentTypes
#' @param init.samplings: list of first sampling times for each Compartment
#' @param max.attempts: number of tries to sample transmission events
#'        
#' @return data frame of transmission events, each made up of: time, 
#'         recipient Type, source Type, and LineageType (not yet implemented)
#'         
#' @keywords internal
.calc.transmission.events <- function(init.conds, popn.rates, init.samplings, 
                                      max.attempts=10) {
  
  # prepare outcome container
  types <- colnames(popn.rates)
  t_events <- data.frame(time=numeric(), r_type=character(), s_type=character(), 
                         stringsAsFactors = FALSE)
  
  for (attempt in 1:max.attempts) {
    
    # current time starts at user-specified time for each epidemic
    current.time <- as.numeric(init.conds['originTime'])
    
    # prepare count vectors
    susceptible <- unlist(init.conds$size)
    susceptible[init.conds$indexType] <- susceptible[init.conds$indexType] - 1
    
    infected <- lapply(init.conds$size, function(x) 0)
    infected[init.conds$indexType] <- 1
    infected <- unlist(infected)
    
    # check sampling times with epidemic start time
    if (any(init.samplings > current.time)) {
      stop ('Not possible to have Compartment initial sampling time(s) precede the start time of the ', 
            'epidemic. Please set the start time of the epidemic further back in time.')
    }
    
    # note, this is forward-time simulation on a reverse-time scale
    # FIXME: need to implement Compartment death events here
    while (current.time > min(init.samplings) && sum(susceptible) > 0) {

      rates <- t(t(infected * popn.rates) * susceptible)
      delta.t <- rexp(1, sum(rates))  # draw waiting time to any event
      
      current.time <- current.time - delta.t
      # transmission event is more recent than the most recent sampling time
      # this event, and any subsequent transmissions, are irrelevant
      if (current.time <= min(init.samplings)) { break }
      
      
      which.evt <- sample(1:length(rates), size=1, prob=as.vector(rates))
      r <- floor((which.evt-1) / nrow(rates)) + 1  # source Type
      s <- ((which.evt-1) %% nrow(rates)) + 1  # recipient Type
      
      # store time, source and recipient types of transmission
      t_events <- rbind(t_events, 
                        list(time=current.time, r_type=types[r], s_type=types[s]),
                        stringsAsFactors=FALSE)
      
      # update counts
      susceptible[r] <- susceptible[r] - 1
      infected[r] <- infected[r] + 1
    }
    

    
    # check if viable # of transmission times @ each unique user-given Compartment 
    # `sampling.time`
    checks <- sapply(unique(init.samplings), function(x) {
      if (sum(init.samplings==x) > length(t_events$time > x)) {
        # not enough viable transmission times at this `sampling.time`
        # regenerate transmission times up to max.attempts
        FALSE
      } else {
        TRUE
      }
    })
    
    if (attempt == max.attempts && any(checks) == FALSE) {
      stop ('Transmission times generated overshoots given `sampling.times` of ', 
            'Compartments. Some options to fix this error include: changing the', 
            'origin time of the "', v.name, '" epidemic to be further back in ', 
            'time or increasing transmission rates.')
    } else if (any(checks) == FALSE) {
      next
    } else {
      break
    }
  
  }
  
  t_events
}


#' .assign.transmission.times
#' 
#' Support function for sim.outer.tree.  Assignment of transmission times specific 
#' to each type.
#'
#' @param infected:  named list of Compartment objects to be assigned transmission times
#' @param events:  list of possible transmission events for the infected
#' @param initial.samplings:  list of respective initial sampling times for each 
#' Compartment
#' @param type:  CompartmentType object
#' 
#' @return  Returning vector of all transmission times, including times used and unused 
#' in this current assignment of transmission times.
#' @keywords internal
.assign.transmission.times <- function(infected, events, initial.samplings, type) {
  
  i.times <- initial.samplings[!is.na(initial.samplings)]  # separate sampled infected comps from 
  u.times <- initial.samplings[is.na(initial.samplings)]  # unsampled infected comps

  # labelled logical vector tracking used and unused transmission times  
  all.t.times <- events$time
  names(all.t.times) <- vector(length=length(events$time))
  
  while (length(i.times) >= 1) {
    # for sampled infected Compartments, assignments based off of Type-specific events 
    # and Type-specific waiting time distribution
    
    earliest.time <- max(i.times, na.rm=T)  # pick Compartment with earliest sampling time (furthest back in time)
    recipients.names <- names(i.times)[which(i.times==earliest.time)]  # group all infected with the same sampling time together
    
    recipients.inds <- sapply(recipients.names, function(x) which(names(infected) == x))
    recipients <- infected[recipients.inds]
    
    # retrieve transmission times with a recipient type equal to the general type of `recipients`
    # weight each transmission time according to the `wait.time.distr` provided for this type
    possibleTimes <- events$time[ which(events$r_type == type$get.name()) ]
    densities <- sapply(possibleTimes, function(x) eval(parse(text=type$get.wait.time.distr())) )
    total.density <- sum(densities)
    weights <- sapply(densities, function(d) d/total.density)
    
    if (length(possibleTimes) == 1) {
      sampledTimes <- possibleTimes[1]
    } else if (length(possibleTimes) < length(recipients)) {  # for case where one of these recipients will be labeled as the root node (ie. no unsampled hosts)
      sampledTimes <- sample(possibleTimes, size=length(recipients)-1, prob=weights, replace=F)
    } else {
      sampledTimes <- sample(possibleTimes, size=length(recipients), prob=weights, replace=F)
    }
    
    for (x in 1:length(sampledTimes)) {
      # now assign `sampledTimes` to `recipients`, uniformly distributed
      ind <- sample.int(length(recipients), 1)
      t.time <- sampledTimes[x]
      recipients[[ind]]$set.branching.time(t.time)
      
      used.t.time.ind <- which(all.t.times == t.time)
      # this transmission time has been assigned
      names(all.t.times)[used.t.time.ind] <- TRUE
      
      recipients <- recipients[-ind]  # remove from remaining list of recipients
      events <- events[ -which(events$time == t.time), ]  # remove transmission event from events
    }
    
    if (length(recipients) != 0) sapply(recipients, function(x) x$set.branching.time(NA))
    
    i.times <- i.times[ -which(i.times==earliest.time) ] 
  }
  
  while (length(u.times) >= 1) {
    # for unsampled infected Compartments, assignments based off of Type-specific events 
    # and uniform transmission time distribution
    u.ind <- sample.int(length(u.times), 1)
    u.name <- names(u.times)[u.ind]
    u.comp <- infected[[ which(names(infected) == u.name) ]]
    
    u.event <- sample(nrow(events), 1)
    t.time <- events[u.event, 'time']
    u.comp$set.branching.time(t.time)
    
    used.t.time.ind <- which(all.t.times == t.time)  # this transmission time has been assigned, set name in logic vector to TRUE
    names(all.t.times)[used.t.time.ind] <- TRUE                             
    
    u.times <- u.times[-u.ind]  # remove used transmission time for unsampled Compartment
    events <- events[-u.event,]  # remove used transmission event from events
  }
  
  all.t.times  
}
