
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
#' outer <- eventlog.from.tree('(((A:1,B:1)B:1,C:1)C:1,D:1)D:1;')
#' outer
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
#' mod <- MODEL$new(settings)
#' 
#' e <- init.branching.events(mod)
#' plot(e)
#' 
#' @seealso eventlog.from.tree, sim.outer.tree
#' @export
init.branching.events <- function(model, eventlog=NA) {
  if (is.na(eventlog)) {
    eventlog <- EventLogger$new()
  }
  
  # store fixed sampling times of the tips for plotting functions
  eventlog$store.fixed.samplings(model$get.fixed.samplings())
  
  # if the user input includes a tree (host tree) then add transmission events
  comps <- model$get.compartments()
  lineages <- model$get.lineages()
  locations <- sapply(lineages, function(l) l$get.location()$get.name())

  . <- sapply(comps, function(comp) {
    branching.time <- comp$get.branching.time()
    
    if (is.numeric(branching.time)) {
      
      if (is.R6(comp$get.source())) {
        source <- comp$get.source()$get.name()
        # find all lineages that are located in this compartment
        my.lines <- which(locations == comp$get.name())
        
        lineage <- lineages[[ my.lines ]]$get.name()
      } else {
        source <- comp$get.source()
        lineage <- NA
      }
      
      # add transmission event to EventLogger object
      eventlog$add.event('transmission',  branching.time, lineage, NA, 
                         comp$get.name(), source)
    }
    
  })
  
  eventlog
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
#' mod <- MODEL$new(settings)
#' 
#' tree <- sim.outer.tree(mod)
#' 
#' @seealso init.branching.events, eventlog.from.tree
#' @export
sim.outer.tree <- function(model) {
  # initialize a new object as return value
  eventlog <- EventLogger$new()
  
  # ptm <- proc.time()   # benchmark start time
  
  # store fixed sampling times of Lineages as specified by model
  eventlog$store.fixed.samplings(model$get.fixed.samplings())
  
  # extract objects from MODEL
  comps <- model$get.compartments()           
  compnames <- model$get.names(comps)
  types <- model$get.types()  # CompartmentType objects
  
  # maps Compartments to CompartmentTypes by index
  indiv.types <- sapply(unlist(comps), function(a){a$get.type()$get.name()})
  
  
  # record population totals and transmission rates for all Types
  popn.totals <- model$get.origin.times()
  popn.rates <- .calc.popn.rates(types, indiv.types)  # matrix of transmission rates
  
  
  # record max sampling times (furthest back in time) of lineages for each Compartment
  storage <- .store.initial.samplings(comps, types, model$get.lineages())
  time.bands <- storage$initial.times
  possible.source.types <- storage$s.types
  
  
  # generate transmission events based on population dynamics and Compartments' 
  # initial sampling times: (time, recipient Type, source Type, virus Type)
  t_events <- .calc.transmission.events(popn.totals, popn.rates, time.bands, 
                                        possible.source.types)
  
  
  # generate unsampled hosts
  model$clear.unsampled()  # reset MODEL - issue #68
  last.indiv <- T  # FIXME: what is the role of this flag?
  
  for (t in types) {
    num.unsampled <- length(which(t_events$r_type == t$get.name())) - 
      length(which(indiv.types == t$get.name()))
    
    if (num.unsampled > 0) {
      if (last.indiv) {
        # this is the reason why unsampled are not generated at the 
        # initialization of the MODEL object
        # FIXME: why do we have one extra unsampled host of the first Type?
        model$generate.unsampled(num.unsampled+1, t)      
        last.indiv <- F
      } else {
        model$generate.unsampled(num.unsampled, t)
      }
    }
  }
  
  # update `time.bands` to include unsampled hosts
  time.bands <- c(time.bands, rep(NA, length(model$get.unsampled.hosts())))
  names(time.bands) <- c(names(time.bands)[nzchar(x=names(time.bands))], 
                         model$get.names(model$get.unsampled.hosts()))
  
  # unsampled infected now calculated and generated, set source population
  sources <- c(comps, model$get.unsampled.hosts())
  sources.names <- model$get.names(sources)
  
  
  
  # FIXME: I think this sapply() can be removed
  #
  # # names(possible.source.types) represents names of types that can be recipients
  # . <- sapply(names(possible.source.types), function(x) {
  #   # for each Type that can be a recipient, assign transmission times to all of its 
  #   # infected compartments
  #   r.indices <- which( sapply(sources, function(y) y$get.type()$get.name() == x) )
  #   r.comps <- sources[r.indices]
  #   r.init.samplings <- time.bands[ which(names(time.bands) %in% sources.names[r.indices]) ]
  #   r.events <- t_events[ which(t_events$r_type == x), ]
  #   
  #   
  #   # no transmission events for case of ie. blood compartment
  #   if (nrow(r.events) != 0) {            
  #     # assigns transmission times but also returns a binary (logical) vector of which 
  #     # transmission times have been used
  #     # will store this binary vector with associated transmission times into the 
  #     # MODEL (specific to each type)
  #     type.comp <- types[[ which(sapply(types, function(t) t$get.name() == x)) ]]
  #     
  #     # FIXME: this is an unused variable
  #     assigned.times <- .assign.transmission.times(r.comps, r.events, r.init.samplings, 
  #                                                  type.comp)
  #   } 
  #   
  # })
  
  
  
  # mark1 <- proc.time() - ptm            # benchmark 1 (generation of transmission events)
  # cat("GTE:", round(mark1[['elapsed']],5), "s\n")
  # ptm1 <- proc.time()
  
  
  # After transmission times are matched with infected Compartments as recipients, 
  # now have to assign source Compartments
  # Order infection times from most recent to furthest back in time
  comps[sapply(sapply(comps, function(x) x$get.branching.time()), is.null)] <- NULL 
  new.order <- order(sapply(comps, function(x) x$get.branching.time()))
  comps <- comps[ new.order ]
  compnames <- compnames[ new.order ]
  
  # each assigned transmission time is associated with an event, which determines what TYPE its source is, just not which Compartment in particular
  
  source.popns = source.popns.names <- setNames(vector(length(types), mode="list"), 
                                                names(types))

  # separate source populations into lists that are as many as the number of distinct 
  # Types in the model
  for (x in 1:length(names(source.popns))) {
    specific.type <- names(source.popns)[x]
    
    # get all source Compartments of this Type
    s.pop.by.type <- sapply(sources, function(y) {
      if (y$get.type()$get.name() == specific.type) {
        
        # include unsampled Compartments as potential sources
        # root no branching time
        if (length(y$get.branching.time()) == 0 
            || is.null(y$get.branching.time())) {y}
        
        # exclude Compartments sampled *after* the recipient
        else if (y$get.branching.time() > comps[[1]]$get.branching.time()) {y}
      } 
    })
    s.pop.by.type[sapply(s.pop.by.type, is.null)] <- NULL    # cleanup
    s.pop.by.type.names <- sapply(s.pop.by.type, function(z) z$get.name())
    source.popns[[specific.type]] <- s.pop.by.type
    source.popns.names[[specific.type]] <- s.pop.by.type.names
  }
  
  
  # assign source compartments for all sampled infected and promoted unsampled 
  # infected Compartments
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
    compnames <- compnames[-r_ind_comps]
    ind_source_popn <- which(source.popns.names[[r_type]] == r_name)
    if (length(ind_source_popn) != 0) {
      # TRUE for all iterations except initial iteration case (bc starting 
      # recipient already removed)
      source.popns[[r_type]][[ind_source_popn]] <- NULL
      source.popns.names[[r_type]] <- source.popns.names[[r_type]][-ind_source_popn]
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
      
      eventlog$add.event('transmission', recipient$get.branching.time(), NA, NA, r_name, s_name)
      #print(eventlog)  ### DEBUGGING ###
    }
    
    # if source is an unsampled infected Compartment, now holds a sampled lineage 
    # we care about (promote us_comp)
    if (source$is.unsampled() 
        && source$get.name() %in% compnames == FALSE
        && length(source$get.branching.time()) != 0) {
      # add us_comp to list `comps` (once first a source, can now be a recipient)
      # AND ONLY if source is not "final" transmission event (b/c it has no source, 
      # therefore can't be recipient)
      comps[[length(comps)+1]] <- source
      compnames[[length(compnames)+1]] <- s_name
    }
    
    # update recipient object `source` attr
    recipient$set.source(source)
    
    # update number of active compartments (excludes recipients, includes promoted us_comps)
    new.order <- order(sapply(comps, function(x) x$get.branching.time()))
    comps <- comps[ new.order ]
    compnames <- compnames[ new.order ]
    numActive <- length(comps)
  }
  
  
  # mark2 <- proc.time() - ptm1               # benchmark 2 (assignment of sources)
  # cat("AST:", round(mark2[['elapsed']],5), "s\n")
  # total <- proc.time() - ptm
  # cat("Total:", round(total[['elapsed']],5), "s\n")
    
  eventlog
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
  
  for (x in types) {                                                      
    for (y in names(types)) {
      # store instrinsic transmission rates for all typeA -> typeB pairs
      rate <- x$get.branching.rate(y)
      popn.rates[x$get.name(), y] <- rate
    }
  }
  
  popn.rates
}



#' .store.initial.samplings
#' 
#' Stores first sampling time of a Lineage for each sampled infected 
#' Compartment.
#'
#' @param infected = list of sampled infected Compartment objects
#' @param types = list of CompartmentType objects
#' @param lineages = list of Lineage objects
#'
#' @return named list of possible Types of `source` with recipient types as names
#' and list of first Lineage sampling times for each Compartment
#' @keywords internal
.store.initial.samplings <- function(infected, types, lineages) {
  
  # Generate dictionary of different types of source that each recipient Type could 
  #  possibly receive a transmission from
  possibleSourceTypes <- lapply(types, function(source) {
    rates <- source$get.branching.rates()  # per recipient Type
    names(rates[rates>0])
  })

  time.bands <- vector()  # initial (earliest) sampling times for each Compartment
  lineage.locations <- sapply(lineages, function(x) {x$get.location()$get.name()})
  
  for (i in infected) {
    # filtered to remove all source -> recipient pairs with branching rates of 0
    recipientType <- i$get.type()$get.name()
    
    if (recipientType %in% names(possibleSourceTypes)) {
      i.native.lineages <- which(lineage.locations == i$get.name())
      single.i.lineages <- lineages[ i.native.lineages ]
      single.i.sampling.times <- sapply(single.i.lineages, function(b) {
        b$get.sampling.time()
      })
      
      # store initial sampling time for this recipient Compartment
      time.bands <- c(time.bands, max(single.i.sampling.times))
      names(time.bands) <- c(names(time.bands)[nzchar(x=names(time.bands))], i$get.name())
      
    } else {
      # this compartment will never be a recipient (ie. example3.yaml 'blood' compartment)
      next
    }
    
  }
  
  list(s.types=possibleSourceTypes, initial.times=time.bands)
}



#' .calc.transmission.events
#' 
#' Generates transmission events only, based on population dynamics of the
#' MODEL
#' 
#' @param popn.totals: (named matrix) numbers of susceptible and infected 
#'        Lineages at time t=0, specific to each LineageType (rows) and 
#'        CompartmentType (columns)
#' @param popn.rates: rates of transmission between different CompartmentTypes
#' @param init.samplings: list of first sampling times for each Compartment
#' @param possible.sources: list of possible Sources that each recipient Type can 
#'        receive a transmission from 
#' @param max.attempts: number of tries to sample transmission events
#'        
#' @return data frame of transmission events, each made up of: time, 
#'         recipient Type, source Type, and LineageType (not yet implemented)
#' @keywords internal
.calc.transmission.events <- function(popn.totals, popn.rates, init.samplings, 
                                      possible.sources, max.attempts=10) {
  
  # prepare outcome container
  t_events <- data.frame(time=numeric(), r_type=character(), s_type=character(), 
                         v_type=character(), stringsAsFactors = FALSE)
  #maxAttempts <- 10
  r.types <- names(possible.sources)
  
  
  for (attempt in 1:max.attempts) {
    
    for (v in 1:nrow(popn.totals)) {
      # TODO: if more than one LineageType, will need to make a copy of the 
      # population totals (as you are modifing the totals while generating 
      # events for each type)
      
      v.name <- rownames(popn.totals)[v]
      virus <- popn.totals[v,]
      
      # current time starts at user-specified time for each epidemic
      current.time <- as.numeric(virus['start'])

      # FIXME: this should be partitioned by source Type      
      num.infected <- 1  # =I, index case
      
      # check sampling times with epidemic start time
      if (any(init.samplings > current.time)) {
        stop ('Not possible to have Compartment initial sampling time(s) precede the start time of the "', 
              v.name, '" epidemic. Please set the start time of the epidemic further back in time.')
      }
      
      
      # note, this is forward-time simulation on a reverse-time scale
      while (current.time > min(init.samplings) && 
             all(virus[2:length(virus)] >= 1)) {  #&& all(virus != 1)

        r <- sample(r.types, 1)  # sample recipient Type
        s <- sample(possible.sources[[r]], 1)  # sample source Type
        
        # waiting time ~ exp(-beta * S * I)
        rate <- popn.rates[s, r] * (virus[s]) * num.infected
        delta.t <- rexp(n=1, rate=rate)
        current.time <- current.time - delta.t
        
        # transmission event is more recent than the most recent sampling time
        # this event, and any subsequent transmissions, are irrelevant
        if (current.time < min(init.samplings)) { break }
        
        # store time, source and recipient types of transmission, and virus type
        t_events <- rbind(t_events, 
                          list(time=current.time, r_type=r, s_type=s, v_type=v.name),
                          stringsAsFactors=FALSE)
        
        # update counts of total population  (number of susceptibles decrease by 1, number of infected increases by 1)
        virus[s] <- virus[s] - 1            
        num.infected <- num.infected + 1
        
      }
    
    }
    
    # check if viable # of transmission times @ each unique user-given Compartment `sampling.time`
    checks <- sapply(unique(init.samplings), function(x) {
      if (length(init.samplings==x) > length(t_events$time > x)) {
        # not enough viable transmission times at this `sampling.time`
        # regenerate transmission times up to maxAttempts
        FALSE
      } else {
        TRUE
      }
    })
    
    if (attempt == maxAttempts && any(checks) == FALSE) {
      stop ('Transmission times generated overshoots given `sampling.times` of Compartments. Some options to fix this error include: changing the origin time of the "', v.name, '" epidemic to be further back in time or increasing transmission rates.')
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
#' Assignment of transmission times specific to each type.
#'
#' @param infected:  list of Compartment objects to be assigned transmission times
#' @param events:  list of possible transmission events for the infected
#' @param initial.samplings:  list of respective initial sampling times for each 
#' Compartment
#' @param type:  CompartmentType object
#' 
#' @return  Returning vector of all transmission times, including times used and unused 
#' in this current assignment of transmission times.
#' @keywords internal
.assign.transmission.times <- function(infected, events, initial.samplings, type) {
  
  infected.names <- sapply(infected, function(x) x$get.name())
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
    
    recipients.inds <- sapply(recipients.names, function(x) which(infected.names == x))
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
    u.comp <- infected[[ which(infected.names == u.name) ]]
    
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
