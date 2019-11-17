
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
  run <- Run$new(model=model)
  
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
      eventlog$add.event('transmission',  branching.time, NA, NA, 
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
#' \code{sim.outer.tree} simulates transmission, migration and transition events 
#' and fixes them to the timeline of lineage sampled events.  
#' 
#' A *transmission* is the passage of one or more Lineages from a source 
#' Compartment to an uninfected recipient Compartment that did not carry any 
#' Lineages.  
#' 
#' A *migration* is the passage of one or more Lineages from a source Compartment
#' to a recipient Compartment that has already been infected (it already carries
#' one or more Lineages).
#' 
#' A *transition* occurs when a Compartment switches to another Type.  This 
#' is useful for handling host death events, for example.
#' 
#' @param model: object of class 'MODEL'
#' @return object of class 'EventLogger'
#' 
#' @examples
#' 
#' #' # load susceptible-infected (SI) compartmental model
#' path <- system.file('extdata', 'structSI.yaml', package='twt')
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
  
  # store fixed sampling times of Lineages as specified by model
  eventlog$store.fixed.samplings(model$get.fixed.samplings())
  
  # extract objects from Run object
  comps <- run$get.compartments()
  types <- run$get.types()  # CompartmentType objects
  
  # maps Compartments to CompartmentTypes by index
  indiv.types <- sapply(unlist(comps), function(a){a$get.type()$get.name()})
  
  # record population totals and transmission rates for all Types
  init.conds <- run$get.initial.conds()
  popn.totals <- init.conds$size  # list of CompartmentType->size
  
  # extract transmission, migration and transition rates as convenient named matrices
  popn.rates <- .get.rate.matrices(types)

  # record max sampling times (furthest back in time) of lineages for each Compartment
  init.samplings <- .get.initial.samplings(comps)
  
  
  # generate transmission events based on population dynamics and Compartments' 
  # (time, recipient Type, source Type)
  events <- .sample.outer.events(
    types, init.conds, popn.rates, 
    split(init.samplings$init.sample, init.samplings$type)
    )
  
  
  # generate unsampled hosts
  .promote.unsampled.hosts(run, events)
  
  
  # expand `time.bands` to include unsampled hosts
  unsampled <- run$get.unsampled.hosts()
  time.bands <- c( init.samplings$init.sample, rep(NA, length(unsampled)) )
  names(time.bands) <- c( row.names(init.samplings), run$get.names(unsampled) )
  
  
  # unsampled infected now calculated and generated, set source population
  infected <- c(comps, unsampled)
  
  # assign transmission times
  recipient.types <- names(apply(popn.rates, 2, sum) > 0)
  
  . <- sapply(recipient.types, function(rtype) {
    r.indices <- which( sapply(infected, function(comp) {
      comp$get.type()$get.name() == rtype 
      }) )
    r.comps <- infected[r.indices]
    r.init.samplings <- time.bands[ which(names(time.bands) %in% names(infected)[r.indices]) ]
    r.events <- events[ which(events$r.type == rtype), ]
    
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
  
  eventlog <- .assign.source.compartments(eventlog, comps, source.popns, events)

  return(run)
}




#' .get.rate.matrices
#' 
#' Stores population rates for each CompartmentType specified by the user.
#' Handle unspecified rate parameters.
#
#' @param types:  list of CompartmentType objects as returned by Run$get.types()
#' @return  list of named matrices of transmission and migration rates for each 
#'          CompartmentType to CompartmentType
#'          
#' @keywords internal
.get.rate.matrices <- function(types) {
  n <- length(types)
  dimnames <- list(names(types), names(types))
  
  # transmission, transition and migration
  t.rates <- matrix(nrow=n, ncol=n, dimnames=dimnames)
  s.rates <- matrix(nrow=n, ncol=n, dimnames=dimnames)
  m.rates <- matrix(nrow=n, ncol=n, dimnames=dimnames)
  
  for (source in types) {
    # row = source
    s.name <- source$get.name()
    for (r.name in names(types)) {
      # column = recipient
      t.rate <- source$get.branching.rate(r.name)
      if (is.null(t.rate)) {
        # user did not specify transmission rate
        t.rate <- 0
      }
      if (t.rate < 0) {
        stop("Error in .get.rate.matrices: detected negative transmission rate for ",
             "CompartmentType ", source$get.name(), " to ", r.name)
      }
      
      # default to zero if null (unspecified)
      t.rates[s.name, r.name] <- ifelse(is.null(t.rate), 0, t.rate)
      
      
      s.rate <- source$get.transition.rate(r.name)
      if (is.null(s.rate)) {
        s.rate <- 0
      }
      if (s.rate < 0) {
        stop("Error in .get.rate.matrices: detected negative transition rate for ",
             "CompartmentType ", source$get.name(), " to ", r.name)
      }
      s.rates[s.name, r.name] <- ifelse(is.null(s.rate), 0, s.rate)
      
      
      m.rate <- source$get.migration.rate(r.name)
      if (is.null(m.rate)) {
        m.rate <- 0
      }
      if (m.rate < 0) {
        stop("Error in .get.rate.matrices: detected negative migration rate for ",
             "CompartmentType ", source$get.name(), " to ", r.name)
      }
      m.rates[s.name, r.name] <- ifelse(is.null(m.rate), 0, m.rate)
    }
  }
  
  # check diagonal of transition rate matrix
  if ( any(diag(s.rates) > 0) ) {
    warning("Detected positive transition rates to self, zeroing out.")
    diag(s.rates) <- 0
  }
  
  list(transmission=t.rates, transition=s.rates, migration=m.rates)
}



#' .get.initial.samplings
#' 
#' Support function for sim.outer.tree.  Retrieves first sampling time of a Lineage 
#' for each sampled infected Compartment.
#'
#' @param compartments:  list of sampled infected Compartment objects
#'
#' @return named vector of first Lineage sampling times for each Compartment
#' @keywords internal
.get.initial.samplings <- function(compartments) {
  data.frame(
    type=sapply(compartments, function(comp) comp$get.type()$get.name()),
    init.sample=sapply(compartments, function(comp) {
      max( sapply(comp$get.lineages(), function(l) l$get.sampling.time()) )
    })
  )
}


#' .scale.contact.rates
#' 
#' Helper function for .sample.outer.events - calculate population rates 
#' given per-contact rate (\beta) and population sizes of source and
#' recipient populations, *i.e.*, \beta*S*I
#' 
#' @param mx:  rate matrix for source->recipient events
#' @param size:  vector of subpopulation sizes
#' @return  matrix of same dimensions as 'mx'
#' 
#' @keywords internal
.scale.contact.rates <- function(mx, susceptible, infected) {
  if (nrow(mx) != ncol(mx)) {
    stop("Error in .scale.contact.rates: argument 'mx' must be square matrix, ",
         "dimension is ", dim(mx))
  }
  if (nrow(mx) != length(infected) || 
      nrow(mx) != length(susceptible)) {
    stop("Error in .scale.contact.rates: dimension of matrix must equal length ",
         "of susceptible/infected vectors.")
  }
  
  temp <- t(apply(mx, 1, function(row) row*susceptible))
  apply(temp, 2, function(col) col*infected)
}


#' .sample.outer.events
#' 
#' Support function for sim.outer.tree.  Samples transmission, transition and
#' migration events based on population dynamics of the Model.
#' 
#' @param types:  list of CompartmentTypes as returned by Run$get.types()
#' @param init.conds: initial conditions, comprising (1) time scale of simulation
#'        (2) initial population size per CompartmentType and (3) Type
#'        of index case
#' @param popn.rates:  names matrices for transmission, transition and migration
#'        rates.
#' @param init.samplings: list of first Lineage sampling times per Compartment,
#'        grouped by CompartmentType
#' @param max.attempts: number of tries to sample events
#'        
#' @return data frame of outer tree events, each made up of: time, 
#'         event type, recipient CompartmentType and source CompartmentType
#'         
#' @keywords internal
.sample.outer.events <- function(types, init.conds, popn.rates, init.samplings, 
                                 max.attempts=10) {
  
  # convert data frame into list
  init.by.type <- split(init.samplings$init.sample, init.samplings$type)
  
  attempt <- 1
  while (attempt < max.attempts) {
    
    # initialize populations at origin (named vectors)
    susceptible <- unlist(init.conds$size)
    infected <- sapply(susceptible, function(x) 0)
    susceptible[init.conds$indexType] <- susceptible[init.conds$indexType] - 1
    infected[init.conds$indexType] <- 1
    
    # prepare outcome containers
    events <- data.frame(time=numeric(), event.type=character(), r.type=character(), 
                         s.type=character(), stringsAsFactors = FALSE)
    counts <- c()
    
    # simulate trajectories in forward time (indexed in reverse, ha ha!)
    current.time <- init.conds$originTime
    
    while (current.time >= 0) {
      # scale per-contact rates
      t.rates <- .scale.contact.rates(popn.rates[['transmission']], susceptible, infected)
      m.rates <- .scale.contact.rates(popn.rates[['migration']], infected, infected)
      
      # scale non-contact rates
      s.rates <- popn.rates[['transition']] * (susceptible+infected)
      
      total.rate <- sum(t.rates, s.rates, m.rates)
      if (total.rate == 0) {
        # nothing can happen, simulation over
        break
      }
      
      # draw waiting time to next event
      wait <- rexp(1, rate=total.rate)
      current.time <- current.time - wait
      if (current.time < 0) {
        # end of simulation
        break
      }
      
      # append counts to outcome container BEFORE event
      counts <- rbind(counts, c(susceptible, infected))
      
      # which event?
      if ( runif(1, max=total.rate) <= sum(t.rates) ) {
        event.type <- 'transmission'
        
        if (length(t.rates) == 1) {
          # only one CompartmentType
          source <- types[[1]]
          recipient <- types[[1]]
        } 
        else {
          row <- sample(1:nrow(t.rates), 1, prob=apply(t.rates, 1, sum))
          source <- types[[row]]
          recipient <- sample(types, 1, prob=t.rates[row, ])[[1]]
        }
        
        # update counts (recipient becomes source type)
        susceptible[recipient$get.name()] <- susceptible[recipient$get.name()] - 1
        infected[source$get.name()] <- infected[source$get.name()] + 1
      }
      else {
        total.s.rate <- sum(s.rates)
        total.m.rate <- sum(m.rates)
        
        if (runif(1, max=total.s.rate+total.m.rate) < total.s.rate) {
          event.type <- 'transition'
          
          if (length(s.rates) == 1) {
            source <- types[[1]]
            recipient <- types[[1]]
          }
          else {
            row <- sample(1:nrow(s.rates), 1, prob=apply(s.rates, 1, sum))
            source <- types[[row]]
            recipient <- sample(types, 1, prob=s.rates[row, ])[[1]]            
          }
          
          # susceptible or infected?
          sus <- susceptible[source$get.name()]
          inf <- infected[source$get.name()]
          if (runif(1, max=sus+inf) < sus) {
            susceptible[source$get.name()] <- sus - 1
            susceptible[recipient$get.name()] <- susceptible[recipient$get.name()] + 1
          } 
          else {
            infected[source$get.name()] <- inf - 1
            infected[recipient$get.name()] <- infected[recipient$get.name()] + 1
          }
        }
        else {
          event.type <- 'migration'
          
          if (length(s.rates) == 1) {
            source <- types[[1]]
            recipient <- types[[1]]
          }
          else {
            row <- sample(1:nrow(m.rates), 1, prob=apply(m.rates, 1, sum))
            source <- types[[row]]
            recipient <- sample(types, 1, prob=m.rates[row, ])[[1]]
          }
          
          # migration does not affect counts
        }
      }
      
      # append event
      events <- rbind(events, data.frame(
        time=current.time, 
        event.type=event.type,
        r.type=recipient$get.name(),
        s.type=source$get.name(),
        stringsAsFactors=FALSE
      ))

    } 
    
    # check that there are enough Compartments of each Type to be sampled
    checks <- sapply(1:length(infected), function(i) {
      ctype <- names(infected)[i]
      count <- infected[ctype]
      samp.times <- init.by.type[[ctype]]
      if (is.null(samp.times)) {
        stop("Error in .sample.outer.events: CompartmentType ", ctype, 
             " not found in init.samplings")
      }
      length(samp.times) <= count
    })
    
    if (all(checks)) {
      break
    }
  }
 
  
  if (attempt == max.attempts && !all(checks)) {
    stop ('Transmission times generated overshoots given `sampling.times` of ', 
          'Compartments. Some options to fix this error include: changing the ', 
          'origin time of the epidemic to be further back in time or increasing ',
          'transmission rates.')
  }
  
  # merge events and counts into a single data frame
  counts <- as.data.frame(counts)
  names(counts) <- c(paste('S.', names(susceptible), sep=''), 
                     paste('I.', names(infected), sep=''))
  
  cbind(events, counts)
}


#' .assign.events
#' 
#' 
.assign.events <- function(run, events) {
  eventlog <- run$get.eventlog()
  eventlog$clear.events()
  
  # start with sampled infected Compartments
  active <- run$get.compartments()
  
  # iterate through events in reverse time (start with most recent)
  for (row in seq(nrow(events), 1, -1)) {
    
    if (length(active) == 1) {
      # reached the root of sampled Compartments
      break
    }
    
    types <- sapply(active, function(comp) comp$get.type()$get.name())
    
    e <- events[row, ]  # go to the next event
    
    # TODO: filter by initial sampling time for transmission events
    #       Need to make sure Compartment is not infected after a Lineage 
    #       was sampled.
    n.active.recipients <- sum(types == e$r.type)
    n.active.sources <- sum(types == e$s.type)
    
    # recipients who were uninfected (susceptible) BEFORE event
    n.recipients <- as.numeric(e[paste('S.', e$r.type, sep='')])
    n.sources <- as.numeric(e[paste('I.', e$s.type, sep='')])
    
    
    if (e$event.type == 'transmission' || e$event.type == 'migration') {
      
      if (runif(1, max=n.recipients) < n.active.recipients) {
        # recipient is an active Compartment
        recipient <- sample(active[types==e$r.type], 1)[[1]]
        
        if (e$event.type == 'transmission') {
          # remove recipient from active list
          active[recipient$get.name()] <- NULL
          types <- sapply(active, function(comp) comp$get.type()$get.name())
          
          # update counts (recipient is no longer infected)
          if (e$s.type == e$r.type) {
            n.active.sources <- n.active.sources - 1
            n.sources <- n.sources - 1
          }
        }
        # migration can *potentially* remove compartment from active list
        # but this requires inner tree simulation (at Lineage level)
        
        if (runif(1, max=n.sources) < n.active.sources) {
          # source is an active Compartment
          source <- recipient
          while (source$get.name() == recipient$get.name()) {
            # make sure source is different Compartment
            source <- sample(active[types==e$s.type], 1)[[1]] 
          }
        }
        else {
          # source is an unsampled Compartment 
          source <- run$generate.unsampled(e$s.type)
          active[[source$get.name()]] <- source
        }
        
        # update event log
        eventlog$add.event(
          time = e$time,
          type = e$event.type,
          comp1 = recipient$get.name(),
          comp2 = source$get.name()
        )
      }
      # otherwise recipient is not a sampled Compartment - ignore
    }
    
    else if (e$event.type == 'transition') {
      
      if (runif(1, max=n.recipients+n.sources) < n.active.recipients) {
        compartment <- sample(active[types==e$r.type], 1)[[1]]
        
        old.type <- compartment$get.type()$get.name()
        compartment$set.type(run$get.types()[[e$s.type]])
        
        eventlog$add.event(
          time = e$time,
          type = 'transition',
          comp1 = compartment$get.name(),
          type1 = old.type,
          type2 = compartment$get.type()$get.name()
        )
      }
      # otherwise ignore transition of unsampled Compartment
    }
    
    else {
      stop("Error in .assign.events: unrecognized event type ", 
           e$event.type)
    }
  }
}



#' .promote.unsampled.hosts
#' Transmission and migration events 
.promote.unsampled.hosts <- function(run, events) {
  last.indiv <- T  # FIXME: what is the role of this flag?
  for (t in types) {
    # difference between number of recipient Compartments and 
    #  number of sampled Compartments of this Type
    num.unsampled <- length(which(events$r.type == t$get.name())) - 
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

