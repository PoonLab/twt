
#' attach.tree
#' 
#' \code{attach.tree} converts a Newick tree string or ape::phylo object
#' into a sequence of transmission events that are stored as an eventlog object,
#' bypassing outer tree simulation
#' 
#' @param tree: either a Newick tree string or an object of class 'phylo' (ape)
#' that represents the transmission tree (history).  Internal node labels must 
#' be present to indicate transmission sources.
#' @param model:  an R6 model of class Model.  The Model should provide the parameters
#'        of a CompartmentType and sampled Compartments.  Currently only one 
#'        CompartmentType is supported by this method, until we can devise a 
#'        procedure for the user to assign unsampled Compartments at internal
#'        nodes to different CompartmentTypes.
#' 
#' @return
#'   An object of class Model
#' 
#' @examples
#' e <- attach.tree('(((A:1,B:1)B:1,C:1)C:1,D:1)D:1;')
#' e  # print contents
#' 
#' 
#' @seealso sim.outer.tree, load.outer.tree
#' DEPRECATED
attach.tree <- function(tree, model) {
  types <- model$get.types()
  if (length(types) > 1) {
    stop("Sorry, attach.tree() only supports Models with one CompartmentType for now.")
  }
  comps <- model$get.compartments()
  
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
    stop("Error in attach.tree(), tree must be rooted")
  }
  if (!is.binary(phy)) {
    stop("Error in attach.tree(), tree must be binary")
  }
  if (is.null(phy$node.label)) {
    stop('Node labels must be present in host transmission tree to indicate sources.')
  }
  
  phy$labels <- c(phy$tip.label, phy$node.label)
  if (all(!grepl("_1$", phy$labels))) {
    # add "_1" suffix to compartment names
    phy$labels <- paste0(phy$labels, "_1")
  }
  else {
    stop("Error in attach.tree(): node labels cannot end with \"_1\"; ",
         "this suffix will be automatically added by this function.")
  }
  
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
      #eventlog$add.event('transmission', phy$heights[s_ind], NA, NA, 
      #            recipient_label, source_label)
    }  # otherwise no transmission on branch
  })
  
  model
}


#' load.outer.tree
#' 
#' `load.outer.tree`` initializes a Run object from a Model object
#' in which the user has explicitly listed transmission events; *i.e.*, 
#' bypassing simulation of the outer tree.
#' 
#' @param model: R6 object of class Model
#' @return R6 object of class Run
#' 
#' @examples
#' # load model
#' path <- system.file('extdata', 'chain.yaml', package='twt')
#' settings <- yaml.load_file(path)
#' mod <- Model$new(settings)
#' 
#' run <- load.outer.tree(mod)
#' plot(run)
#' 
#' @seealso attach.tree, sim.outer.tree
#' @export
load.outer.tree <- function(model) {
  
  # initialize Run object from Model
  run <- Run$new(model)
  eventlog <- run$get.eventlog()
  
  
  # store fixed sampling times of the tips for plotting functions
  eventlog$store.fixed.samplings(run$get.fixed.samplings())
  
  # if the user input includes a tree (host tree) then add transmission events
  comps <- run$get.compartments()
  lineages <- run$get.lineages()
  locations <- sapply(lineages, function(l) l$get.location()$get.name())

  . <- sapply(comps, function(comp) {
    branching.time <- comp$get.branching.time()
    
    if (is.numeric(branching.time)) {
      s.name <- comp$get.source()
      if (is.element('R6', class(s.name))) {
        # used for unit testing
        source <- s.name
      }
      else {
        source <- comps[[s.name]]
        if (is.null(source)) {
          source <- comps[[paste(s.name, "1", sep="_")]]
          if (is.null(source)) {
            stop("Model missing Compartment ", s.name, " named as source for ", 
                 comp$get.name())
          }
        }  
      }
      
      # add transmission event to EventLogger object
      eventlog$add.event(
        type = 'transmission', 
        time = branching.time, 
        line1 = NA, 
        line2 = NA, 
        comp1 = comp$get.name(), 
        comp2 = source$get.name(),
        type1 = comp$get.type()$get.name(),
        type2 = source$get.type()$get.name()
        )
    }
    else {
      
      if (branching.time == 'NA' && comp$get.source() == 'NA') {
        # the index case will have no branching time specified
      } 
      else {
        stop ("Error in load.outer.tree(): cannot parse Compartment ",
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
#'    inner tree simulation (`attach.tree`).
#' 2. User manually inputs a host transmission tree into YAML format under 
#'    Compartments header (`load.outer.tree`).
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
#' @param model:  R6 object of class 'Model'
#' @param max.attempts:  maximum number of attempts to simulate outer events.
#'        Defaults to 5.
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
#' @seealso load.outer.tree
#' @export
sim.outer.tree <- function(model, max.attempts=5) {
  attempt <- 1
  
  while (attempt <= max.attempts) {
    run <- Run$new(model)
    
    # initialize a new object as return value
    eventlog <- run$get.eventlog()
    
    # store fixed sampling times of Lineages as specified by model
    eventlog$store.fixed.samplings(model$get.fixed.samplings())
    
    # sample events based on population dynamics and rates
    events <- .sample.outer.events(run)
    
    # assign events only to history of sampled Compartments
    resolved <- .assign.events(run, events)
    if (resolved) {
      break
    }
    
    attempt <- attempt + 1  # else try again
    warning("Failed to resolve outer tree, starting attempt ", attempt, immediate.=TRUE)
  }
  
  if (attempt == max.attempts) {
    stop("sim.outer.tree() failed to resolve an outer tree. ",
         "Try increasing the total population size or ")
  }

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
.get.rate.matrices <- function(types, current.time) {
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
      t.rate <- source$get.branching.rate(as.character(current.time), r.name)
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



#' .scale.contact.rates
#' 
#' Helper function for .sample.outer.events - calculate population rates 
#' given per-contact rate (\eqn{\beta}) and population sizes of source and
#' recipient populations, *i.e.*, \eqn{\beta S I}
#' 
#' @param mx:  rate matrix for source->recipient events
#' @param susceptible:  named vector, number of uninfected Compartments per Type
#' @param infected:  named vector, number of infected Compartments per Type
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
#' @param run:  R6 object of class Run
#' @param max.attempts: number of tries to sample events
#'        
#' @return data frame of outer tree events, each made up of: time, 
#'         event type, recipient CompartmentType and source CompartmentType
#'         
#' @keywords internal
.sample.outer.events <- function(run, max.attempts=10, chunk.size=1e4) {
  # extract objects from Run object
  comps <- run$get.compartments()
  types <- run$get.types()  # CompartmentType objects
  
  # record population totals and transmission rates for all Types
  init.conds <- run$get.initial.conds()
  popn.totals <- init.conds$size  # list of CompartmentType->size
  
  # extract transmission, migration and transition rates as convenient named matrices
  popn.rates <- .get.rate.matrices(types)
  
  attempt <- 1
  while (attempt < max.attempts) {
    
    # initialize populations at origin (named vectors)
    susceptible <- unlist(init.conds$size)
    infected <- sapply(susceptible, function(x) 0)
    susceptible[init.conds$indexType] <- susceptible[init.conds$indexType] - 1
    infected[init.conds$indexType] <- 1
    
    # pre-allocate outcome containers
    events <- data.frame(
      time=numeric(chunk.size), 
      event.type=character(chunk.size), 
      r.type=character(chunk.size), 
      s.type=character(chunk.size), 
      stringsAsFactors = FALSE
      )
    counts <- matrix(NA, nrow=chunk.size, ncol=length(c(susceptible, infected)))
    row.num <- 1
    
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
        
        # update counts
        r.name <- recipient$get.name()
        susceptible[r.name] <- susceptible[r.name] - 1
        infected[r.name] <- infected[r.name] + 1
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
          s.name <- source$get.name()
          r.name <- recipient$get.name()
          sus <- susceptible[s.name]
          inf <- infected[s.name]
          if (runif(1, max=sus+inf) < sus) {
            susceptible[s.name] <- sus - 1
            susceptible[r.name] <- susceptible[r.name] + 1
          } 
          else {
            infected[s.name] <- inf - 1
            infected[r.name] <- infected[r.name] + 1
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
      
      # append counts to outcome container after event
      #counts <- rbind(counts, c(susceptible, infected))
      counts[row.num, ] <- c(susceptible, infected)
      
      # append event
      events[row.num, ] <- c(current.time, event.type, recipient$get.name(), source$get.name())
      
      # go to next row
      row.num <- row.num + 1
      if (row.num > nrow(counts)) {
        # time to allocate more space
        counts <- rbind(counts, matrix(NA, nrow=chunk.size, ncol=ncol(counts)))
        events <- rbind(events, data.frame(
          time=numeric(chunk.size), 
          event.type=character(chunk.size), 
          r.type=character(chunk.size), 
          s.type=character(chunk.size), 
          stringsAsFactors = FALSE
        ))
      }
    } 
    
    # trim unused rows
    counts <- counts[!is.na(counts[,1]), ]
    events <- events[events$event.type != '', ]
    
    # check that there are enough Compartments of each Type (ignoring
    # sampling times)
    all.types <- sapply(comps, function(comp) comp$get.type()$get.name())
    checks <- sapply(1:length(infected), function(i) {
      ctype <- names(infected)[i]
      count <- infected[ctype]
      sum(all.types==ctype) <= count
    })
    
    if (all(checks)) {
      break
    }
    
    warning(".sample.outer.events failed one or more checks with ", 
            paste(infected, collapse=','), "\nAttempt ", attempt, 
            " of ", max.attempts, immediate.=TRUE)
    attempt <- attempt + 1
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
  run$set.counts(cbind(time=as.double(events$time), counts))  # store for plotting
  
  result <- cbind(events, counts)
  result$time <- as.double(result$time)
  result$event.type <- as.factor(result$event.type)
  result
}


#' .assign.events
#' 
#' Assignment of outer events proceeds in reverse (coalescent) time,
#' starting from the sampled Compartments.
#' @param run:  R6 object of class Run
#' @param events:  data frame from .sample.outer.events()
#' 
#' @return TRUE if events converge to index case, otherwise FALSE
#'         to trigger sim.outer.tree() to re-run analysis
#' @examples
#' # reproduce issue #103
#' set.seed(5)
#' run <- Run$new(model)
#' eventlog <- run$get.eventlog()
#' eventlog$store.fixed.samplings(model$get.fixed.samplings())
#' events <- .sample.outer.events(run)
#' @keywords internal
.assign.events <- function(run, events) {
  # TODO: access population dynamics through Run$get.counts() instead of 
  #       cbind'ed events
  eventlog <- run$get.eventlog()
  eventlog$clear.events()
  
  # start with sampled infected Compartments
  active <- run$get.compartments()
  types <- sapply(active, function(comp) comp$get.type()$get.name())
  samp.times <- sapply(active, function(comp) comp$get.sampling.time())
  
  # iterate through events in reverse time (start with most recent)
  for (row in seq(nrow(events), 1, -1)) {
    
    if (length(active) == 1) {
      # reached the root of sampled Compartments
      break
    }
    
    # go to the next event
    e <- events[row, ]  
    
    n.active.recipients <- sum(types==e$r.type)
    n.active.sources <- sum(types==e$s.type)
    
    if ( is.element(e$event.type, c('transmission', 'migration')) ) {
      
      # total numbers of infected Compartments of respective Types
      n.recipients <- as.numeric(e[paste('I.', e$r.type, sep='')])
      n.sources <- as.numeric(e[paste('I.', e$s.type, sep='')])
      
      # check if recipient is an active Compartment
      if (runif(1, max=n.recipients) < n.active.recipients) {
        # infection must precede first Lineage sampling time
        
        # FIXME: this is slow
        #eligible <- Filter(function(comp) is.na(comp$get.sampling.time()) || 
        #                     e$time > comp$get.sampling.time(), 
        #                   active[types==e$r.type])
        eligible <- active[which(types==e$r.type & (is.na(samp.times) || samp.times < e$time))]
          
        if (length(eligible) == 0) {
          stop("Error in .assign.events(): No eligible compartments for transmission event ", e)
        }
        recipient <- sample(eligible, 1)[[1]]
        
        # recipient cannot be its own source
        if (e$r.type == e$s.type) n.active.sources <- n.active.sources - 1
      }
      else {
        # recipient is NOT an active Compartment
        next
      }
      
        
      if (e$r.type == e$s.type) {
        # recipient cannot be its own source
        n.sources <- n.sources - 1
      }
      
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

        types <- c(types, source$get.type()$get.name())
        names(types)[length(types)] <- source$get.name()
        #types[[source$get.name()]] <- source$get.type()
        
        samp.times <- c(samp.times, NA)
        names(samp.times)[length(samp.times)] <- source$get.name()
        # if recipient is sampled, upgrade the source
        #if ( !recipient$is.unsampled() ) source$set.unsampled(FALSE)
      }
      
      # update event log
      eventlog$add.event(
        time = e$time,
        type = e$event.type,
        comp1 = recipient$get.name(),
        comp2 = source$get.name(),
        type1 = e$r.type,  # recipient (derived)
        type2 = e$s.type  # source (ancestral)
      )
      
      if (e$event.type == 'transmission') {
        # update recipient Compartment
        recipient$set.branching.time(e$time)          
        recipient$set.source(source)
        
        # remove recipient from active list
        active[recipient$get.name()] <- NULL
        types <- types[-which(names(types)==recipient$get.name())]
        samp.times <- samp.times[-which(names(samp.times)==recipient$get.name())]
        #types[recipient$get.name()] <- NULL
      }
    }
    
    else if (e$event.type == 'transition') {
      # either infected or uninfected Compartments may transition
      # look ahead to determine which
      if (row > 1) {
        next.e <- events[row-1, ]
        key <- paste0('S.', e$r.type)
        if (e[key] == next.e[key]) {
          # transitioning Compartment was infected
          n.recipients <- as.numeric(e[paste('I.', e$r.type, sep='')])
          
          if (runif(1, max=n.recipients) < n.active.recipients) {
            # transitioning Compartment is a sampled one
            compartment <- sample(active[types==e$r.type], 1)[[1]]
            
            # change the CompartmentType to ancestral ("source") state
            old.type <- compartment$get.type()$get.name()
            compartment$set.type(run$get.types()[[e$s.type]])
            
            eventlog$add.event(
              time = e$time,
              type = 'transition',
              comp1 = compartment$get.name(),
              type1 = old.type,
              type2 = compartment$get.type()$get.name()
            )
            
            # update <types> vector
            types[compartment$get.name()] <- compartment$get.type()$get.name()
          }
          # otherwise ignore transition of unsampled Compartment
        }
        # otherwise transitioning Compartment was not infected
      }
      # otherwise ignore first event (no more transmissions)
    }
    
    else {
      stop("Error in .assign.events: unrecognized event type ", 
           e$event.type)
    }
  }
  
  return(length(active) == 1)
}


