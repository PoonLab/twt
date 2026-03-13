#' sim.dynamics
#'
#' Simulate population trajectories for all compartments forward in time.
#' Results are recorded as a data frame or exported to a CSV file with the 
#' following information:
#'   time:  numeric, time of event in simulation time
#'   event:  str, type of event (birth, death, migration or transmission)
#'   from.comp:  str, Compartment originating the event
#'   to.comp:  str, for migration or transmission, the Compartment added to
#'   source:  str, for transmission only, compartment of the infecting host
#' 
#' @param mod:  R6 object of class Model
#' @param logfile:  str, optional path to write event log (default: NULL)
#' @param counted:  bool, if FALSE then do not record compartment sizes over 
#'                  time to reduce memory consumption.  These counts can be 
#'                  regenerated later from logged events by get.counts(). 
#'                  (default: TRUE)
#' @param max.attempts:  int, number of attempts to simulate trajectories 
#'                       (default: 3)
#' @param chunk.size:  int, number of rows to allocate to data.frame that 
#'                     stores events (default: 1e4)
#' @return data.frame
#' @examples
#' require(twt)
#' settings <- yaml.load_file("examples/SIRS_serial.yaml")
#' mod <- Model$new(settings)
#' sim.dynamics(mod, logfile="eventlog.csv")
#'
#' @export
sim.dynamics <- function(mod, logfile=NA, counted=TRUE, max.attempts=3, 
                              chunk.size=1e4) {
  # check that input is a Model object
  if ( !is.R6(mod) | !is.element("Model", class(mod)) ) {
    stop("Error: input `mod` must be an R6 object of class `Model`")
  }
  
  # unpack the Model object
  params <- mod$get.parameters()
  cnames <- mod$get.compartments()
  k <- length(cnames)
  sampling <- mod$get.sampling()
  
  # initialize model parameters in new environment
  envir <- .init.model(mod)
  
  attempt <- 1
  while (attempt <= max.attempts) {
    
    # initialize Compartments and set rate matrices
    rates <- .update.rates(mod, envir, reset=TRUE)
    
    # are there other stopping criteria?
    targets <- NULL
    if (sampling$mode == "compartment") {
      targets <- sampling$targets
    }
    
    # pre-allocate containers for logging events
    events <- data.frame(
      time=numeric(chunk.size), event=character(chunk.size),
      from.comp=character(chunk.size), to.comp=character(chunk.size),
      source=character(chunk.size),
      stringsAsFactors=FALSE
    )
    # check if we are recording compartment sizes over time
    if (counted) {
      for (cn in cnames) { events[[cn]] <- integer(chunk.size) }
    }
    
    row.num <- 1
    
    if (!is.na(logfile)) {
      # open file to write event log
      conn <- file(logfile, open='w')
      writeLines(text="time,event,from.comp,to.comp,source", con=conn)
    }
    
    ### MAIN SIMULATION LOOP ###
    cur.time <- params$simTime
    while (cur.time >= 0) {
      
      rate.sums <- sapply(rates, sum)
      total.rate <- sum(rate.sums)
      if (total.rate == 0) {
        break  # nothing can occur, end simulation
      }
      
      # sample waiting time to next event
      wait <- rexp(1, total.rate)
      cur.time <- cur.time - wait
      if (cur.time <= 0) {
        break  # end simulation
      }
      
      # which event?
      event <- sample(names(rates), 1, prob=rate.sums/total.rate)
      to.comp <- NA  # default values
      source <- NA
      if (event == 'birth') {  
        from.comp <- sample(cnames, 1, prob=rates[['birth']])
        .plus.one(from.comp, envir)
        
      } else if (event == 'death') {  
        from.comp <- sample(cnames, 1, prob=rates[['death']])
        .minus.one(from.comp, envir)
        
      } else if (event == 'migration') { 
        from.comp <- sample(cnames, 1, prob=apply(rates[['migration']], 1, sum))
        to.comp <- sample(cnames, 1, prob=rates[['migration']][from.comp, ])
        .minus.one(from.comp, envir)
        .plus.one(to.comp, envir)
        
      } else if (event == 'transmission') {
        # from compartment
        from.comp <- sample(cnames, 1, prob=apply(
          rates[['transmission']], 1, sum))
        
        # to compartment
        probs.given.from <- apply(
          rates[['transmission']][from.comp, , ], 1, sum)
        to.comp <- sample(cnames, 1, prob=probs.given.from)
        
        # infected by (source)
        source <- sample(cnames, 1, 
                         prob=rates[['transmission']][from.comp, to.comp, ])
        
        if (!mod$get.infected(from.comp)) {
          # recipient leaves source (uninfected) compartment
          .minus.one(from.comp, envir)
          .plus.one(to.comp, envir)  # enters destination (infected) compartment
        }
        # otherwise it is superinfection and compartment sizes do not change!
      } else {
        stop("This shouldn't be possible! Aughhhh!")
      }
      
      # log event
      row <- list(time=params$simTime-cur.time, event=event, 
                  from.comp=from.comp, to.comp=to.comp, source=source)
      if (counted) {
        counts <- sapply(cnames, function(cn) {
          eval(parse(text=cn), envir=envir)
        })
        row <- c(row, counts)
      }
      events[row.num, ] <- row
      
      row.num <- row.num + 1
      if (row.num > nrow(events)) {
        # allocate more space
        chunk <- data.frame(
          time=numeric(chunk.size), event=character(chunk.size), 
          from.comp=character(chunk.size), to.comp=character(chunk.size), 
          source=character(chunk.size),
          stringsAsFactors = FALSE
        )
        if (counted) {
          for (cn in cnames) { chunk[[cn]] <- integer(chunk.size) }
        }
        events <- rbind(events, chunk)
      }        
      
      if (!is.na(logfile)) {
        writeLines(text=paste(event, sep=","), con=conn)
        flush(conn)
      }
      
      # check stopping criteria (serial sampling)
      if (event=='migration' & sampling$mode == 'compartment') {
        sampled <- sapply(names(targets), function(cn) {
          eval(parse(text=cn), envir=envir)
        })
        if (all(sampled >= targets)) {
          break  # no need to simulate further
        }
      }
      
      # a change in compartment size could affect any of these rates
      rates <- .update.rates(mod, envir)
    }
    # end while, cur.time exceeds limit
    
    # check sample sizes
    if (is.null(targets)) {
      break
    } else {
      sampled <- sapply(names(targets), function(cn) {
        eval(parse(text=cn), envir=envir)
      })
      if (all(sampled >= targets)) { break }
    }
    
    message("Failed sample size requirements (attempt ", attempt, "/", 
                max.attempts, ")")
    attempt <- attempt + 1
    
    if (!is.na(logfile)) { close(conn) }
  }
  
  if (attempt > max.attempts) {
    warning("All attempts failed with insufficient sample size. You may ",
            "need to increase simulation time (simTime), or the sampling ",
            "or transmission rates in the model.")
  }

  # generate return value
  dynamics <- list()
  dynamics$events <- events[events$event != '', ]  # discard unused rows
  dynamics$model <- mod
  dynamics$logfile <- logfile
  dynamics$is.counted <- counted
  class(dynamics) <- c("dynamics")
  return(dynamics)
}


#' .init.model
#' Create a new environment and instantiate the model parameters
#' @param mod:  R6 object of class Model
#' @return  R environment
#' @keywords internal
#' @noRd
.init.model <- function(mod) {
  params <- mod$get.parameters()
  envir <- new.env()
  for (key in names(params)) {
    eval(parse(text=paste(key, "<-", params[[key]])), envir=envir)
  }
  return(envir)
}


#' .update.rates
#' 
#' @param mod:  R6 object of class Model
#' @param envir:  R environment
#' @param reset:  bool, if TRUE then set compartments back to initial sizes
#'                defined in Model object
#' @return  list with rate vectors and matrices, keyed by event type
#' @keywords internal
#' @noRd
.update.rates <- function(mod, envir, reset=FALSE) {
  
  if (reset) {
    # initialize compartments and set initial sizes
    init.sizes <- mod$get.init.sizes()
    for (cn in mod$get.compartments()) {
      eval(parse(text=paste(cn, "<-", init.sizes[[cn]])), envir=envir)
    }    
  }
  
  rates <- list()
  
  rates[['birth']] <- sapply(mod$get.birth.rates(), function(x) {
    eval(parse(text=x), envir=envir)
  })
  
  rates[['death']] <- sapply(mod$get.death.rates(), function(x) {
    eval(parse(text=x), envir=envir)
  })
  
  rates[['migration']] <- apply(
    mod$get.migration.rates(), MARGIN=c(1,2), 
    function(x) { eval(parse(text=x), envir=envir) })
  
  rates[['transmission']] <- apply(
    mod$get.transmission.rates(), MARGIN=c(1,2,3),
    function(x) { eval(parse(text=x), envir=envir) })
  
  return(rates)
}


#' Convenience functions
#' @keywords internal
#' @noRd
.plus.one <- function(comp, envir) {
  eval(parse(text=paste(comp, "<-", comp, "+1")), envir=envir)
}

#' @keywords internal
#' @noRd
.minus.one <- function(comp, envir) {
  eval(parse(text=paste(comp, "<-", comp, "-1")), envir=envir)
}


#' get.counts
#' 
#' Convert an event log into population size trajectories for every 
#' compartment in the model.  The user may choose to omit these counts to 
#' reduce the size of the eventlog `data.frame`.
#' 
#' @param dynamics:  S3 object of class 'dynamics'
#' @return dynamics updated with compartment sizes
#' 
#' @export
get.counts <- function(dynamics) {
  if (dynamics$is.counted) {
    warning("dynamics object is already counted!")
    return (dynamics)
  }
  eventlog <- dynamics$events
  mod <- dynamics$model
  
  # prepare data frame to track compartment sizes
  n <- nrow(eventlog) + 1
  counts <- data.frame(time=numeric(n))
  cnames <- mod$get.compartments()
  for (cn in cnames) { counts[[cn]] <- integer(n) }
  counts$time <- 0
  counts[1, cnames] <- mod$get.init.sizes()

  for (i in 1:nrow(eventlog)) {
    event <- eventlog$event[i]
    from.comp <- eventlog$from.comp[i]
    to.comp <- eventlog$to.comp[i]
    
    counts[i+1, ] <- counts[i, ]
    counts[i+1, 1] <- as.numeric(eventlog$time[i])
    if (event=="birth") {
      counts[i+1, from.comp] <- counts[i+1, from.comp] + 1
    } else if (event=="death") {
      counts[i+1, from.comp] <- counts[i+1, from.comp] - 1
    } else if (
      event=="migration" | 
      (event=="transmission" & !mod$get.infected(from.comp))
      ) {
      counts[i+1, from.comp] <- counts[i+1, from.comp] - 1
      counts[i+1, to.comp] <- counts[i+1, to.comp] + 1
    }
  }
  
  dynamics$events <- cbind(eventlog, counts[-1, -1])
  dynamics$is.counted <- TRUE
  return(dynamics)
}


#' Generic S3 print function for dynamics objects
#' @export
#' @noRd
print.dynamics <- function(obj) {
  cat("twt Dynamics object\n")
  cat("  Counted: ", obj$is.counted, "\n")
  cat(" ", nrow(obj$events), "events:")
  print(table(obj$events$event))
  t.range <- range(obj$events$time)
  cat("  Time range: ", t.range[1], "->", t.range[2], "\n")
}


#' plot.dynamics
#' 
#' S3 method for visualizing the simulated population dynamics from 
#' sim.dynamics().
#' @param counts:  S3 object of class `dynamics`
#' @param pal:  character, palette for drawing lines; defaults to Dark2
#' @param xlab:  character, label for x-axis (default: 'Time')
#' @param ylab:  character, label for y-axis (default: 'Frequency')
#' @param lwd:  numeric, line width for plotting
#' @param bty:  character, box type (default: 'n')
#' @param ylim:  numeric, y-axis limits (defaults to whole range)
#' @param mar:  numeric, plot margins (defaults to c(5,5,1,5))
#' @param ...:  additional arguments passed to base plot() function
#' 
#' @export
plot.dynamics <- function(dynamics, pal=NA, xlab="Time", ylab="Frequency", 
                            lwd=2, bty='n', ylim=NA, mar=c(5,5,1,5), ...) {
  if (!dynamics$is.counted) {
    warning("plot.dynamics requires counts, calling get.counts()")
    dynamics <- get.counts(dynamics)
  }
  
  # unpack the input object
  mod <- dynamics$model
  events <- dynamics$events
  
  # extract compartment size trajectories
  cnames <- mod$get.compartments()
  counts <- events[c('time', cnames)]
  row0 <- c(time=0, mod$get.init.sizes())
  counts <- rbind(row0, counts)
  counts$time <- as.numeric(counts$time)
  
  k <- length(cnames)
  if (all(is.na(pal))) {
    pal <- hcl.colors(n=k, palette="Dark2")  # generate colour palette
  }
  
  if (any(is.na(ylim))) {
    ylim <- range(counts[,2:(k+1)])  # default range
  }
  
  par(mar=mar)  # plot margins
  plot(counts$time, counts[,2], ylim=ylim, col=pal[1],
       type='s', xlab=xlab, ylab=ylab, lwd=lwd, bty=bty, ...)
  if (k > 1) {
    for (i in 2:(k+1)) {
      lines(counts$time, counts[,i], type='s', col=pal[i-1], lwd=lwd)
      text(1.05*max(counts$time), counts[nrow(counts),i], label=names(counts)[i],
           xpd=NA, adj=0, col=pal[i-1], cex=0.7)
    }
  }
}
