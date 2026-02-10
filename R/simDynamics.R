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
#' @param logfile:  str, path to write event log; if not specified, then 
#'                  this function will return a data.frame
#' @param max.attempts:  int, number of attempts to simulate trajectories
#' @param chunk.size:  int, number of rows to allocate to data.frame that 
#'                     stores events
#' @return data.frame if user does not specify logfile
#' @examples
#' require(twt)
#' settings <- yaml.load_file("examples/SIRS_serial.yaml")
#' mod <- Model$new(settings)
#' sim.dynamics(mod, logfile="eventlog.csv")
#'
#' @export
sim.dynamics <- function(mod, logfile=NULL, max.attempts=3, 
                              chunk.size=1e4) {
  # unpack the Model object
  params <- mod$get.parameters()
  cnames <- mod$get.compartments()
  k <- length(cnames)
  sampling <- mod$get.sampling()
  
  # initialize model parameters in new environment
  e <- .init.model(mod)
  
  attempt <- 1
  while (attempt <= max.attempts) {
    
    # initialize Compartments and set rate matrices
    rates <- .update.rates(mod, e, reset=TRUE)
    
    # are there other stopping criteria?
    targets <- NULL
    if (sampling$mode == "compartment") {
      targets <- sampling$targets
    }
    
    if (is.null(logfile)) {
      # pre-allocate containers for logging events
      events <- data.frame(
        time=numeric(chunk.size), event=character(chunk.size),
        from.comp=character(chunk.size), to.comp=character(chunk.size),
        source=character(chunk.size),
        stringsAsFactors=FALSE
      )
      row.num <- 1
    } else {
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
        .plus.one(from.comp, e)
        
      } else if (event == 'death') {  
        from.comp <- sample(cnames, 1, prob=rates[['death']])
        .minus.one(from.comp, e)
        
      } else if (event == 'migration') { 
        from.comp <- sample(cnames, 1, prob=apply(rates[['migration']], 1, sum))
        to.comp <- sample(cnames, 1, prob=rates[['migration']][from.comp, ])
        .minus.one(from.comp, e)
        .plus.one(to.comp, e)
        
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
          .minus.one(from.comp, e)
          .plus.one(to.comp, e)  # enters destination (infected) compartment
        }
        # otherwise it is superinfection and compartment sizes do not change!
      } else {
        stop("This shouldn't be possible! Aughhhh!")
      }
      
      # log event
      if (is.null(logfile)) {
        events[row.num, ] <- c(
          params$simTime - cur.time, event, from.comp, to.comp, source)
        
        row.num <- row.num+1
        if (row.num > nrow(events)) {
          # allocate more space
          events <- rbind(events, data.frame(
            time=numeric(chunk.size), event=character(chunk.size), 
            from.comp=character(chunk.size), to.comp=character(chunk.size), 
            source=character(chunk.size),
            stringsAsFactors = FALSE
          ))
        }        
      } else {
        writeLines(text=paste(
          params$simTime-cur.time, event, from.comp, to.comp, source, 
          sep=","), con=conn)
        flush(conn)
      }
      
      # check stopping criteria
      if (event=='migration') {
        sampled <- sapply(names(targets), function(cn) {
          eval(parse(text=cn), envir=e)
        })
        if (all(sampled >= targets)) {
          break  # no need to simulate further
        }
      }
      
      # a change in compartment size could affect any of these rates
      rates <- .update.rates(mod, e)
    }
    # end while, cur.time exceeds limit
    
    # check sample sizes
    if (is.null(targets)) {
      break
    } else {
      sampled <- sapply(names(targets), function(cn) {
        eval(parse(text=cn), envir=e)
      })
      if (all(sampled >= targets)) { break }
    }
    
    message("Failed sample size requirements (attempt ", attempt, "/", 
                max.attempts, ")")
    attempt <- attempt + 1
    if (!is.null(logfile)) { close(conn) }
  }
  
  if (attempt > max.attempts) {
    warning("All attempts failed with insufficient sample size. You may ",
            "need to increase simulation time (simTime), or the sampling ",
            "or transmission rates in the model.")
  }
  
  # discard unused rows and return
  if (is.null(logfile)) {
    return(events[events$event != '', ])
  } else {
    close(conn)
    return(NULL)
  }
}


#' .init.model
#' Create a new environment and instantiate the model parameters
#' @param mod:  R6 object of class Model
#' @return  R environment
#' @keywords internal
.init.model <- function(mod) {
  params <- mod$get.parameters()
  e <- new.env()
  for (key in names(params)) {
    eval(parse(text=paste(key, "<-", params[[key]])), envir=e)
  }
  return(e)
}


#' .update.rates
#' 
#' @param mod:  R6 object of class Model
#' @param envir:  R environment
#' @param reset:  bool, if TRUE then set compartments back to initial sizes
#'                defined in Model object
#' @return  list with rate vectors and matrices, keyed by event type
#' @keywords internal
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
.plus.one <- function(comp, envir) {
  eval(parse(text=paste(comp, "<-", comp, "+1")), envir=envir)
}
#' @keywords internal
.minus.one <- function(comp, envir) {
  eval(parse(text=paste(comp, "<-", comp, "-1")), envir=envir)
}


#' get.counts
#' 
#' Convert an event log into population size trajectories for every 
#' compartment in the model.  We omit these to keep the event log compact.
#' @param eventlog:  data.frame or character, result from sim.dynamics()
#' @param mod:  R6 object of class Model
#' @param chunk.size:  integer, used to allocate new rows for output data frame
#' @return data.frame of compartment sizes over time
#' @examples
#' counts <- get.counts("eventlog.csv", mod)
#' plot(counts)  # S3 method
#' 
#' @export
get.counts <- function(eventlog, mod, chunk.size=100) {
  # prepare output container
  cnames <- mod$get.compartments()
  counts <- data.frame(
    time = numeric(chunk.size)
  )
  for (cn in cnames) {
    counts[[cn]] <- integer(chunk.size)
  }
  counts$time[1] <- 0
  counts[1, cnames] <- mod$get.init.sizes()
  row.num <- 2
  
  if (is.character(eventlog)) {
    # input is a CSV file
    con <- file(eventlog, open='r')
    while (length(line <- readLines(con, n=1, warn=FALSE)) > 0) {
      items <- strsplit(line, ",")[[1]]
      time <- suppressWarnings(as.numeric(items[1]))
      if (is.na(time)) {
        next  # probably header line
      }
      event <- items[2]
      from.comp <- items[3]
      to.comp <- items[4]
      
      counts[row.num, ] <- counts[row.num-1, ]
      counts[row.num, 1] <- time
      if (event=="birth") {
        counts[row.num, from.comp] <- counts[row.num, from.comp] + 1
        
      } else if (event=="death") {
        counts[row.num, from.comp] <- counts[row.num, from.comp] - 1
        
      } else if (
        event=='migration' | 
        (event=='transmission' & !mod$get.infected(from.comp))
        ) {
        # not a superinfection
        counts[row.num, from.comp] <- counts[row.num, from.comp] - 1
        counts[row.num, to.comp] <- counts[row.num, to.comp] + 1          
      }
      row.num <- row.num + 1
    }
    close(con)
  } else if (is.data.frame(eventlog)) {
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
  } else {
    stop("Unrecognized type for `eventlog`")
  }
  
  # remove unused rows
  idx <- (counts$time==0)
  idx[1] <- FALSE
  counts <- counts[!idx, ]
  
  class(counts) <- c('twt.counts', 'data.frame')
  return(counts)
}


#' plot.twt.counts
#' 
#' S3 method for visualizing the simulated population dynamics from 
#' simulate.dynamics().
#' @param counts:  S3 object of class twt.counts
#' @param pal:  character, palette for drawing lines; defaults to Dark2
#' @param xlab:  character, label for x-axis (default: 'Time')
#' @param ylab:  character, label for y-axis (default: 'Frequency')
#' @param lwd:  numeric, line width for plotting
#' @param ...:  additional arguments passed to base plot() function
#' 
#' @export
plot.twt.counts <- function(counts, pal=NA, xlab="Time", ylab="Frequency", 
                            lwd=2, bty='n', ...) {
  k <- ncol(counts)-1
  cnames <- names(counts)[2:(k+1)]
  
  if (all(is.na(pal))) {
    pal <- hcl.colors(n=k, palette="Dark2")
  }
  
  par(mar=c(5,5,1,5))
  plot(counts$time, counts[,2], ylim=range(counts[,2:(k+1)]), col=pal[1],
       type='s', xlab=xlab, ylab=ylab, lwd=lwd, bty=bty, ...)
  if (k > 1) {
    for (i in 2:(k+1)) {
      lines(counts$time, counts[,i], type='s', col=pal[i-1], lwd=lwd)
      text(1.05*max(counts$time), counts[nrow(counts),i], label=names(counts)[i],
           xpd=NA, adj=0, col=pal[i-1], cex=0.7)
    }
  }
}
