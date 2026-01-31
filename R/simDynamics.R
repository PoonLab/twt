#' update.rates
#' 
#' @param mod:  R6 object of class Model
#' @param envir:  R environment
#' @param reset:  bool, if TRUE then set compartments back to initial sizes
#'                defined in Model object
#' @return  list with rate vectors and matrices, keyed by event type
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
    mod$get.transmission.rates(), MARGIN=c(1,2),
    function(x) { eval(parse(text=x), envir=envir) })
  
  return(rates)
}


#' sim.dynamics
#'
#' Simulate population trajectories for all compartments forward in time.
#' Results are recorded as a data frame or exported to a CSV file with the 
#' following information:
#'   time:  numeric, time of event in simulation time
#'   event:  str, type of event (birth, death, migration or transmission)
#'   src:  str, Compartment originating the event
#'   dest:  str, for migration or transmission, the Compartment added to
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

  # instantiate model parameters
  e <- new.env()
  for (key in names(params)) {
    eval(parse(text=paste(key, "<-", params[[key]])), envir=e)
  }
  
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
        src=character(chunk.size), dest=character(chunk.size),
        stringsAsFactors=FALSE
      )
      row.num <- 1
    } else {
      # open file to write event log
      conn <- file(logfile, open='w')
      writeLines(text="time,event,src,dest", con=conn)
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
      if (event == 'birth') {  
        src <- sample(cnames, 1, prob=rates[['birth']])
        dest <- NA
        .plus.one(src, e)
        
      } else if (event == 'death') {  
        src <- sample(cnames, 1, prob=rates[['death']])
        dest <- NA
        .minus.one(src, e)
        
      } else if (event == 'migration') { 
        src <- sample(cnames, 1, prob=apply(rates[['migration']], 1, sum))
        dest <- sample(cnames, 1, prob=rates[['migration']][src,])
        .minus.one(src, e)
        .plus.one(dest, e)
        
      } else if (event == 'transmission') {
        src <- sample(cnames, 1, prob=apply(rates[['transmission']], 1, sum))
        dest <- sample(cnames, 1, prob=rates[['transmission']][src,])
        .minus.one(src, e)
        .plus.one(dest, e)
        
      } else {
        stop("This shouldn't be possible! Aughhhh!")
      }
      
      # log event
      if (is.null(logfile)) {
        events[row.num, ] <- c(params$simTime - cur.time, event, src, dest)
        row.num <- row.num+1
        if (row.num > nrow(events)) {
          # allocate more space
          events <- rbind(events, data.frame(
            time=numeric(chunk.size), event=character(chunk.size), 
            src=character(chunk.size), dest=character(chunk.size), 
            stringsAsFactors = FALSE
          ))
        }        
      } else {
        writeLines(
          text=paste(params$simTime-cur.time, event, src, dest, sep=","),
          con=conn)
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


#' Convenience functions
.plus.one <- function(comp, envir) {
  eval(parse(text=paste(comp, "<-", comp, "+1")), envir=envir)
}
.minus.one <- function(comp, envir) {
  eval(parse(text=paste(comp, "<-", comp, "-1")), envir=envir)
}


#' get.counts
#' 
#' Convert an event log into population size trajectories for every 
#' compartment in the model.  We omit these to keep the event log compact.
#' 
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
      src <- items[3]
      dest <- items[4]
      
      counts[row.num, ] <- counts[row.num-1, ]
      counts[row.num, 1] <- time
      if (event=="birth") {
        counts[row.num, src] <- counts[row.num, src] + 1
      } else if (event=="death") {
        counts[row.num, src] <- counts[row.num, src] - 1
      } else {
        # note src = dest if superinfection
        counts[row.num, src] <- counts[row.num, src] - 1
        counts[row.num, dest] <- counts[row.num, dest] + 1
      }
      row.num <- row.num + 1
    }
    close(con)
  } else if (is.data.frame(eventlog)) {
    for (i in 1:nrow(eventlog)) {
      event <- eventlog$event[i]
      src <- eventlog$src[i]
      dest <- eventlog$dest[i]
      
      counts[i+1, ] <- counts[i, ]
      counts[i+1, 1] <- eventlog$time[i]
      if (event=="birth") {
        counts[i+1, src] <- counts[i+1, src] + 1
      } else if (event=="death") {
        counts[i+1, src] <- counts[i+1, src] - 1
      } else {
        counts[i+1, src] <- counts[i+1, src] - 1
        counts[i+1, dest] <- counts[i+1, dest] + 1
      }
    }
  } else {
    stop("Unrecognized type for `eventlog`")
  }
  
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
                            lwd=2, ...) {
  k <- ncol(counts)-1
  cnames <- names(counts)[2:(k+1)]
  
  if (all(is.na(pal))) {
    pal <- hcl.colors(n=k, palette="Dark2")
  }
  
  par(mar=c(5,5,1,1))
  plot(counts$time, counts[,2], ylim=range(counts[,2:(k+1)]), col=pal[1],
       type='s', xlab=xlab, ylab=ylab, lwd=lwd, ...)
  if (k > 1) {
    for (i in 3:(k+1)) {
      lines(counts$time, counts[,i], type='s', col=pal[i-1], lwd=lwd)
    }
  }
}
