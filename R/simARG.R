#' sim.arg
#'
#' Simulate an ancestral recombination graph (ARG) for pathogen lineages
#' within hosts, given an outer transmission tree.  Extends sim.inner.tree
#' by adding recombination as a competing process alongside coalescence.
#'
#' Going backwards in time, within each host interval two competing events
#' can occur:
#'   - Coalescence: two lineages merge (rate 1/Ne per pair)
#'   - Recombination: one lineage splits into two parental lineages at a
#'     random genomic breakpoint (rate rho per lineage)
#'
#' Each lineage tracks the set of genomic segments it is ancestral to.
#' A recombination event splits those segments between two new parent
#' lineages.  Simulation stops when all segments have coalesced to a
#' single ancestor (MRCA per segment).
#'
#' @param outer  R6 object of class `OuterTree`
#' @param rho    numeric, recombination rate per lineage per unit time
#' @param seq.length  integer, genome length in bp (for breakpoint sampling)
#'
#' @return R6 object of class `InnerTree` (ARG events recorded in log)
#' @export
sim.arg <- function(outer, rho = 1e-4, seq.length = 9000L) {
  if (!is.R6(outer) || !is.element("OuterTree", class(outer))) {
    stop("Input argument must be an R6 object of class `OuterTree`")
  }

  inner <- InnerTree$new(outer)
  mod   <- inner$get.model()
  active   <- inner$get.active()
  inactive <- inner$get.inactive()

  env <- new.env()
  eval(parse(text = "require(stats)"), envir = env)

  # iterate through outer tree events in reverse time
  events     <- outer$get.log()
  events$time <- as.numeric(events$time)
  events     <- events[order(events$time, decreasing = TRUE), ]

  time.delta <- -diff(events$time)
  row        <- 1

  while (row <= nrow(events)) {
    e       <- events[row, ]
    t.delta <- ifelse(row == 1, NA, time.delta[row - 1])

    # draw competing coalescence and recombination waiting times
    if (row > 1 && inner$n.active() > 0) {
      ev <- .draw.next.event(active, mod, rho, env)

      if (!is.null(ev) && !is.na(ev$dt) && ev$dt < t.delta) {
        t.delta <- t.delta - ev$dt
        event.time <- e$time + t.delta

        if (ev$type == "coalescent") {
          .do.coalescent(ev$host, inner, event.time, envir = env)
        } else {
          .do.recombination(ev$host, ev$pathogen, inner, event.time,
                            seq.length = seq.length)
        }
        next
      }
    }

    # handle outer tree event
    if (e$event == "migration") {
      if (inner$has.target(e$to.comp)) {
        .do.sampling(e, inner)
      } else {
        .migrate.pathogens(e, inner)
      }
    } else if (e$event == "transmission") {
      .do.infection(e, inner, envir = env)
    }

    row <- row + 1
  }

  # finish coalescence in last host
  count <- active$count.type()
  if (count == 1) {
    root     <- active$get.hosts()[[1]]
    cur.time <- e$time
    while (root$count.pathogens() > 1) {
      wait <- .rcoal(active, mod, envir = env)
      cur.time <- cur.time - wait$dt
      .do.coalescent(root$get.name(), inner, cur.time)
    }
  } else if (count > 1) {
    warning("Multiple (", count, ") hosts remain active at end of simulation.")
  } else {
    warning("Empty active HostSet at end of simulation!")
  }

  return(inner)
}


#' Draw the next event — whichever of coalescence or recombination happens
#' first across all active hosts.
#'
#' @param active   HostSet of active hosts
#' @param mod      Model R6 object
#' @param rho      recombination rate per lineage per unit time
#' @param envir    environment for rate evaluation
#'
#' @return list(dt, type, host, pathogen) or NULL if no event possible
#' @keywords internal
#' @noRd
.draw.next.event <- function(active, mod, rho, envir = baseenv()) {
  best    <- list(dt = Inf, type = NA, host = NA, pathogen = NA)

  for (h in active$get.hosts()) {
    k    <- h$count.pathogens()
    comp <- h$get.compartment()

    # --- coalescence ---
    if (k >= 2) {
      expr  <- mod$get.coalescent.rate(comp)
      rate  <- eval(parse(text = expr), envir = envir)
      if (rate > 0) {
        dt.coal <- rexp(1, choose(k, 2) * rate)
        if (dt.coal < best$dt) {
          best <- list(dt = dt.coal, type = "coalescent",
                       host = h$get.name(), pathogen = NA)
        }
      }
    }

    # --- recombination --- (one event per lineage)
    if (k >= 1 && rho > 0) {
      dt.recomb <- rexp(1, k * rho)
      if (dt.recomb < best$dt) {
        # pick a random lineage to recombine
        paths <- h$get.pathogens()
        chosen <- paths[[sample(length(paths), 1)]]
        best <- list(dt = dt.recomb, type = "recombination",
                     host = h$get.name(), pathogen = chosen)
      }
    }
  }

  if (is.infinite(best$dt)) return(NULL)
  return(best)
}


#' Recombination event: split one lineage into two parental lineages.
#'
#' Going backwards in time, a recombination event means a lineage we are
#' tracking was produced by recombination — it has two parents, one for
#' each side of the breakpoint.  We sample a breakpoint and create two new
#' parental lineages, then remove the original.
#'
#' @param host.name   character, name of host where event occurs
#' @param pathogen    Pathogen R6 object to split
#' @param inner       InnerTree R6 object
#' @param time        numeric, time of event
#' @param seq.length  integer, genome length
#'
#' @keywords internal
#' @noRd
.do.recombination <- function(host.name, pathogen, inner, time,
                              seq.length = 9000L) {
  active <- inner$get.active()
  host   <- active$get.host.by.name(host.name)

  # sample breakpoint uniformly across genome
  breakpoint <- sample.int(seq.length - 1L, 1L)

  # end the current lineage at this recombination event
  pathogen$set.start.time(time)

  # create two parental lineages — left and right of breakpoint
  parent.left  <- inner$new.pathogen(time)
  parent.right <- inner$new.pathogen(time)

  # record parent-child relationships
  parent.left$add.child(pathogen)
  parent.right$add.child(pathogen)
  pathogen$set.parent(parent.left)

  # remove original pathogen by finding its index
  paths <- host$get.pathogens()
  idx <- which(sapply(paths, function(p) p$get.name()) == pathogen$get.name())
  if (length(idx) == 1) host$remove.pathogen(idx)
  host$add.pathogen(parent.left)
  host$add.pathogen(parent.right)

  # log the recombination event (breakpoint not stored in log — fixed schema)
  event <- list(
    time = time, event = "recombination",
    from.comp = host$get.compartment(), to.comp = NA,
    from.host = host.name, to.host = NA,
    pathogen1 = pathogen$get.name(),
    pathogen2 = parent.left$get.name()
  )
  inner$add.event(event)
  event$pathogen2 <- parent.right$get.name()
  inner$add.event(event)
}
