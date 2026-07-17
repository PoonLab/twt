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
  breakpoints <- list()  # named list: pathogen name -> breakpoint position

  # iterate through outer tree events in reverse time
  events     <- outer$get.log()
  events$time <- as.numeric(events$time)
  events     <- events[order(events$time, decreasing = TRUE), ]

  time.delta <- -diff(events$time)
  row        <- 1

  while (row <= nrow(events)) {
    e       <- events[row, ]
    remaining <- ifelse(row == 1, NA, time.delta[row - 1])

    # draw competing coalescence and recombination waiting times
    if (row > 1 && inner$n.active() > 0) {
      repeat {
        if (inner$n.active() > 500) {
          warning("Active lineage count exceeded 500 -- aborting inner simulation")
          break
        }
        ev <- .draw.next.event(active, mod, rho, env)
        if (is.null(ev) || is.na(ev$dt) || ev$dt >= remaining) break

        remaining  <- remaining - ev$dt
        event.time <- e$time + remaining

        if (ev$type == "coalescent") {
          .do.coalescent(ev$host, inner, event.time, envir = env)
        } else {
          bp <- .do.recombination(ev$host, ev$pathogen, inner, event.time,
                                  seq.length = seq.length)
          breakpoints[[bp$child]] <- bp$position
        }
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

  return(list(inner = inner, breakpoints = breakpoints))
}


#' Draw the next event using the Gillespie algorithm.
#'
#' Computes all event rates into a single vector, draws one exponential
#' waiting time from the total rate, then selects which event occurred
#' using the broken-stick (normalized rates) method.
#'
#' Events tracked:
#'   - Coalescence: one entry per host with 2+ lineages, rate = choose(k,2)/Ne
#'   - Recombination: one entry per host, rate = k * rho (per-lineage rate)
#'
#' @param active   HostSet of active hosts
#' @param mod      Model R6 object
#' @param rho      recombination rate per lineage per unit time
#' @param envir    environment for rate evaluation
#'
#' @return list(dt, type, host, pathogen) or NULL if no events possible
#' @keywords internal
#' @noRd
.draw.next.event <- function(active, mod, rho, envir = baseenv()) {
  hosts      <- active$get.hosts()
  rates      <- numeric(0)
  event.list <- list()

  for (h in hosts) {
    k    <- h$count.pathogens()
    comp <- h$get.compartment()

    # coalescence rate for this host (requires 2+ lineages)
    if (k >= 2) {
      expr <- mod$get.coalescent.rate(comp)
      rate <- eval(parse(text = expr), envir = envir)
      if (rate > 0) {
        rates      <- c(rates, choose(k, 2) * rate)
        event.list <- c(event.list, list(
          list(type = "coalescent", host = h$get.name(), pathogen = NA)
        ))
      }
    }

    # recombination rate for this host (all lineages combined)
    if (k >= 1 && rho > 0) {
      rates      <- c(rates, k * rho)
      event.list <- c(event.list, list(
        list(type = "recombination", host = h$get.name(), pathogen = NULL)
      ))
    }
  }

  if (length(rates) == 0 || sum(rates) == 0) return(NULL)

  # single exponential draw from total rate
  total.rate <- sum(rates)
  dt         <- rexp(1, total.rate)

  # broken-stick selection: choose event proportional to its rate
  u       <- runif(1) * total.rate
  chosen  <- which(cumsum(rates) >= u)[1]
  ev      <- event.list[[chosen]]

  # for recombination, choose a random lineage from the host now
  if (ev$type == "recombination") {
    h     <- active$get.host.by.name(ev$host)
    paths <- h$get.pathogens()
    ev$pathogen <- paths[[sample(length(paths), 1)]]
  }

  return(list(dt = dt, type = ev$type, host = ev$host, pathogen = ev$pathogen))
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
  pathogen$set.breakpoint(breakpoint)

  # end the current lineage at this recombination event
  pathogen$set.start.time(time)

  # create two parental lineages — left and right of breakpoint
  parent.left  <- inner$new.pathogen(time)
  parent.right <- inner$new.pathogen(time)

  # record parent-child relationships (recombination has two parents)
  parent.left$add.child(pathogen)
  parent.right$add.child(pathogen)
  pathogen$add.parent(parent.left)
  pathogen$add.parent(parent.right)

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

  return(list(child = pathogen$get.name(), position = breakpoint,
              left = parent.left$get.name(), right = parent.right$get.name()))
}


#' resolve.arg
#'
#' Resolve an ARG into local trees per genomic segment by traversing
#' the event log.  For each segment, at each recombination node the
#' appropriate parent lineage is followed (left if position <= breakpoint,
#' right otherwise).
#'
#' @param arg.result  list returned by sim.arg ($inner, $breakpoints)
#' @param seq.length  integer, total genome length in bp
#'
#' @return list:
#'   - segments: data.frame(start, end)
#'   - local.trees: list of phylo objects, one per segment
#'   - breakpoints: sorted numeric vector
#' @export
resolve.arg <- function(arg.result, seq.length = 9000L) {
  inner      <- arg.result$inner
  bps        <- sort(unlist(arg.result$breakpoints))
  log        <- inner$get.log()
  log$time   <- as.numeric(log$time)

  # genomic segments
  starts <- c(1L, as.integer(bps) + 1L)
  ends   <- c(as.integer(bps), as.integer(seq.length))

  coal.rows   <- log[log$event == "coalescent", ]
  recomb.rows <- log[log$event == "recombination", ]
  samp.rows   <- log[log$event == "sampling", ]
  trans.rows  <- log[log$event == "transmission", ]

  # tip names (sampled pathogens)
  tips <- unique(samp.rows$pathogen1)

  # node creation time from pathogen registry
  all.paths <- inner$get.all.pathogens()
  node.created <- setNames(
    sapply(all.paths, function(p) p$get.end.time()),
    names(all.paths)
  )

  # for each segment, trace from tips to root
  local.trees <- vector("list", length(starts))
  for (i in seq_along(starts)) {
    pos <- (starts[i] + ends[i]) / 2
    # parent map for this segment
    parent.map <- setNames(coal.rows$pathogen1, coal.rows$pathogen2)
    for (r in seq_len(nrow(trans.rows))) {
      p2 <- trans.rows$pathogen2[r]
      p1 <- trans.rows$pathogen1[r]
      if (!is.na(p2) && !is.na(p1)) parent.map[p2] <- p1
    }
    recomb.children <- unique(recomb.rows$pathogen1)
    for (child in recomb.children) {
      child.rows <- recomb.rows[recomb.rows$pathogen1 == child, ]
      if (nrow(child.rows) < 2) next
      bp <- bps[names(bps) == child]
      if (length(bp) == 0) bp <- bps[1]
      parent <- if (pos <= bp) child.rows$pathogen2[1] else child.rows$pathogen2[2]
      parent.map[child] <- parent
    }

    # trace tips to root
    all.nodes <- character(0)
    edges     <- data.frame(parent=character(0), child=character(0),
                             stringsAsFactors=FALSE)
    for (tip in tips) {
      cur <- tip
      while (!is.na(parent.map[cur]) && !is.null(parent.map[cur])) {
        par <- parent.map[cur]
        edges <- rbind(edges, data.frame(parent=par, child=cur,
                                          stringsAsFactors=FALSE))
        all.nodes <- c(all.nodes, cur, par)
        cur <- par
      }
      all.nodes <- c(all.nodes, cur)
    }

    all.nodes <- unique(all.nodes)
    edges     <- unique(edges)
    root.node <- all.nodes[!all.nodes %in% edges$child]
    if (length(root.node) > 1) root.node <- root.node[1]

    get.children <- function(node) edges$child[edges$parent == node]

    get.time <- function(node) {
      t <- node.created[node]
      if (is.na(t) || is.null(t)) return(0)
      as.numeric(t)
    }

    to.newick <- function(node, parent.t = NULL) {
      children <- get.children(node)
      t <- get.time(node)
      bl <- if (!is.null(parent.t)) max(0, t - parent.t) else 0
      if (length(children) == 0) {
        return(paste0(node, ":", round(bl, 6)))
      } else {
        subtrees <- sapply(children, to.newick, parent.t = t)
        return(paste0("(", paste(subtrees, collapse=","), ")",
                      node, ":", round(bl, 6)))
      }
    }

    nwk <- paste0(to.newick(root.node), ";")
    phy <- tryCatch(ape::read.tree(text=nwk), error=function(e) NULL)

    local.trees[[i]] <- list(
      start=starts[i], end=ends[i],
      newick=nwk, phylo=phy
    )
  }

  return(list(
    segments    = data.frame(start=starts, end=ends, stringsAsFactors=FALSE),
    local.trees = local.trees,
    breakpoints = bps
  ))
}
