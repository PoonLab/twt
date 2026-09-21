#' .compute.recomb.free.hosts
#'
#' Precompute, from the already-generated outer event log, which hosts can
#' safely skip recombination entirely: a host whose founding transmission
#' is a complete bottleneck (exactly one lineage enters), which never
#' receives a superinfection, and which has at most one directly-sampled
#' tip of its own can never produce an observable effect from internal
#' recombination -- everything that happens inside it either dies with it
#' or gets coalesced down to the single lineage that exits it anyway.
#'
#' This evaluates each host's bottleneck-size expression earlier than the
#' main loop otherwise would -- a deliberate, opt-in tradeoff (see
#' skip.recomb.free in sim.arg). It changes RNG call order relative to
#' skip.recomb.free=FALSE, so results from the same seed will differ, even
#' though the underlying distribution doesn't change.
#'
#' @keywords internal
#' @noRd
.compute.recomb.free.hosts <- function(events, mod, inner, envir=baseenv()) {
  is.infected <- mod$get.infected()

  trans      <- events[events$event == "transmission", ]
  is.si      <- is.infected[trans$from.comp]
  si.rows    <- trans[is.si, ]
  found.rows <- trans[!is.si, ]  # each host's own founding transmission

  si.count <- table(si.rows$to.host)

  mig <- events[events$event == "migration", ]
  if (nrow(mig) > 0) {
    is.target  <- vapply(mig$to.comp, inner$has.target, logical(1))
    samp.count <- table(mig$from.host[is.target])
  } else {
    samp.count <- table(character(0))
  }

  hosts    <- unique(found.rows$to.host)
  profiles <- vector("list", length(hosts))
  names(profiles) <- hosts

  for (h in hosts) {
    f.row <- found.rows[found.rows$to.host == h, ]
    if (nrow(f.row) != 1) {
      profiles[[h]] <- list(recomb.free = FALSE, bottleneck.size = NA)
      next
    }
    expr   <- mod$get.bottleneck.size(f.row$to.comp)
    b.size <- eval(parse(text = expr), envir = envir)

    si <- if (h %in% names(si.count))   si.count[[h]]   else 0
    ns <- if (h %in% names(samp.count)) samp.count[[h]] else 0

    profiles[[h]] <- list(
      recomb.free     = (b.size == 1 && si == 0 && ns <= 1),
      bottleneck.size = b.size
    )
  }
  profiles
}


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
#' @param skip.recomb.free  logical, if TRUE precompute which hosts can
#'   never have recombination observed in their sampled descendants
#'   (complete bottleneck, no superinfection, at most one sampled tip) and
#'   skip drawing recombination events for those hosts entirely. Default
#'   FALSE preserves exact prior behaviour/RNG stream.
#'
#' @return R6 object of class `InnerTree` (ARG events recorded in log)
#' @export
sim.arg <- function(outer, rho = 1e-4, seq.length = 9000L,
                     skip.recomb.free = FALSE) {
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

  host.profiles <- NULL
  if (skip.recomb.free) {
    host.profiles <- .compute.recomb.free.hosts(events, mod, inner, envir = env)
  }

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
        ev <- .draw.next.event(active, mod, rho, env, host.profiles = host.profiles)
        if (is.null(ev) || is.na(ev$dt) || ev$dt >= remaining) break

        remaining  <- remaining - ev$dt
        event.time <- e$time + remaining

        if (ev$type == "coalescent") {
          .do.coalescent(ev$host, inner, event.time, envir = env)
        } else {
          host.obj <- active$get.host.by.name(ev$host)
          p.size.expr <- mod$get.pop.size(host.obj$get.compartment())
          p.size.val <- eval(parse(text = p.size.expr), envir = env)
          bp <- .do.recombination(ev$host, ev$pathogen, inner, event.time,
                                  seq.length = seq.length, p.size = p.size.val)
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
      .do.infection(e, inner, envir = env, host.profiles = host.profiles)
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
#' @param host.profiles  optional named list from .compute.recomb.free.hosts;
#'   hosts flagged recomb.free=TRUE have their recombination rate zeroed.
#'
#' @return list(dt, type, host, pathogen) or NULL if no events possible
#' @keywords internal
#' @noRd
.draw.next.event <- function(active, mod, rho, envir = baseenv(),
                              host.profiles = NULL) {
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

    # recombination rate for this host (all lineages combined) -- skipped
    # entirely for hosts flagged recomb-free by skip.recomb.free (see
    # .compute.recomb.free.hosts): recombination inside such a host can
    # never be observed in the sampled tree, so drawing it would only
    # waste simulation time (and RNG draws) without changing any output.
    is.recomb.free <- !is.null(host.profiles) &&
      isTRUE(host.profiles[[h$get.name()]]$recomb.free)

    if (k >= 1 && rho > 0 && !is.recomb.free) {
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
                              seq.length = 9000L, p.size = NULL) {
  active <- inner$get.active()
  host   <- active$get.host.by.name(host.name)

  # initialize this host's FIXED lineage pool (sample
  # recombination parents from a fixed pool of p.size lineages, only
  # activating a previously-inactive one when chosen, instead of always
  # creating a brand-new lineage de novo. Total pool size never changes.)
  if (!host$is.pool.initialized()) {
    if (is.null(p.size)) {
      stop(".do.recombination: p.size must be provided to initialize ",
           "host's lineage pool on first use")
    }
    host$init.pool(p.size)
  }

  # ensure the recombining pathogen already occupies a pool slot -- if
  # this is the first pool-tracked event involving it, assign one now.
  # Can't exceed p.size though: if every slot's taken, this "new"
  # lineage must actually be the same individual as one we're already
  # tracking (needs a coalescent merge, not done yet) -- fail loudly
  # instead of quietly breaking the active <= pool.size invariant.
  if (is.na(pathogen$get.slot.id())) {
    if (host$count.active.slots() >= host$get.pool.size()) {
      stop(sprintf(
        "Host %s's lineage pool is full (%d/%d slots), can't assign a new one. ",
        host$get.name(), host$count.active.slots(), host$get.pool.size()),
        "Increase p.size, decrease rho, or add coalescent merging to the pool.")
    }
    host$activate.new.slot(pathogen)
  }
  own.slot <- pathogen$get.slot.id()

  # sample breakpoint uniformly across genome
  breakpoint <- sample.int(seq.length - 1L, 1L)
  pathogen$set.breakpoint(breakpoint)
  pathogen$set.start.time(time)

  # LEFT parent: continues in the SAME slot as the child
  parent.left <- inner$new.pathogen(time)
  parent.left$set.slot.id(own.slot)
  host$activate.slot(own.slot, parent.left)

  # RIGHT parent: sample from the fixed pool (excluding own slot)
  if (host$get.pool.size() <= 1) {
    # degenerate case: pool size 1, nothing else to sample from
    parent.right <- parent.left
  } else {
    draw <- host$sample.other.slot(exclude.slot.id = own.slot)
    if (draw$active) {
      parent.right <- draw$pathogen
    } else {
      parent.right <- inner$new.pathogen(time)
      host$activate.new.slot(parent.right)
    }
  }

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
  already.present <- any(sapply(host$get.pathogens(), function(p) {
    p$get.name() == parent.right$get.name()
  }))
  if (!already.present) host$add.pathogen(parent.right)

  # log the recombination event (breakpoint not stored in log -- fixed schema)
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
  # bp.by.child keeps names (per-child lookup); bps is deduped positions
  # only, for segment boundaries -- unique() strips names, so don't use
  # bps for per-child lookup
  bp.by.child <- unlist(arg.result$breakpoints)
  bps        <- sort(unique(bp.by.child))
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
  # precompute everything that does NOT depend on segment position --
  # these were being recomputed fresh inside the per-segment loop below
  # (unvectorized trans.rows loop, re-filtering recomb.rows per child per
  # segment) even though none of it references pos. With thousands of
  # segments this was the dominant cost (confirmed via Rprof: [.data.frame
  # subsetting was ~66% of total resolve.arg time). Only the recombination
  # parent CHOICE (left vs right) genuinely depends on pos.
  base.parent.map <- setNames(coal.rows$pathogen1, coal.rows$pathogen2)
  trans.p1 <- trans.rows$pathogen1
  trans.p2 <- trans.rows$pathogen2
  trans.valid <- !is.na(trans.p1) & !is.na(trans.p2)
  base.parent.map[trans.p2[trans.valid]] <- trans.p1[trans.valid]

  recomb.by.child <- split(recomb.rows, recomb.rows$pathogen1)
  recomb.children <- names(recomb.by.child)

  local.trees <- vector("list", length(starts))
  for (i in seq_along(starts)) {
    pos <- (starts[i] + ends[i]) / 2
    parent.map <- base.parent.map
    for (child in recomb.children) {
      child.rows <- recomb.by.child[[child]]
      if (nrow(child.rows) < 2) next
      bp <- bp.by.child[child]
      if (length(bp) == 0 || is.na(bp)) bp <- bps[1]
      parent <- if (pos <= bp) child.rows$pathogen2[1] else child.rows$pathogen2[2]
      parent.map[child] <- parent
    }

    # trace tips to root -- preallocate buffers and fill by index rather
    # than growing vectors with c() (still O(n^2) even without rbind/
    # data.frame -- each c() call still copies the whole vector so far).
    # A chain from any tip can't revisit a node (monotonic trace upward),
    # so length(all.paths) is a safe upper bound on any single chain's
    # length; multiply by n.tips for a safe total upper bound, then
    # truncate to what was actually used.
    max.steps <- length(all.paths) * length(tips)
    edge.parents  <- character(max.steps)
    edge.children <- character(max.steps)
    all.nodes.buf <- character(max.steps * 2L + length(tips))
    n.edges <- 0L
    n.nodes <- 0L
    for (tip in tips) {
      cur <- tip
      while (!is.na(parent.map[cur]) && !is.null(parent.map[cur])) {
        par <- parent.map[cur]
        n.edges <- n.edges + 1L
        edge.parents[n.edges]  <- par
        edge.children[n.edges] <- cur
        n.nodes <- n.nodes + 1L; all.nodes.buf[n.nodes] <- cur
        n.nodes <- n.nodes + 1L; all.nodes.buf[n.nodes] <- par
        cur <- par
      }
      n.nodes <- n.nodes + 1L
      all.nodes.buf[n.nodes] <- cur
    }
    edges <- data.frame(parent=edge.parents[seq_len(n.edges)],
                         child=edge.children[seq_len(n.edges)],
                         stringsAsFactors=FALSE)
    all.nodes <- all.nodes.buf[seq_len(n.nodes)]

    all.nodes <- unique(all.nodes)
    edges     <- unique(edges)
    root.node <- all.nodes[!all.nodes %in% edges$child]
    if (length(root.node) > 1) root.node <- root.node[1]

    # O(1) child lookup via precomputed split, instead of scanning the
    # full edges data.frame on every call -- to.newick() below calls
    # this once per node in the tree, so the previous O(n) scan made
    # the whole traversal O(n^2) (confirmed dominant cost via Rprof).
    children.map <- split(edges$child, edges$parent)
    get.children <- function(node) {
      ch <- children.map[[node]]
      if (is.null(ch)) character(0) else ch
    }

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
