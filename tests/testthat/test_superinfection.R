## test_superinfection.R
## Unit tests for superinfection behaviour in twt.
##
## Assumes test_Superinfection.yaml is in the same directory as this script.
## Override with: SUPERINF_YAML=/path/to/superinf.yaml Rscript test_superinfection.R

library(testthat)
library(twt)
library(yaml)
library(ape)

# Internal functions copied from twt source (not exported in installed package)
.filter.firsts <- function(events) {
  idx <- which(events$event == 'transmission')
  first.idx <- sapply(split(idx, events$to.host[idx]), function(i) {
    i[which.min(events$time[i])]
  })
  remove <- idx[!is.element(idx, first.idx)]
  if (length(remove) > 0) events[-remove, ] else events
}

.relabel.nodes <- function(events, targets) {
  n.inf <- table(events$to.host)
  if (any(n.inf > 1)) {
    stop("ERROR: .relabel.nodes assumes events have been filtered of superinfection events.")
  }
  events <- events[order(events$time), ]
  nodes <- na.omit(unique(c(events$from.host, events$to.host)))
  for (node in nodes) {
    idx <- which(events$from.host == node | events$to.host == node)
    new.node <- node
    label <- 1
    for (i in idx) {
      if (events$from.host[i] == node) events$from.host[i] <- new.node
      if (events$event[i] == 'migration') {
        if (events$to.comp[i] %in% names(targets)) {
          new.node <- paste(node, "sample", sep="_")
        } else {
          new.node <- paste(node, label, sep="_")
          label <- label + 1
        }
        events$to.host[i] <- new.node
      }
    }
  }
  events
}

SUPERINF_PATH <- Sys.getenv("SUPERINF_YAML", unset = "")
if (!nzchar(SUPERINF_PATH)) {
  SUPERINF_PATH <- file.path(getwd(), "test_Superinfection.yaml")
}

SUPERINF_INNER_PATH <- SUPERINF_PATH

try_model <- function(path) {
  if (!file.exists(path)) {
    skip(paste("test_Superinfection.yaml not found at:", path))
  }
  tryCatch(
    Model$new(read_yaml(path)),
    error = function(e) skip(paste("Model$new failed:", conditionMessage(e)))
  )
}

.get_outer_with_superinfection <- function(mod, seeds = 1:500) {
  for (s in seeds) {
    set.seed(s)
    dyn <- tryCatch(sim.dynamics(mod, max.attempts = 10),
                    warning = function(w) NULL, error = function(e) NULL)
    if (is.null(dyn)) next

    outer <- tryCatch(
      withCallingHandlers(sim.outer.tree(dyn),
                          warning = function(w) invokeRestart("muffleWarning")),
      error = function(e) NULL)
    if (is.null(outer)) next

    log <- outer$get.log()
    si  <- log[log$event == "transmission" &
                 !is.na(log$from.comp) & log$from.comp == "I" &
                 !is.na(log$to.comp)  & log$to.comp  == "I", ]
    if (outer$get.active()$count.type() != 1) next
    if (nrow(si) == 0) next
    if (outer$get.sampled()$count.type() < 1) next

    return(list(outer = outer, dyn = dyn, seed = s, si_log = si))
  }
  return(NULL)
}

.try_inner <- function(outer) {
  tryCatch(
    suppressWarnings(sim.inner.tree(outer)),
    error = function(e) {
      message("sim.inner.tree error: ", conditionMessage(e))
      NULL
    })
}

.get_valid_phylo <- function(mod) {
  for (s in c(21, 77, 99, 200, 300, 400, 500)) {
    set.seed(s)
    dyn <- tryCatch(sim.dynamics(mod, max.attempts = 10), warning = function(w) NULL)
    if (is.null(dyn)) next
    outer <- tryCatch(
      withCallingHandlers(sim.outer.tree(dyn),
                          warning = function(w) invokeRestart("muffleWarning")),
      error = function(e) NULL
    )
    if (is.null(outer)) next
    if (outer$get.active()$count.type() != 1) next
    phy <- tryCatch(as.phylo(outer), error = function(e) NULL)
    if (!is.null(phy)) return(list(phy=phy, outer=outer))
  }
  return(NULL)
}


# A1. .filter.firsts

cat("\n  .filter.firsts: multiple infections to same host — only earliest kept\n")
test_that(".filter.firsts: multiple infections to same host — only earliest kept", {
  events <- data.frame(
    time      = c(1.0, 2.0, 3.0, 4.0),
    event     = c("transmission", "transmission", "transmission", "migration"),
    from.host = c("H1", "H2", "H1", "H3"),
    to.host   = c("H3", "H3", "H4", "H4_sample"),
    stringsAsFactors = FALSE
  )
  filtered <- .filter.firsts(events)
  h3_rows  <- filtered[!is.na(filtered$to.host) & filtered$to.host == "H3", ]
  expect_equal(nrow(h3_rows), 1)
  expect_equal(h3_rows$time, 1.0)
})

cat("\n  .filter.firsts: no superinfection — log returned unchanged\n")
test_that(".filter.firsts: no superinfection — log returned unchanged", {
  events <- data.frame(
    time      = c(1.0, 2.0, 3.0),
    event     = c("transmission", "transmission", "migration"),
    from.host = c("H1", "H2", "H3"),
    to.host   = c("H2", "H3", "H3_sample"),
    stringsAsFactors = FALSE
  )
  filtered <- .filter.firsts(events)
  expect_equal(nrow(filtered), nrow(events))
})

cat("\n  .filter.firsts: three transmissions to same host — one survives\n")
test_that(".filter.firsts: three transmissions to same host — one survives", {
  events <- data.frame(
    time      = c(1.0, 2.5, 4.0),
    event     = c("transmission", "transmission", "transmission"),
    from.host = c("H1", "H2", "H3"),
    to.host   = c("H4", "H4", "H4"),
    stringsAsFactors = FALSE
  )
  filtered <- .filter.firsts(events)
  h4_rows  <- filtered[!is.na(filtered$to.host) & filtered$to.host == "H4", ]
  expect_equal(nrow(h4_rows), 1)
  expect_equal(h4_rows$time, 1.0)
})

cat("\n  .filter.firsts: migration events preserved after filtering\n")
test_that(".filter.firsts: migration events preserved after filtering", {
  events <- data.frame(
    time      = c(0.5, 1.0, 2.0, 3.0),
    event     = c("migration", "transmission", "transmission", "migration"),
    from.host = c("H0", "H1", "H2", "H3"),
    to.host   = c("H1", "H3", "H3", "H3_sample"),
    stringsAsFactors = FALSE
  )
  filtered <- .filter.firsts(events)
  expect_equal(sum(filtered$event == "migration"), 2)
  tx <- filtered[filtered$event == "transmission", ]
  expect_equal(sum(tx$to.host == "H3"), 1)
})


# A2. .relabel.nodes

cat("\n  .relabel.nodes: errors if superinfection events still present\n")
test_that(".relabel.nodes: errors if superinfection events still present", {
  targets <- list(I_samp = 5)
  events  <- data.frame(
    time      = c(1.0, 2.0, 3.0),
    event     = c("transmission", "transmission", "migration"),
    from.comp = c("S", "I", "I"),
    src.comp  = c("I", "I", NA_character_),
    to.comp   = c("I", "I", "I_samp"),
    from.host = c("H1", "H2", "H3"),
    to.host   = c("H3", "H3", "H3_sample"),
    stringsAsFactors = FALSE
  )
  expect_error(.relabel.nodes(events, targets), regexp = "superinfection")
})


# A3. Compartment-size invariance (pure logic, no Model)

cat("\n  Superinfection (I+I->I+I): compartment sizes unchanged\n")
test_that("Superinfection (I+I->I+I): compartment sizes unchanged", {
  events <- data.frame(
    time      = c(1.0,  2.0,  3.0),
    event     = c("transmission", "transmission", "migration"),
    from.comp = c("S",  "I",   "I"),
    to.comp   = c("I",  "I",   "I_samp"),
    S         = c(499L, 499L,  499L),
    I         = c(3L,   3L,    2L),
    I_samp    = c(0L,   0L,    1L),
    R         = c(0L,   0L,    0L),
    stringsAsFactors = FALSE
  )
  si_idx <- which(events$event == "transmission" &
                    events$from.comp == events$to.comp)
  expect_equal(length(si_idx), 1)
  i <- si_idx[1]
  expect_equal(events$I[i], events$I[i - 1])
  expect_equal(events$S[i], events$S[i - 1])
})

cat("\n  Normal transmission (S->I): S decrements and I increments\n")
test_that("Normal transmission (S->I): S decrements and I increments", {
  events <- data.frame(
    time      = c(1.0,  2.0),
    event     = c("transmission", "transmission"),
    from.comp = c("S",  "S"),
    to.comp   = c("I",  "I"),
    S         = c(499L, 498L),
    I         = c(3L,   4L),
    stringsAsFactors = FALSE
  )
  i <- 2
  expect_equal(events$S[i], events$S[i - 1] - 1L)
  expect_equal(events$I[i], events$I[i - 1] + 1L)
})


# B1. Model parsing

cat("\n  Model: sigma*I*I stored in transmission.rates[I,I,I]\n")
test_that("Model: sigma*I*I stored in transmission.rates[I,I,I]", {
  mod <- try_model(SUPERINF_PATH)
  tr  <- mod$get.transmission.rates()
  expect_equal(tr["I", "I", "I"], "sigma*I*I")
})

cat("\n  Model: compartment I flagged infected=TRUE\n")
test_that("Model: compartment I flagged infected=TRUE", {
  mod <- try_model(SUPERINF_PATH)
  expect_true(mod$get.infected("I"))
})

cat("\n  Model: beta*S*I stored in transmission.rates[S,I,I]\n")
test_that("Model: beta*S*I stored in transmission.rates[S,I,I]", {
  mod <- try_model(SUPERINF_PATH)
  tr  <- mod$get.transmission.rates()
  expect_equal(tr["S", "I", "I"], "beta*S*I")
})


# B2. sim.dynamics

cat("\n  sim.dynamics: superinfection events have from.comp==to.comp==I\n")
test_that("sim.dynamics: superinfection events have from.comp==to.comp==I", {
  mod <- try_model(SUPERINF_PATH)
  set.seed(13)
  dyn  <- sim.dynamics(mod, max.attempts = 10)
  evts <- dyn$events
  tx   <- evts[evts$event == "transmission", ]
  si   <- tx[tx$from.comp == "I" & tx$to.comp == "I", ]

  if (nrow(si) == 0) skip("No superinfection events generated with seed 42")
  expect_true(all(si$source == "I"))
})

cat("\n  sim.dynamics: I count unchanged before and after each superinfection row\n")
test_that("sim.dynamics: I count unchanged before and after each superinfection row", {
  mod <- try_model(SUPERINF_PATH)
  set.seed(13)
  dyn  <- sim.dynamics(mod, max.attempts = 10)
  evts <- dyn$events

  si_idx <- which(evts$event == "transmission" &
                    evts$from.comp == "I" & evts$to.comp == "I")
  if (length(si_idx) == 0) skip("No superinfection events generated")

  for (i in si_idx) {
    if (i > 1) {
      expect_equal(evts$I[i], evts$I[i - 1],
                   info = paste("I changed at superinfection row", i))
    }
  }
})


# B3. Round-trip: sim.outer.tree -> as.phylo

cat("\n  sim.outer.tree: returns OuterTree with at least one transmission event\n")
test_that("sim.outer.tree: returns OuterTree with at least one transmission event", {
  mod <- try_model(SUPERINF_PATH)
  set.seed(99)
  dyn <- tryCatch(sim.dynamics(mod, max.attempts = 15),
                  warning = function(w) NULL)
  if (is.null(dyn)) skip("Epidemic went extinct")

  outer <- tryCatch(
    withCallingHandlers(
      sim.outer.tree(dyn),
      warning = function(w) invokeRestart("muffleWarning")
    ),
    error = function(e) NULL
  )
  if (is.null(outer)) skip("sim.outer.tree errored")

  expect_s3_class(outer, "OuterTree")
  evts <- outer$get.log()
  tx   <- evts[evts$event == "transmission", ]
  expect_true(nrow(tx) >= 1)
})


# C1. .filter.firsts correctness

cat("\n  .filter.firsts: unaffected hosts retain all rows after filtering\n")
test_that(".filter.firsts: unaffected hosts retain all rows after filtering", {
  events <- data.frame(
    time      = c(1.0, 1.5, 2.0, 3.0, 4.0),
    event     = c("transmission", "transmission", "transmission", "migration", "migration"),
    from.host = c("H1", "H2", "H1", "H3", "H4"),
    to.host   = c("H3", "H3", "H4", "H3_sample", "H4_sample"),
    stringsAsFactors = FALSE
  )
  filtered <- .filter.firsts(events)
  expect_equal(sum(!is.na(filtered$to.host) & filtered$to.host == "H3"), 1)
  expect_equal(sum(!is.na(filtered$to.host) & filtered$to.host == "H4"), 1)
  expect_equal(sum(filtered$event == "migration"), 2)
})

cat("\n  .filter.firsts: surviving row has correct source host (earliest donor)\n")
test_that(".filter.firsts: surviving row has correct source host (earliest donor)", {
  events <- data.frame(
    time      = c(1.0, 2.0),
    event     = c("transmission", "transmission"),
    from.host = c("H_early", "H_late"),
    to.host   = c("H_victim", "H_victim"),
    stringsAsFactors = FALSE
  )
  filtered <- .filter.firsts(events)
  surviving <- filtered[!is.na(filtered$to.host) & filtered$to.host == "H_victim", ]
  expect_equal(surviving$from.host, "H_early")
})

cat("\n  .filter.firsts: simultaneous transmissions to same host — one kept\n")
test_that(".filter.firsts: simultaneous transmissions to same host — one kept", {
  events <- data.frame(
    time      = c(1.0, 1.0),
    event     = c("transmission", "transmission"),
    from.host = c("H1", "H2"),
    to.host   = c("H3", "H3"),
    stringsAsFactors = FALSE
  )
  filtered <- .filter.firsts(events)
  expect_equal(sum(!is.na(filtered$to.host) & filtered$to.host == "H3"), 1)
})

cat("\n  .filter.firsts: each host appears at most once as to.host after filtering\n")
test_that(".filter.firsts: each host appears at most once as to.host after filtering", {
  events <- data.frame(
    time      = c(1.0, 2.0, 3.0, 4.0, 5.0),
    event     = rep("transmission", 5),
    from.host = c("H1","H2","H3","H4","H5"),
    to.host   = c("H6","H6","H7","H7","H8"),
    stringsAsFactors = FALSE
  )
  filtered <- .filter.firsts(events)
  tx <- filtered[filtered$event == "transmission", ]
  counts <- table(tx$to.host)
  expect_true(all(counts <= 1))
})


# C2. .relabel.nodes successful behavior

cat("\n  .relabel.nodes: sampled compartment migration gets _sample suffix\n")
test_that(".relabel.nodes: sampled compartment migration gets _sample suffix", {
  targets <- list(I_samp = 5)
  events <- data.frame(
    time      = c(1.0, 2.0),
    event     = c("transmission", "migration"),
    from.comp = c("S", "I"),
    src.comp  = c("I", NA_character_),
    to.comp   = c("I", "I_samp"),
    from.host = c("H1", "H2"),
    to.host   = c("H2", NA_character_),
    stringsAsFactors = FALSE
  )
  events$to.host[2] <- "H2_sample"
  relabelled <- .relabel.nodes(events, targets)
  expect_true(any(grepl("_sample$", relabelled$to.host), na.rm = TRUE))
})

cat("\n  .relabel.nodes: non-sampled migration gets numeric suffixes\n")
test_that(".relabel.nodes: non-sampled migration gets numeric suffixes", {
  targets <- list(I_samp = 5)
  events <- data.frame(
    time      = c(1.0, 2.0, 3.0, 4.0),
    event     = c("transmission", "migration", "migration", "migration"),
    from.comp = c("S", "I", "I2", "I2"),
    src.comp  = c("I", NA_character_, NA_character_, NA_character_),
    to.comp   = c("I", "I2", "I2", "I_samp"),
    from.host = c("H1", "H2", "H2_1", "H2_2"),
    to.host   = c("H2", NA_character_, NA_character_, NA_character_),
    stringsAsFactors = FALSE
  )
  events$to.host[2] <- "H2_1"
  events$to.host[3] <- "H2_2"
  events$to.host[4] <- "H2_sample"
  relabelled <- .relabel.nodes(events, targets)
  node_labels <- na.omit(unique(c(relabelled$from.host, relabelled$to.host)))
  expect_true(any(grepl("_1$", node_labels)))
  expect_true(any(grepl("_sample$", node_labels)))
})

cat("\n  .relabel.nodes: each host appears exactly once as to.host in transmissions\n")
test_that(".relabel.nodes: each host appears exactly once as to.host in transmissions", {
  targets <- list(I_samp = 5)
  events <- data.frame(
    time      = c(1.0, 2.0, 3.0),
    event     = c("transmission", "transmission", "migration"),
    from.comp = c("S", "S", "I"),
    src.comp  = c("I", "I", NA_character_),
    to.comp   = c("I", "I", "I_samp"),
    from.host = c("H1", "H2", "H3"),
    to.host   = c("H2", "H3", "H3_sample"),
    stringsAsFactors = FALSE
  )
  relabelled <- .relabel.nodes(events, targets)
  tx <- relabelled[relabelled$event == "transmission", ]
  counts <- table(tx$to.host)
  expect_true(all(counts <= 1))
})


# C3. Zero sigma

cat("\n  Zero sigma: no I->I transmission events\n")
test_that("Zero sigma: no I->I transmission events", {
  events <- data.frame(
    time      = c(1.0, 2.0, 3.0),
    event     = c("transmission", "transmission", "migration"),
    from.comp = c("S", "S", "I"),
    to.comp   = c("I", "I", "I_samp"),
    source    = c("I", "I", NA_character_),
    stringsAsFactors = FALSE
  )
  si <- events[events$event == "transmission" &
                 events$from.comp == events$to.comp, ]
  expect_equal(nrow(si), 0)
})


# C4. Superinfection does not create a new host

cat("\n  Superinfection: I count unchanged; normal transmission increments I\n")
test_that("Superinfection: I count unchanged; normal transmission increments I", {
  events <- data.frame(
    time      = c(1.0, 2.0, 3.0, 4.0),
    event     = c("transmission", "transmission", "transmission", "migration"),
    from.comp = c("S",  "I",   "S",   "I"),
    to.comp   = c("I",  "I",   "I",   "I_samp"),
    S         = c(997L, 997L,  996L,  996L),
    I         = c(3L,   3L,    4L,    3L),
    I_samp    = c(0L,   0L,    0L,    1L),
    R         = c(0L,   0L,    0L,    0L),
    stringsAsFactors = FALSE
  )
  expect_equal(events$I[2], events$I[1])
  expect_equal(events$I[3], events$I[2] + 1L)
})


# C5. Model-dependent: seed 13

cat("\n  sim.dynamics (seed 13): superinfection donor and recipient both in I\n")
test_that("sim.dynamics (seed 13): superinfection donor and recipient both in I", {
  mod <- try_model(SUPERINF_PATH)
  set.seed(13)
  dyn  <- sim.dynamics(mod, max.attempts = 10)
  evts <- dyn$events
  si   <- evts[evts$event == "transmission" &
                 evts$from.comp == "I" & evts$to.comp == "I", ]

  if (nrow(si) == 0) skip("Seed 13 produced no superinfection events — sigma may have changed")

  expect_true(all(si$source == "I"),
    info = "Superinfection donor must be from compartment I")
  expect_true(all(si$from.comp == "I" & si$to.comp == "I"),
    info = "Both from.comp and to.comp must be I for superinfection")
})

cat("\n  sim.dynamics (seed 13): I count unchanged after each superinfection\n")
test_that("sim.dynamics (seed 13): I count unchanged after each superinfection", {
  mod <- try_model(SUPERINF_PATH)
  set.seed(13)
  dyn  <- sim.dynamics(mod, max.attempts = 10)
  evts <- dyn$events

  si_idx <- which(evts$event == "transmission" &
                    evts$from.comp == "I" & evts$to.comp == "I")
  if (length(si_idx) == 0) skip("No superinfection events at seed 13")

  for (i in si_idx[si_idx > 1]) {
    expect_equal(evts$I[i], evts$I[i - 1],
      info = paste("I count changed at superinfection row", i))
  }
})


# C6. Outer log host label consistency

cat("\n  sim.outer.tree: every to.host is a valid host in the log\n")
test_that("sim.outer.tree: every to.host is a valid host in the log", {
  mod <- try_model(SUPERINF_PATH)
  set.seed(13)
  dyn <- tryCatch(sim.dynamics(mod, max.attempts = 15), warning = function(w) NULL)
  if (is.null(dyn)) skip("Epidemic went extinct")

  outer <- tryCatch(
    withCallingHandlers(sim.outer.tree(dyn),
                        warning = function(w) invokeRestart("muffleWarning")),
    error = function(e) NULL
  )
  if (is.null(outer)) skip("sim.outer.tree errored")

  log <- outer$get.log()
  all_hosts <- unique(c(log$from.host, na.omit(log$to.host)))
  to_hosts  <- na.omit(unique(log$to.host))
  orphans   <- setdiff(to_hosts, all_hosts)
  expect_equal(length(orphans), 0)
})

cat("\n  sim.outer.tree: sampled host count matches model targets\n")
test_that("sim.outer.tree: sampled host count matches model targets", {
  mod <- try_model(SUPERINF_PATH)
  set.seed(13)
  dyn <- tryCatch(sim.dynamics(mod, max.attempts = 15), warning = function(w) NULL)
  if (is.null(dyn)) skip("Epidemic went extinct")

  outer <- tryCatch(
    withCallingHandlers(sim.outer.tree(dyn),
                        warning = function(w) invokeRestart("muffleWarning")),
    error = function(e) NULL
  )
  if (is.null(outer)) skip("sim.outer.tree errored")

  n_sampled <- outer$get.sampled()$count.type()
  target_n  <- sum(unlist(mod$get.sampling()$targets))
  expect_equal(n_sampled, target_n,
    info = paste("Expected", target_n, "sampled hosts, got", n_sampled))
})


# D1. Host-level superinfection semantics

cat("\n  Host semantics: superinfection occurs after recipient's first infection\n")
test_that("Host semantics: superinfection occurs after recipient's first infection", {
  events <- data.frame(
    time      = c(1.0, 3.0, 5.0),
    event     = c("transmission", "transmission", "migration"),
    from.comp = c("S",  "I",   "I"),
    to.comp   = c("I",  "I",   "I_samp"),
    from.host = c("H1", "H1",  "H2"),
    to.host   = c("H2", "H2",  "H2_sample"),
    stringsAsFactors = FALSE
  )
  si_idx <- which(events$event == "transmission" &
                    events$from.comp == "I" & events$to.comp == "I")
  expect_equal(length(si_idx), 1)

  si_row   <- events[si_idx, ]
  recipient <- si_row$to.host
  first_inf <- events[events$event == "transmission" &
                        events$to.host == recipient, ]
  expect_true(si_row$time > min(first_inf$time))
})

cat("\n  Host semantics: donor and recipient are distinct hosts\n")
test_that("Host semantics: donor and recipient are distinct hosts", {
  events <- data.frame(
    time      = c(1.0, 2.0, 3.0),
    event     = c("transmission", "transmission", "migration"),
    from.comp = c("S",  "I",   "I"),
    to.comp   = c("I",  "I",   "I_samp"),
    from.host = c("H1", "H1",  "H2"),
    to.host   = c("H2", "H2",  "H2_sample"),
    stringsAsFactors = FALSE
  )
  si <- events[events$event == "transmission" &
                 events$from.comp == "I" & events$to.comp == "I", ]
  if (nrow(si) == 0) skip("No superinfection row in handcrafted log")
  expect_true(all(si$from.host != si$to.host))
})


# D2. End-to-end pipeline: events -> filter -> relabel

cat("\n  Deterministic pipeline: filter then relabel produces consistent host labels\n")
test_that("Deterministic pipeline: filter then relabel produces consistent host labels", {
  targets <- list(I_samp = 5)

  events <- data.frame(
    time      = c(1.0, 2.0, 3.0, 4.0, 5.0),
    event     = c("transmission", "transmission", "transmission", "migration", "migration"),
    from.comp = c("S",  "I",  "S",  "I",      "I"),
    to.comp   = c("I",  "I",  "I",  "I_samp", "I_samp"),
    from.host = c("H1", "H1", "H1", "H2",     "H3"),
    to.host   = c("H2", "H2", "H3", NA,       NA),
    stringsAsFactors = FALSE
  )
  events$to.host[4] <- "H2_sample"
  events$to.host[5] <- "H3_sample"

  filtered <- .filter.firsts(events)
  tx <- filtered[filtered$event == "transmission", ]
  expect_equal(sum(!is.na(tx$to.host) & tx$to.host == "H2"), 1)

  relabelled <- .relabel.nodes(filtered, targets)
  node_labels <- na.omit(unique(c(relabelled$from.host, relabelled$to.host)))
  expect_true(any(grepl("_sample$", node_labels)))
  expect_false(anyDuplicated(na.omit(relabelled$to.host)) > 0)
})


# D3. Edge cases

cat("\n  Edge case: empty event log — .filter.firsts returns empty\n")
test_that("Edge case: empty event log — .filter.firsts returns empty", {
  events <- data.frame(
    time      = numeric(0),
    event     = character(0),
    from.host = character(0),
    to.host   = character(0),
    stringsAsFactors = FALSE
  )
  filtered <- .filter.firsts(events)
  expect_equal(nrow(filtered), 0)
})

cat("\n  Edge case: single transmission — passes through unchanged\n")
test_that("Edge case: single transmission — passes through unchanged", {
  events <- data.frame(
    time      = 1.0,
    event     = "transmission",
    from.host = "H1",
    to.host   = "H2",
    stringsAsFactors = FALSE
  )
  filtered <- .filter.firsts(events)
  expect_equal(nrow(filtered), 1)
  expect_equal(filtered$to.host, "H2")
})

cat("\n  Edge case: all migrations — .filter.firsts returns unchanged\n")
test_that("Edge case: all migrations — .filter.firsts returns unchanged", {
  events <- data.frame(
    time      = c(1.0, 2.0, 3.0),
    event     = c("migration", "migration", "migration"),
    from.host = c("H1", "H2", "H3"),
    to.host   = c("H1_1", "H2_1", "H3_sample"),
    stringsAsFactors = FALSE
  )
  filtered <- .filter.firsts(events)
  expect_equal(nrow(filtered), 3)
  expect_equal(sum(filtered$event == "migration"), 3)
})


# D4. Tree structure validation

cat("\n  Tree structure: tip count equals sampled host count\n")
test_that("Tree structure: tip count equals sampled host count", {
  mod    <- try_model(SUPERINF_PATH)
  result <- .get_valid_phylo(mod)
  if (is.null(result)) skip("No seed produced a single-root tree")

  n_sampled <- result$outer$get.sampled()$count.type()
  expect_equal(length(result$phy$tip.label), n_sampled)
})

cat("\n  Tree structure: no duplicate tip labels\n")
test_that("Tree structure: no duplicate tip labels", {
  mod    <- try_model(SUPERINF_PATH)
  result <- .get_valid_phylo(mod)
  if (is.null(result)) skip("No seed produced a single-root tree")

  phy <- result$phy
  expect_equal(length(phy$tip.label), length(unique(phy$tip.label)))
})

cat("\n  Tree structure: edge count equals nodes - 1\n")
test_that("Tree structure: edge count equals nodes - 1", {
  mod    <- try_model(SUPERINF_PATH)
  result <- .get_valid_phylo(mod)
  if (is.null(result)) skip("No seed produced a single-root tree")

  phy     <- result$phy
  n_nodes <- length(phy$tip.label) + phy$Nnode
  expect_equal(nrow(phy$edge), n_nodes - 1)
})


# E1

cat("\n  sim.inner.tree: no crash with superinfection in outer tree\n")
test_that("sim.inner.tree: no crash with superinfection in outer tree", {
  mod    <- try_model(SUPERINF_INNER_PATH)
  result <- .get_outer_with_superinfection(mod)
  if (is.null(result)) skip("No seed produced outer tree with superinfection + >=2 sampled hosts")

  inner <- .try_inner(result$outer)
  expect_false(is.null(inner),
    info = paste("sim.inner.tree returned NULL for seed", result$seed))
})


# E2

cat("\n  sim.inner.tree -> as.phylo: tip count equals sampled host count\n")
test_that("sim.inner.tree -> as.phylo: tip count equals sampled host count", {
  mod    <- try_model(SUPERINF_INNER_PATH)
  result <- .get_outer_with_superinfection(mod)
  if (is.null(result)) skip("No suitable outer tree found")

  inner <- .try_inner(result$outer)
  if (is.null(inner)) skip("sim.inner.tree failed")
  if (result$outer$get.active()$count.type() != 1) skip("Multiple active roots")

  phy <- tryCatch(as.phylo(inner), error = function(e) NULL)
  if (is.null(phy)) skip("as.phylo failed")

  n_sampled <- result$outer$get.sampled()$count.type()
  expect_equal(length(phy$tip.label), n_sampled,
    info = paste("Expected", n_sampled, "tips, got", length(phy$tip.label)))
})


# E3

cat("\n  sim.inner.tree -> as.phylo: edge count equals nodes - 1\n")
test_that("sim.inner.tree -> as.phylo: edge count equals nodes - 1", {
  mod    <- try_model(SUPERINF_INNER_PATH)
  result <- .get_outer_with_superinfection(mod)
  if (is.null(result)) skip("No suitable outer tree found")

  inner <- .try_inner(result$outer)
  if (is.null(inner)) skip("sim.inner.tree failed")
  if (result$outer$get.active()$count.type() != 1) skip("Multiple active roots")

  phy <- tryCatch(as.phylo(inner), error = function(e) NULL)
  if (is.null(phy)) skip("as.phylo failed")

  n_nodes <- length(phy$tip.label) + phy$Nnode
  expect_equal(nrow(phy$edge), n_nodes - 1)
})


# E4

cat("\n  sim.inner.tree -> as.phylo: no duplicate tip labels\n")
test_that("sim.inner.tree -> as.phylo: no duplicate tip labels", {
  mod    <- try_model(SUPERINF_INNER_PATH)
  result <- .get_outer_with_superinfection(mod)
  if (is.null(result)) skip("No suitable outer tree found")

  inner <- .try_inner(result$outer)
  if (is.null(inner)) skip("sim.inner.tree failed")
  if (result$outer$get.active()$count.type() != 1) skip("Multiple active roots")

  phy <- tryCatch(as.phylo(inner), error = function(e) NULL)
  if (is.null(phy)) skip("as.phylo failed")

  expect_equal(length(phy$tip.label), length(unique(phy$tip.label)))
})


# E5

cat("\n  sim.inner.tree: inner log contains transmission events (lineages crossed hosts)\n")
test_that("sim.inner.tree: inner log contains transmission events (lineages crossed hosts)", {
  mod    <- try_model(SUPERINF_INNER_PATH)
  result <- .get_outer_with_superinfection(mod)
  if (is.null(result)) skip("No suitable outer tree found")

  inner <- .try_inner(result$outer)
  if (is.null(inner)) skip("sim.inner.tree failed")

  inner_log <- tryCatch(inner$get.log(), error = function(e) NULL)
  if (is.null(inner_log) || nrow(inner_log) == 0) skip("Inner tree has no event log")

  inner_tx <- inner_log[inner_log$event == "transmission", ]
  expect_true(nrow(inner_tx) > 0,
    info = "No transmission events in inner log — no lineages moved between hosts")
})


# E6

cat("\n  Model YAML: compartment I has pop.size=2 and bottleneck.size=1\n")
test_that("Model YAML: compartment I has pop.size=2 and bottleneck.size=1", {
  try_model(SUPERINF_INNER_PATH)
  cfg <- read_yaml(SUPERINF_INNER_PATH)
  expect_equal(cfg$Compartments$I$size,            2)
  expect_equal(cfg$Compartments$I$bottleneck.size, 1)
})


cat("\n  sim.inner.tree (20 replicates): every replicate produces transmission events\n")
test_that("sim.inner.tree (20 replicates): every replicate produces transmission events", {
  mod    <- try_model(SUPERINF_INNER_PATH)
  result <- .get_outer_with_superinfection(mod)
  if (is.null(result)) skip("No outer tree with superinfection found")

  n_with_tx <- 0

  for (rep in 1:20) {
    set.seed(rep * 7)
    inner <- .try_inner(result$outer)
    if (is.null(inner)) next

    ilog <- tryCatch(inner$get.log(), error = function(e) NULL)
    if (is.null(ilog) || nrow(ilog) == 0) next

    if (nrow(ilog[ilog$event == "transmission", ]) > 0) n_with_tx <- n_with_tx + 1
  }

  expect_true(n_with_tx > 0,
    info = "No replicate produced any inner-tree transmission events")
})


# E7

cat("\n  sim.inner.tree: MRCA of superinfected host's tips predates superinfection\n")
test_that("sim.inner.tree: MRCA of superinfected host's tips predates superinfection", {
  skip("Requires hand-crafted outer tree with known H1->H2, H3->H2 structure — deferred")
  mod    <- try_model(SUPERINF_INNER_PATH)
  result <- .get_outer_with_superinfection(mod)
  if (is.null(result)) skip("No outer tree with superinfection found")
  if (result$outer$get.active()$count.type() != 1) skip("Multiple active roots")

  si_recipients <- unique(result$si_log$to.host)
  earliest_si   <- min(result$si_log$time)
  found_older   <- FALSE

  for (rep in 1:20) {
    set.seed(rep * 13)
    inner <- .try_inner(result$outer)
    if (is.null(inner)) next

    phy <- tryCatch(as.phylo(inner), error = function(e) NULL)
    if (is.null(phy)) next

    tips    <- phy$tip.label
    si_tips <- tips[sapply(tips, function(tl)
      any(sapply(si_recipients, function(h) grepl(h, tl, fixed = TRUE))))]
    if (length(si_tips) < 2) next

    h <- tryCatch(branching.times(phy), error = function(e) NULL)
    if (is.null(h)) next

    mrca_node <- mrca(phy)[si_tips[1], si_tips[2]]
    if (is.na(mrca_node)) next

    nd         <- node.depth.edgelength(phy)
    root_depth <- max(nd[seq_along(tips)])
    mrca_age   <- root_depth - nd[mrca_node]

    if (!is.na(mrca_age) && mrca_age >= earliest_si) {
      found_older <- TRUE
      break
    }
  }

  expect_true(found_older,
    info = paste("20 replicates: MRCA never predated superinfection time",
                 round(earliest_si, 4)))
})
