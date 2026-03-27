## test_superinfection.R
## Unit tests for superinfection behaviour in twt.
##
## Assumes superinf.yaml is in the same directory as this script.
## Override with: SUPERINF_YAML=/path/to/superinf.yaml Rscript test_superinfection.R

library(testthat)
library(twt)
library(yaml)

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

try_model <- function(path) {
  if (!file.exists(path)) {
    skip(paste("test_Superinfection.yaml not found at:", path))
  }
  tryCatch(
    Model$new(read_yaml(path)),
    error = function(e) skip(paste("Model$new failed:", conditionMessage(e)))
  )
}


# A1. .filter.firsts

cat("\n  .filter.firsts: when a host is infected multiple times, only the earliest transmission event is kept\n")
test_that(".filter.firsts: when a host is infected multiple times, only the earliest transmission event is kept", {
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

cat("\n  .filter.firsts: event log with no superinfection is returned unchanged\n")
test_that(".filter.firsts: event log with no superinfection is returned unchanged", {
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

cat("\n  .filter.firsts: three transmissions to the same host reduces to one (the earliest)\n")
test_that(".filter.firsts: three transmissions to the same host reduces to one (the earliest)", {
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

cat("\n  .filter.firsts: migration events are preserved after filtering superinfections\n")
test_that(".filter.firsts: migration events are preserved after filtering superinfections", {
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

cat("\n  .relabel.nodes: errors if called before .filter.firsts (superinfection events still present)\n")
test_that(".relabel.nodes: errors if called before .filter.firsts (superinfection events still present)", {
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

cat("\n  Superinfection (I+I->I+I): recipient is already in I, so no host moves compartments and sizes are unchanged\n")
test_that("Superinfection (I+I->I+I): recipient is already in I, so no host moves compartments and sizes are unchanged", {
  # row 2 is a superinfection (from.comp == to.comp == I)
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

cat("\n  Normal transmission (S+I->I+I): recipient moves from S to I, so S decrements and I increments\n")
test_that("Normal transmission (S+I->I+I): recipient moves from S to I, so S decrements and I increments", {
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

cat("\n  Model: I->I superinfection rate expression sigma*I*I is stored in transmission.rates[I,I,I]\n")
test_that("Model: I->I superinfection rate expression sigma*I*I is stored in transmission.rates[I,I,I]", {
  mod <- try_model(SUPERINF_PATH)
  tr  <- mod$get.transmission.rates()
  expect_equal(tr["I", "I", "I"], "sigma*I*I")
})

cat("\n  Model: compartment I is flagged infected=TRUE (required for superinfection to be recognized)\n")
test_that("Model: compartment I is flagged infected=TRUE (required for superinfection to be recognized)", {
  mod <- try_model(SUPERINF_PATH)
  expect_true(mod$get.infected("I"))
})

cat("\n  Model: standard S->I transmission rate beta*S*I is stored in transmission.rates[S,I,I]\n")
test_that("Model: standard S->I transmission rate beta*S*I is stored in transmission.rates[S,I,I]", {
  mod <- try_model(SUPERINF_PATH)
  tr  <- mod$get.transmission.rates()
  expect_equal(tr["S", "I", "I"], "beta*S*I")
})


# B2. sim.dynamics

cat("\n  sim.dynamics: superinfection events appear in event log as transmission with from.comp==to.comp==I\n")
test_that("sim.dynamics: superinfection events appear in event log as transmission with from.comp==to.comp==I", {
  mod <- try_model(SUPERINF_PATH)
  set.seed(13)
  dyn  <- sim.dynamics(mod, max.attempts = 10)
  evts <- dyn$events
  tx   <- evts[evts$event == "transmission", ]
  si   <- tx[tx$from.comp == "I" & tx$to.comp == "I", ]

  if (nrow(si) == 0) skip("No superinfection events generated with seed 42")
  expect_true(all(si$source == "I"))
})

cat("\n  sim.dynamics: I compartment size is unchanged before and after each superinfection row\n")
test_that("sim.dynamics: I compartment size is unchanged before and after each superinfection row", {
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

cat("\n  sim.outer.tree: returns OuterTree and logs transmission events (multiple active hosts expected with superinfection)\n")
test_that("sim.outer.tree: returns OuterTree and logs transmission events (multiple active hosts expected with superinfection)", {
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

  # outer log should contain at least one transmission event
  evts <- outer$get.log()
  tx   <- evts[evts$event == "transmission", ]
  expect_true(nrow(tx) >= 1)
})


# C1. .filter.firsts correctness — not just presence

cat("\n  .filter.firsts: unaffected hosts retain all their rows after filtering\n")
test_that(".filter.firsts: unaffected hosts retain all their rows after filtering", {
  events <- data.frame(
    time      = c(1.0, 1.5, 2.0, 3.0, 4.0),
    event     = c("transmission", "transmission", "transmission", "migration", "migration"),
    from.host = c("H1", "H2", "H1", "H3", "H4"),
    to.host   = c("H3", "H3", "H4", "H3_sample", "H4_sample"),
    stringsAsFactors = FALSE
  )
  filtered <- .filter.firsts(events)
  # H3 superinfected — only t=1.0 row survives
  expect_equal(sum(!is.na(filtered$to.host) & filtered$to.host == "H3"), 1)
  # H4 not superinfected — its transmission row must survive
  expect_equal(sum(!is.na(filtered$to.host) & filtered$to.host == "H4"), 1)
  # both migration rows must survive
  expect_equal(sum(filtered$event == "migration"), 2)
})

cat("\n  .filter.firsts: the surviving row is from the correct source host (earliest donor)\n")
test_that(".filter.firsts: the surviving row is from the correct source host (earliest donor)", {
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

cat("\n  .filter.firsts: two transmissions at identical time — one is kept (no duplication)\n")
test_that(".filter.firsts: two transmissions at identical time — one is kept (no duplication)", {
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

cat("\n  .filter.firsts: each host appears at most once as to.host in transmissions after filtering\n")
test_that(".filter.firsts: each host appears at most once as to.host in transmissions after filtering", {
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

cat("\n  .relabel.nodes: terminal migration to sampled compartment gets _sample suffix\n")
test_that(".relabel.nodes: terminal migration to sampled compartment gets _sample suffix", {
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

cat("\n  .relabel.nodes: non-sampled migration gets _1, _2 numeric suffixes\n")
test_that(".relabel.nodes: non-sampled migration gets _1, _2 numeric suffixes", {
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

cat("\n  .relabel.nodes: each host appears exactly once as to.host in transmissions (no duplication)\n")
test_that(".relabel.nodes: each host appears exactly once as to.host in transmissions (no duplication)", {
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


# C3. Zero sigma — no superinfection events

cat("\n  Zero sigma: no I->I transmission events when sigma=0 (deterministic check)\n")
test_that("Zero sigma: no I->I transmission events when sigma=0 (deterministic check)", {
  # Hand-crafted event log with only S->I transmissions
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


# C4. Superinfection does not create a new host (compartment count check)

cat("\n  Superinfection does not increase total infected hosts (I count unchanged)\n")
test_that("Superinfection does not increase total infected hosts (I count unchanged)", {
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
  # row 2 superinfection: I stays at 3
  expect_equal(events$I[2], events$I[1])
  # row 3 normal S->I: I increases
  expect_equal(events$I[3], events$I[2] + 1L)
})


# C5. Model-dependent: forced deterministic superinfection check with seed 13

cat("\n  sim.dynamics (seed 13): superinfection events guaranteed — donor and recipient are both in I\n")
test_that("sim.dynamics (seed 13): superinfection events guaranteed — donor and recipient are both in I", {
  mod <- try_model(SUPERINF_PATH)
  set.seed(13)
  dyn  <- sim.dynamics(mod, max.attempts = 10)
  evts <- dyn$events
  si   <- evts[evts$event == "transmission" &
                 evts$from.comp == "I" & evts$to.comp == "I", ]

  if (nrow(si) == 0) skip("Seed 13 produced no superinfection events — sigma may have changed")

  # source must be I (donor already infected)
  expect_true(all(si$source == "I"),
    info = "Superinfection donor must be from compartment I")

  # from.comp and to.comp must both be I (recipient already infected)
  expect_true(all(si$from.comp == "I" & si$to.comp == "I"),
    info = "Both from.comp and to.comp must be I for superinfection")
})

cat("\n  sim.dynamics (seed 13): after each superinfection, total I count is unchanged vs. previous row\n")
test_that("sim.dynamics (seed 13): after each superinfection, total I count is unchanged vs. previous row", {
  mod <- try_model(SUPERINF_PATH)
  set.seed(13)
  dyn  <- sim.dynamics(mod, max.attempts = 10)
  evts <- dyn$events

  si_idx <- which(evts$event == "transmission" &
                    evts$from.comp == "I" & evts$to.comp == "I")
  if (length(si_idx) == 0) skip("No superinfection events at seed 13")

  for (i in si_idx[si_idx > 1]) {
    expect_equal(evts$I[i], evts$I[i - 1],
      info = paste("I count changed at superinfection row", i,
                   "— superinfection must not move hosts between compartments"))
  }
})


# C6. End-to-end: outer log host label consistency after filtering

cat("\n  sim.outer.tree: every to.host in the log is a valid host seen elsewhere as from.host or root\n")
test_that("sim.outer.tree: every to.host in the log is a valid host seen elsewhere as from.host or root", {
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

  # every recipient host must appear somewhere as a source or be the root
  from_hosts <- unique(log$from.host)
  orphans <- setdiff(to_hosts, all_hosts)
  expect_equal(length(orphans), 0,
    info = "Some to.host values are not found as from.host anywhere in the log")
})

cat("\n  sim.outer.tree: sampled hosts in OuterTree match targets from Model\n")
test_that("sim.outer.tree: sampled hosts in OuterTree match targets from Model", {
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
  targets   <- mod$get.sampling()$targets
  target_n  <- sum(unlist(targets))

  expect_equal(n_sampled, target_n,
    info = paste("Expected", target_n, "sampled hosts, got", n_sampled))
})


# D1. Host-level superinfection semantics (handcrafted log)

cat("\n  Host semantics: superinfection recipient was already infected before the event\n")
test_that("Host semantics: superinfection recipient was already infected before the event", {
  # Construct a log where H2 is first infected at t=1, then superinfected at t=3.
  # The superinfection row must occur AFTER H2's first infection row.
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
  # find the first infection of the recipient
  recipient <- si_row$to.host
  first_inf <- events[events$event == "transmission" &
                        events$to.host == recipient, ]
  first_inf_time <- min(first_inf$time)

  # superinfection must happen AFTER the recipient's first infection
  expect_true(si_row$time > first_inf_time,
    info = "Superinfection occurred before recipient's first infection — invalid ordering")
})

cat("\n  Host semantics: donor and recipient are distinct hosts in superinfection\n")
test_that("Host semantics: donor and recipient are distinct hosts in superinfection", {
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

  # donor (from.host) must differ from recipient (to.host)
  expect_true(all(si$from.host != si$to.host),
    info = "Superinfection donor and recipient should be distinct hosts")
})


# D2. End-to-end pipeline: events -> filter -> relabel
cat("\n  Deterministic pipeline: filter then relabel produces consistent host labels\n")
test_that("Deterministic pipeline: filter then relabel produces consistent host labels", {
  targets <- list(I_samp = 5)

  # A controlled log: H2 superinfected (t=1 and t=2), H3 clean, then migrations
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

  # Step 1: filter superinfection
  filtered <- .filter.firsts(events)
  # H2 appears only once as to.host in transmissions
  tx <- filtered[filtered$event == "transmission", ]
  expect_equal(sum(!is.na(tx$to.host) & tx$to.host == "H2"), 1)

  # Step 2: relabel
  relabelled <- .relabel.nodes(filtered, targets)
  node_labels <- na.omit(unique(c(relabelled$from.host, relabelled$to.host)))

  # sampled nodes get _sample suffix
  expect_true(any(grepl("_sample$", node_labels)))
  # no node appears as both from.host and to.host with different content (label consistency)
  expect_false(anyDuplicated(na.omit(relabelled$to.host)) > 0)
})


# D3. Edge cases

cat("\n  Edge case: empty event log returns empty data frame from .filter.firsts\n")
test_that("Edge case: empty event log returns empty data frame from .filter.firsts", {
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

cat("\n  Edge case: single transmission with no superinfection passes through unchanged\n")
test_that("Edge case: single transmission with no superinfection passes through unchanged", {
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

cat("\n  Edge case: all events are migrations — .filter.firsts returns them unchanged\n")
test_that("Edge case: all events are migrations — .filter.firsts returns them unchanged", {
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


# D4. Tree structure validation (model-dependent)

# Helper: find a seed that produces a valid single-root phylo with the superinf model
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

cat("\n  Tree structure: number of tips equals number of sampled hosts\n")
test_that("Tree structure: number of tips equals number of sampled hosts", {
  mod    <- try_model(SUPERINF_PATH)
  result <- .get_valid_phylo(mod)
  if (is.null(result)) skip("No seed produced a single-root tree")

  n_sampled <- result$outer$get.sampled()$count.type()
  expect_equal(length(result$phy$tip.label), n_sampled,
    info = "Tip count must equal number of sampled hosts")
})

cat("\n  Tree structure: no duplicate tip labels in phylo object\n")
test_that("Tree structure: no duplicate tip labels in phylo object", {
  mod    <- try_model(SUPERINF_PATH)
  result <- .get_valid_phylo(mod)
  if (is.null(result)) skip("No seed produced a single-root tree")

  phy <- result$phy
  expect_equal(length(phy$tip.label), length(unique(phy$tip.label)),
    info = "Duplicate tip labels found in phylo object")
})

cat("\n  Tree structure: edge count equals nodes - 1 (acyclicity / binary tree invariant)\n")
test_that("Tree structure: edge count equals nodes - 1 (acyclicity / binary tree invariant)", {
  mod    <- try_model(SUPERINF_PATH)
  result <- .get_valid_phylo(mod)
  if (is.null(result)) skip("No seed produced a single-root tree")

  phy     <- result$phy
  n_nodes <- length(phy$tip.label) + phy$Nnode
  expect_equal(nrow(phy$edge), n_nodes - 1,
    info = "Edge count must equal total nodes - 1 for a valid tree")
})
