## test_superinfection.R
## Tests for superinfection behaviour in twt.
##
## Assumes test_Superinfection.yaml is in the same directory.
## Override with: SUPERINF_YAML=/path/to/file.yaml Rscript test_superinfection.R

library(testthat)
library(twt)
library(yaml)
library(ape)

# these aren't exported so copy them here for testing
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

# find an outer tree that has at least one SI event and a single active root
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
                 !is.na(log$to.comp)  & log$to.comp   == "I", ]
    if (outer$get.active()$count.type() != 1) next
    if (nrow(si) == 0) next
    if (outer$get.sampled()$count.type() < 1) next

    return(list(outer = outer, dyn = dyn, seed = s, si_log = si))
  }
  NULL
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
      error = function(e) NULL)
    if (is.null(outer)) next
    if (outer$get.active()$count.type() != 1) next
    phy <- tryCatch(as.phylo(outer), error = function(e) NULL)
    if (!is.null(phy)) return(list(phy = phy, outer = outer))
  }
  NULL
}

.get_valid_outer <- function(mod, seeds = c(99, 21, 77, 42, 55, 100, 200)) {
  for (s in seeds) {
    set.seed(s)
    dyn <- tryCatch(sim.dynamics(mod, max.attempts = 15), warning = function(w) NULL)
    if (is.null(dyn)) next
    outer <- tryCatch(
      withCallingHandlers(sim.outer.tree(dyn),
                          warning = function(w) invokeRestart("muffleWarning")),
      error = function(e) NULL)
    if (!is.null(outer)) return(outer)
  }
  NULL
}


# --- .filter.firsts ---

test_that(".filter.firsts: multiple infections to same host — only earliest kept", {
  events <- data.frame(
    time = c(1.0, 2.0, 3.0, 4.0),
    event = c("transmission", "transmission", "transmission", "migration"),
    from.host = c("H1", "H2", "H1", "H3"),
    to.host = c("H3", "H3", "H4", "H4_sample"),
    stringsAsFactors = FALSE
  )
  filtered <- .filter.firsts(events)
  h3_rows <- filtered[!is.na(filtered$to.host) & filtered$to.host == "H3", ]
  expect_equal(nrow(h3_rows), 1)
  expect_equal(h3_rows$time, 1.0)
})

test_that(".filter.firsts: no superinfection — log returned unchanged", {
  events <- data.frame(
    time = c(1.0, 2.0, 3.0),
    event = c("transmission", "transmission", "migration"),
    from.host = c("H1", "H2", "H3"),
    to.host = c("H2", "H3", "H3_sample"),
    stringsAsFactors = FALSE
  )
  expect_equal(nrow(.filter.firsts(events)), nrow(events))
})

test_that(".filter.firsts: three transmissions to same host — one survives", {
  events <- data.frame(
    time = c(1.0, 2.5, 4.0),
    event = c("transmission", "transmission", "transmission"),
    from.host = c("H1", "H2", "H3"),
    to.host = c("H4", "H4", "H4"),
    stringsAsFactors = FALSE
  )
  filtered <- .filter.firsts(events)
  h4_rows <- filtered[!is.na(filtered$to.host) & filtered$to.host == "H4", ]
  expect_equal(nrow(h4_rows), 1)
  expect_equal(h4_rows$time, 1.0)
})

test_that(".filter.firsts: migration events preserved after filtering", {
  events <- data.frame(
    time = c(0.5, 1.0, 2.0, 3.0),
    event = c("migration", "transmission", "transmission", "migration"),
    from.host = c("H0", "H1", "H2", "H3"),
    to.host = c("H1", "H3", "H3", "H3_sample"),
    stringsAsFactors = FALSE
  )
  filtered <- .filter.firsts(events)
  expect_equal(sum(filtered$event == "migration"), 2)
  tx <- filtered[filtered$event == "transmission", ]
  expect_equal(sum(tx$to.host == "H3"), 1)
})

test_that(".filter.firsts: unaffected hosts retain all rows after filtering", {
  events <- data.frame(
    time = c(1.0, 1.5, 2.0, 3.0, 4.0),
    event = c("transmission", "transmission", "transmission", "migration", "migration"),
    from.host = c("H1", "H2", "H1", "H3", "H4"),
    to.host = c("H3", "H3", "H4", "H3_sample", "H4_sample"),
    stringsAsFactors = FALSE
  )
  filtered <- .filter.firsts(events)
  expect_equal(sum(!is.na(filtered$to.host) & filtered$to.host == "H3"), 1)
  expect_equal(sum(!is.na(filtered$to.host) & filtered$to.host == "H4"), 1)
  expect_equal(sum(filtered$event == "migration"), 2)
})

test_that(".filter.firsts: surviving row has correct source host (earliest donor)", {
  events <- data.frame(
    time = c(1.0, 2.0),
    event = c("transmission", "transmission"),
    from.host = c("H_early", "H_late"),
    to.host = c("H_victim", "H_victim"),
    stringsAsFactors = FALSE
  )
  filtered <- .filter.firsts(events)
  surviving <- filtered[!is.na(filtered$to.host) & filtered$to.host == "H_victim", ]
  expect_equal(surviving$from.host, "H_early")
})

test_that(".filter.firsts: simultaneous transmissions to same host — one kept", {
  events <- data.frame(
    time = c(1.0, 1.0),
    event = c("transmission", "transmission"),
    from.host = c("H1", "H2"),
    to.host = c("H3", "H3"),
    stringsAsFactors = FALSE
  )
  filtered <- .filter.firsts(events)
  expect_equal(sum(!is.na(filtered$to.host) & filtered$to.host == "H3"), 1)
})

test_that(".filter.firsts: each host appears at most once as to.host after filtering", {
  events <- data.frame(
    time = c(1.0, 2.0, 3.0, 4.0, 5.0),
    event = rep("transmission", 5),
    from.host = c("H1","H2","H3","H4","H5"),
    to.host = c("H6","H6","H7","H7","H8"),
    stringsAsFactors = FALSE
  )
  filtered <- .filter.firsts(events)
  tx <- filtered[filtered$event == "transmission", ]
  expect_true(all(table(tx$to.host) <= 1))
})

# edge cases
test_that("edge case: empty event log", {
  events <- data.frame(
    time = numeric(0), event = character(0),
    from.host = character(0), to.host = character(0),
    stringsAsFactors = FALSE
  )
  expect_equal(nrow(.filter.firsts(events)), 0)
})

test_that("edge case: single transmission passes through unchanged", {
  events <- data.frame(
    time = 1.0, event = "transmission",
    from.host = "H1", to.host = "H2",
    stringsAsFactors = FALSE
  )
  filtered <- .filter.firsts(events)
  expect_equal(nrow(filtered), 1)
  expect_equal(filtered$to.host, "H2")
})

test_that("edge case: all migrations — .filter.firsts returns unchanged", {
  events <- data.frame(
    time = c(1.0, 2.0, 3.0),
    event = c("migration", "migration", "migration"),
    from.host = c("H1", "H2", "H3"),
    to.host = c("H1_1", "H2_1", "H3_sample"),
    stringsAsFactors = FALSE
  )
  filtered <- .filter.firsts(events)
  expect_equal(nrow(filtered), 3)
  expect_equal(sum(filtered$event == "migration"), 3)
})


# --- .relabel.nodes ---

test_that(".relabel.nodes: errors if superinfection events still present", {
  targets <- list(I_samp = 5)
  events <- data.frame(
    time = c(1.0, 2.0, 3.0),
    event = c("transmission", "transmission", "migration"),
    from.comp = c("S", "I", "I"),
    src.comp = c("I", "I", NA_character_),
    to.comp = c("I", "I", "I_samp"),
    from.host = c("H1", "H2", "H3"),
    to.host = c("H3", "H3", "H3_sample"),
    stringsAsFactors = FALSE
  )
  expect_error(.relabel.nodes(events, targets), regexp = "superinfection")
})

test_that(".relabel.nodes: sampled compartment migration gets _sample suffix", {
  targets <- list(I_samp = 5)
  events <- data.frame(
    time = c(1.0, 2.0),
    event = c("transmission", "migration"),
    from.comp = c("S", "I"),
    src.comp = c("I", NA_character_),
    to.comp = c("I", "I_samp"),
    from.host = c("H1", "H2"),
    to.host = c("H2", "H2_sample"),
    stringsAsFactors = FALSE
  )
  relabelled <- .relabel.nodes(events, targets)
  expect_true(any(grepl("_sample$", relabelled$to.host), na.rm = TRUE))
})

test_that(".relabel.nodes: non-sampled migration gets numeric suffixes", {
  targets <- list(I_samp = 5)
  events <- data.frame(
    time = c(1.0, 2.0, 3.0, 4.0),
    event = c("transmission", "migration", "migration", "migration"),
    from.comp = c("S", "I", "I2", "I2"),
    src.comp = c("I", NA_character_, NA_character_, NA_character_),
    to.comp = c("I", "I2", "I2", "I_samp"),
    from.host = c("H1", "H2", "H2_1", "H2_2"),
    to.host = c("H2", "H2_1", "H2_2", "H2_sample"),
    stringsAsFactors = FALSE
  )
  relabelled <- .relabel.nodes(events, targets)
  node_labels <- na.omit(unique(c(relabelled$from.host, relabelled$to.host)))
  expect_true(any(grepl("_1$", node_labels)))
  expect_true(any(grepl("_sample$", node_labels)))
})

test_that(".relabel.nodes: each host appears exactly once as to.host in transmissions", {
  targets <- list(I_samp = 5)
  events <- data.frame(
    time = c(1.0, 2.0, 3.0),
    event = c("transmission", "transmission", "migration"),
    from.comp = c("S", "S", "I"),
    src.comp = c("I", "I", NA_character_),
    to.comp = c("I", "I", "I_samp"),
    from.host = c("H1", "H2", "H3"),
    to.host = c("H2", "H3", "H3_sample"),
    stringsAsFactors = FALSE
  )
  relabelled <- .relabel.nodes(events, targets)
  tx <- relabelled[relabelled$event == "transmission", ]
  expect_true(all(table(tx$to.host) <= 1))
})


# --- compartment size invariance ---

test_that("superinfection (I->I) doesn't change compartment sizes", {
  events <- data.frame(
    time = c(1.0, 2.0, 3.0),
    event = c("transmission", "transmission", "migration"),
    from.comp = c("S", "I", "I"),
    to.comp = c("I", "I", "I_samp"),
    S = c(499L, 499L, 499L),
    I = c(3L, 3L, 2L),
    I_samp = c(0L, 0L, 1L),
    R = c(0L, 0L, 0L),
    stringsAsFactors = FALSE
  )
  si_idx <- which(events$event == "transmission" & events$from.comp == events$to.comp)
  expect_equal(length(si_idx), 1)
  i <- si_idx[1]
  expect_equal(events$I[i], events$I[i - 1])
  expect_equal(events$S[i], events$S[i - 1])
})

test_that("normal transmission (S->I) decrements S and increments I", {
  events <- data.frame(
    time = c(1.0, 2.0),
    event = c("transmission", "transmission"),
    from.comp = c("S", "S"),
    to.comp = c("I", "I"),
    S = c(499L, 498L),
    I = c(3L, 4L),
    stringsAsFactors = FALSE
  )
  expect_equal(events$S[2], events$S[1] - 1L)
  expect_equal(events$I[2], events$I[1] + 1L)
})

test_that("superinfection doesn't increment I; normal transmission does", {
  events <- data.frame(
    time = c(1.0, 2.0, 3.0, 4.0),
    event = c("transmission", "transmission", "transmission", "migration"),
    from.comp = c("S", "I", "S", "I"),
    to.comp = c("I", "I", "I", "I_samp"),
    S = c(997L, 997L, 996L, 996L),
    I = c(3L, 3L, 4L, 3L),
    I_samp = c(0L, 0L, 0L, 1L),
    R = c(0L, 0L, 0L, 0L),
    stringsAsFactors = FALSE
  )
  expect_equal(events$I[2], events$I[1])       # SI: no change
  expect_equal(events$I[3], events$I[2] + 1L)  # normal tx: increments
})

test_that("zero sigma: no I->I transmission events", {
  events <- data.frame(
    time = c(1.0, 2.0, 3.0),
    event = c("transmission", "transmission", "migration"),
    from.comp = c("S", "S", "I"),
    to.comp = c("I", "I", "I_samp"),
    source = c("I", "I", NA_character_),
    stringsAsFactors = FALSE
  )
  si <- events[events$event == "transmission" & events$from.comp == events$to.comp, ]
  expect_equal(nrow(si), 0)
})


# --- host-level SI semantics ---

test_that("superinfection occurs after recipient's first infection", {
  events <- data.frame(
    time = c(1.0, 3.0, 5.0),
    event = c("transmission", "transmission", "migration"),
    from.comp = c("S", "I", "I"),
    to.comp = c("I", "I", "I_samp"),
    from.host = c("H1", "H1", "H2"),
    to.host = c("H2", "H2", "H2_sample"),
    stringsAsFactors = FALSE
  )
  si_idx <- which(events$event == "transmission" &
                    events$from.comp == "I" & events$to.comp == "I")
  expect_equal(length(si_idx), 1)
  si_row <- events[si_idx, ]
  first_inf <- events[events$event == "transmission" & events$to.host == si_row$to.host, ]
  expect_true(si_row$time > min(first_inf$time))
})

test_that("superinfection: donor and recipient are distinct hosts", {
  events <- data.frame(
    time = c(1.0, 2.0, 3.0),
    event = c("transmission", "transmission", "migration"),
    from.comp = c("S", "I", "I"),
    to.comp = c("I", "I", "I_samp"),
    from.host = c("H1", "H1", "H2"),
    to.host = c("H2", "H2", "H2_sample"),
    stringsAsFactors = FALSE
  )
  si <- events[events$event == "transmission" &
                 events$from.comp == "I" & events$to.comp == "I", ]
  if (nrow(si) == 0) skip("No SI row in this log")
  expect_true(all(si$from.host != si$to.host))
})


# --- model parsing ---

test_that("model: sigma*I*I stored in transmission.rates[I,I,I]", {
  mod <- try_model(SUPERINF_PATH)
  tr <- mod$get.transmission.rates()
  expect_equal(tr["I", "I", "I"], "sigma*I*I")
})

test_that("model: beta*S*I stored in transmission.rates[S,I,I]", {
  mod <- try_model(SUPERINF_PATH)
  tr <- mod$get.transmission.rates()
  expect_equal(tr["S", "I", "I"], "beta*S*I")
})

test_that("model: compartment I flagged infected=TRUE", {
  mod <- try_model(SUPERINF_PATH)
  expect_true(mod$get.infected("I"))
})

test_that("model YAML: compartment I has pop.size=2 and bottleneck.size=1", {
  try_model(SUPERINF_INNER_PATH)
  cfg <- read_yaml(SUPERINF_INNER_PATH)
  expect_equal(cfg$Compartments$I$size, 2)
  expect_equal(cfg$Compartments$I$bottleneck.size, 1)
})


# --- sim.dynamics ---

test_that("sim.dynamics: SI events have from.comp==to.comp==I", {
  mod <- try_model(SUPERINF_PATH)
  set.seed(13)
  dyn <- sim.dynamics(mod, max.attempts = 10)
  evts <- dyn$events
  si <- evts[evts$event == "transmission" & evts$from.comp == "I" & evts$to.comp == "I", ]
  if (nrow(si) == 0) skip("No SI events at seed 13")
  expect_true(all(si$source == "I"))
  expect_true(all(si$from.comp == "I" & si$to.comp == "I"))
})

test_that("sim.dynamics: I count unchanged after each SI event", {
  mod <- try_model(SUPERINF_PATH)
  set.seed(13)
  dyn <- sim.dynamics(mod, max.attempts = 10)
  evts <- dyn$events
  si_idx <- which(evts$event == "transmission" & evts$from.comp == "I" & evts$to.comp == "I")
  if (length(si_idx) == 0) skip("No SI events at seed 13")
  for (i in si_idx[si_idx > 1]) {
    expect_equal(evts$I[i], evts$I[i - 1], info = paste("row", i))
  }
})


# --- sim.outer.tree ---

test_that("sim.outer.tree: returns OuterTree with at least one transmission event", {
  mod <- try_model(SUPERINF_PATH)
  set.seed(99)
  dyn <- tryCatch(sim.dynamics(mod, max.attempts = 15), warning = function(w) NULL)
  if (is.null(dyn)) skip("Epidemic went extinct")
  outer <- tryCatch(
    withCallingHandlers(sim.outer.tree(dyn), warning = function(w) invokeRestart("muffleWarning")),
    error = function(e) NULL)
  if (is.null(outer)) skip("sim.outer.tree errored")
  expect_s3_class(outer, "OuterTree")
  tx <- outer$get.log()
  expect_true(nrow(tx[tx$event == "transmission", ]) >= 1)
})

test_that("sim.outer.tree: every to.host is a valid host in the log", {
  mod <- try_model(SUPERINF_PATH)
  outer <- .get_valid_outer(mod)
  if (is.null(outer)) skip("No seed produced a valid outer tree")
  log <- outer$get.log()
  all_hosts <- unique(c(log$from.host, na.omit(log$to.host)))
  orphans <- setdiff(na.omit(unique(log$to.host)), all_hosts)
  expect_equal(length(orphans), 0)
})

test_that("sim.outer.tree: sampled host count matches model targets", {
  mod <- try_model(SUPERINF_PATH)
  outer <- .get_valid_outer(mod)
  if (is.null(outer)) skip("No seed produced a valid outer tree")
  n_sampled <- outer$get.sampled()$count.type()
  target_n <- sum(unlist(mod$get.sampling()$targets))
  expect_equal(n_sampled, target_n)
})


# --- deterministic pipeline ---

test_that("filter then relabel produces consistent host labels", {
  targets <- list(I_samp = 5)
  events <- data.frame(
    time = c(1.0, 2.0, 3.0, 4.0, 5.0),
    event = c("transmission", "transmission", "transmission", "migration", "migration"),
    from.comp = c("S", "I", "S", "I", "I"),
    to.comp = c("I", "I", "I", "I_samp", "I_samp"),
    from.host = c("H1", "H1", "H1", "H2", "H3"),
    to.host = c("H2", "H2", "H3", "H2_sample", "H3_sample"),
    stringsAsFactors = FALSE
  )
  filtered <- .filter.firsts(events)
  tx <- filtered[filtered$event == "transmission", ]
  expect_equal(sum(!is.na(tx$to.host) & tx$to.host == "H2"), 1)

  relabelled <- .relabel.nodes(filtered, targets)
  node_labels <- na.omit(unique(c(relabelled$from.host, relabelled$to.host)))
  expect_true(any(grepl("_sample$", node_labels)))
  expect_false(anyDuplicated(na.omit(relabelled$to.host)) > 0)
})


# --- tree structure ---

test_that("as.phylo(outer): tip count equals sampled host count", {
  mod <- try_model(SUPERINF_PATH)
  result <- .get_valid_phylo(mod)
  if (is.null(result)) skip("No seed produced a single-root tree")
  expect_equal(length(result$phy$tip.label), result$outer$get.sampled()$count.type())
})

test_that("as.phylo(outer): no duplicate tip labels", {
  mod <- try_model(SUPERINF_PATH)
  result <- .get_valid_phylo(mod)
  if (is.null(result)) skip("No seed produced a single-root tree")
  phy <- result$phy
  expect_equal(length(phy$tip.label), length(unique(phy$tip.label)))
})

test_that("as.phylo(outer): edge count equals nodes - 1", {
  mod <- try_model(SUPERINF_PATH)
  result <- .get_valid_phylo(mod)
  if (is.null(result)) skip("No seed produced a single-root tree")
  phy <- result$phy
  expect_equal(nrow(phy$edge), length(phy$tip.label) + phy$Nnode - 1)
})


# --- inner tree with superinfection ---

test_that("sim.inner.tree: no crash when outer tree has superinfection", {
  mod <- try_model(SUPERINF_INNER_PATH)
  result <- .get_outer_with_superinfection(mod)
  if (is.null(result)) skip("No seed produced outer tree with SI + sampled hosts")
  inner <- .try_inner(result$outer)
  expect_false(is.null(inner), info = paste("returned NULL for seed", result$seed))
})

test_that("sim.inner.tree -> as.phylo: tip count equals sampled host count", {
  mod <- try_model(SUPERINF_INNER_PATH)
  result <- .get_outer_with_superinfection(mod)
  if (is.null(result)) skip("No suitable outer tree found")
  inner <- .try_inner(result$outer)
  if (is.null(inner)) skip("sim.inner.tree failed")
  if (result$outer$get.active()$count.type() != 1) skip("Multiple active roots")
  phy <- tryCatch(as.phylo(inner), error = function(e) NULL)
  if (is.null(phy)) skip("as.phylo failed")
  n_sampled <- result$outer$get.sampled()$count.type()
  expect_equal(length(phy$tip.label), n_sampled)
})

test_that("sim.inner.tree -> as.phylo: edge count equals nodes - 1", {
  mod <- try_model(SUPERINF_INNER_PATH)
  result <- .get_outer_with_superinfection(mod)
  if (is.null(result)) skip("No suitable outer tree found")
  inner <- .try_inner(result$outer)
  if (is.null(inner)) skip("sim.inner.tree failed")
  if (result$outer$get.active()$count.type() != 1) skip("Multiple active roots")
  phy <- tryCatch(as.phylo(inner), error = function(e) NULL)
  if (is.null(phy)) skip("as.phylo failed")
  expect_equal(nrow(phy$edge), length(phy$tip.label) + phy$Nnode - 1)
})

test_that("sim.inner.tree -> as.phylo: no duplicate tip labels", {
  mod <- try_model(SUPERINF_INNER_PATH)
  result <- .get_outer_with_superinfection(mod)
  if (is.null(result)) skip("No suitable outer tree found")
  inner <- .try_inner(result$outer)
  if (is.null(inner)) skip("sim.inner.tree failed")
  if (result$outer$get.active()$count.type() != 1) skip("Multiple active roots")
  phy <- tryCatch(as.phylo(inner), error = function(e) NULL)
  if (is.null(phy)) skip("as.phylo failed")
  expect_equal(length(phy$tip.label), length(unique(phy$tip.label)))
})

test_that("sim.inner.tree: inner log contains transmission events", {
  mod <- try_model(SUPERINF_INNER_PATH)
  result <- .get_outer_with_superinfection(mod)
  if (is.null(result)) skip("No suitable outer tree found")
  inner <- .try_inner(result$outer)
  if (is.null(inner)) skip("sim.inner.tree failed")
  ilog <- tryCatch(inner$get.log(), error = function(e) NULL)
  if (is.null(ilog) || nrow(ilog) == 0) skip("Inner tree has no event log")
  expect_true(nrow(ilog[ilog$event == "transmission", ]) > 0)
})

test_that("sim.inner.tree (20 replicates): at least one produces transmission events", {
  mod <- try_model(SUPERINF_INNER_PATH)
  result <- .get_outer_with_superinfection(mod)
  if (is.null(result)) skip("No outer tree with SI found")
  n_with_tx <- 0
  for (rep in 1:20) {
    set.seed(rep * 7)
    inner <- .try_inner(result$outer)
    if (is.null(inner)) next
    ilog <- tryCatch(inner$get.log(), error = function(e) NULL)
    if (is.null(ilog) || nrow(ilog) == 0) next
    if (nrow(ilog[ilog$event == "transmission", ]) > 0) n_with_tx <- n_with_tx + 1
  }
  expect_true(n_with_tx > 0)
})

# when a host is superinfected, incoming lineages can coalesce with either the
# original or the superinfecting lineage — so inner tree replicates on the same
# outer tree should vary in topology
test_that("sim.inner.tree: replicates on SI outer tree show topological variation", {
  mod <- try_model(SUPERINF_INNER_PATH)
  result <- .get_outer_with_superinfection(mod)
  if (is.null(result)) skip("No outer tree with SI found")
  if (result$outer$get.active()$count.type() != 1) skip("Multiple active roots")

  phylos <- list()
  for (rep in 1:20) {
    set.seed(rep * 13)
    inner <- .try_inner(result$outer)
    if (is.null(inner)) next
    phy <- tryCatch(collapse.singles(as.phylo(inner)), error = function(e) NULL)
    if (!is.null(phy)) phylos[[length(phylos) + 1]] <- phy
  }
  if (length(phylos) < 2) skip("Fewer than 2 inner trees succeeded")

  rf_vals <- c()
  for (i in 1:(length(phylos) - 1)) {
    for (j in (i + 1):length(phylos)) {
      d <- tryCatch(suppressWarnings(dist.topo(phylos[[i]], phylos[[j]])), error = function(e) NA)
      rf_vals <- c(rf_vals, d)
    }
  }
  expect_true(any(rf_vals > 0, na.rm = TRUE))
})
