# Tests for the transmission bottleneck with several active lineages
# (issue 185). Mock objects keep these independent of the real simulator.

MockPathogen <- R6::R6Class("MockPathogen",
  public = list(
    name = NULL, slot = NA_integer_, parent = NULL, children = list(),
    start = NA_real_,
    initialize = function(name) self$name <- name,
    get.name = function() self$name,
    set.start.time = function(t) self$start <- t,
    get.slot.id = function() self$slot,
    set.slot.id = function(id) self$slot <- id,
    add.child = function(p) self$children <- c(self$children, list(p)),
    set.parent = function(p) self$parent <- p
  ))

MockHost <- R6::R6Class("MockHost",
  public = list(
    paths = list(),
    get.name = function() "H1",
    get.compartment = function() "I",
    count.pathogens = function() length(self$paths),
    add.pathogen = function(p) {
      self$paths <- c(self$paths, list(p))
      invisible(NULL)
    },
    remove.pathogen = function(i) {
      p <- self$paths[[i]]
      self$paths <- self$paths[-i]
      p
    },
    sample.pathogen = function(remove = FALSE) {
      i <- if (length(self$paths) == 1) 1L else sample.int(length(self$paths), 1)
      if (remove) self$remove.pathogen(i) else self$paths[[i]]
    },
    activate.slot = function(id, p) invisible(NULL),
    deactivate.slot = function(id) invisible(NULL)
  ))

MockInner <- R6::R6Class("MockInner",
  public = list(
    host = NULL, n = 0L, events = 0L,
    initialize = function(host) self$host <- host,
    get.active = function() {
      h <- self$host
      list(get.host.by.name = function(name) h)
    },
    new.pathogen = function(time) {
      self$n <- self$n + 1L
      MockPathogen$new(paste0("anc", self$n))
    },
    add.event = function(e) {
      self$events <- self$events + 1L
      invisible(NULL)
    }
  ))

make.case <- function(k) {
  host <- MockHost$new()
  for (i in seq_len(k)) host$add.pathogen(MockPathogen$new(paste0("p", i)))
  list(host = host, inner = MockInner$new(host), original = host$paths)
}

# does lineage p end up in the host, or descend into a lineage that does?
has.ancestor.in.host <- function(p, host) {
  repeat {
    if (any(vapply(host$paths, identical, logical(1), p))) return(TRUE)
    if (is.null(p$parent)) return(FALSE)
    p <- p$parent
  }
}

all.have.ancestor <- function(case) {
  all(vapply(case$original, has.ancestor.in.host, logical(1),
             host = case$host))
}


test_that(".do.bottleneck gives every active lineage an ancestor in the host", {
  set.seed(1)
  bad <- character(0)  # describes any case that breaks the rule
  for (size in 2:4) {
    for (k in 1:6) {
      for (rep in 1:20) {
        case <- make.case(k)
        twt:::.do.bottleneck(case$host, case$inner, 1.0, size)
        n <- case$host$count.pathogens()
        if (n < 1 || n > min(k, size) || !all.have.ancestor(case)) {
          bad <- c(bad, sprintf("bottleneck %d, %d lineages, repeat %d",
                                size, k, rep))
        }
      }
    }
  }
  expect_equal(bad, character(0))
})


test_that(".do.bottleneck: mean number of distinct ancestors is b(1-(1-1/b)^k)", {
  set.seed(42)
  b <- 3
  k <- 4
  n.anc <- replicate(4000, {
    case <- make.case(k)
    twt:::.do.bottleneck(case$host, case$inner, 1.0, b)
    case$host$count.pathogens()
  })
  expect_equal(mean(n.anc), b * (1 - (1 - 1/b)^k), tolerance = 0.03)
})


test_that(".do.coalescent with a bottleneck of one merges everything into one lineage", {
  case <- make.case(5)
  twt:::.do.coalescent("H1", case$inner, 1.0, bottleneck = TRUE,
                       preset.bottleneck.size = 1)
  expect_equal(case$host$count.pathogens(), 1)
  expect_true(all.have.ancestor(case))
})


test_that(".do.coalescent hands bottlenecks of two or more to .do.bottleneck", {
  set.seed(3)
  counts <- integer(100)
  bad <- integer(0)  # repeats where a lineage lost its ancestor
  for (rep in 1:100) {
    case <- make.case(5)
    twt:::.do.coalescent("H1", case$inner, 1.0, bottleneck = TRUE,
                         preset.bottleneck.size = 3)
    counts[rep] <- case$host$count.pathogens()
    if (!all.have.ancestor(case)) bad <- c(bad, rep)
  }
  expect_equal(bad, integer(0))
  expect_true(all(counts <= 3))
  # the old pairwise loop always left exactly 3; independent picks sometimes
  # leave fewer because two lineages share a transmitted lineage
  expect_true(any(counts < 3))
})
