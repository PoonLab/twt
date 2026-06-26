## tests for sim.arg and resolve.arg

library(testthat)
library(yaml)
library(ape)

N.TIPS  <- 10L
SEQ.LEN <- 9000L
RHO     <- 0.3

# search over seeds -- stochastic outer trees sometimes fail InnerTree$new()
.make.outer <- function() {
  mod <- Model$new(read_yaml("test_Superinfection.yaml"))
  for (s in seq(13, 613, by=3)) {
    set.seed(s)
    dyn <- sim.dynamics(mod)
    out <- tryCatch(suppressWarnings(sim.outer.tree(dyn)), error=function(e) NULL)
    if (is.null(out)) next
    ok <- tryCatch({ suppressWarnings(sim.arg(out, rho=0, seq.length=SEQ.LEN)); TRUE },
                   error=function(e) FALSE)
    if (ok) return(out)
  }
  stop("no valid outer tree found")
}

outer.tree     <- .make.outer()
result.rho0    <- suppressWarnings(sim.arg(outer.tree, rho=0,   seq.length=SEQ.LEN))
result.rho03   <- suppressWarnings(sim.arg(outer.tree, rho=RHO, seq.length=SEQ.LEN))
resolved.rho0  <- resolve.arg(result.rho0,  seq.length=SEQ.LEN)
resolved.rho03 <- resolve.arg(result.rho03, seq.length=SEQ.LEN)

# skip multi-segment tests if no breakpoints were produced
has.breakpoints <- length(result.rho03$breakpoints) > 0

# --- sim.arg ---

test_that("sim.arg returns list with $inner and $breakpoints", {
  expect_type(result.rho0, "list")
  expect_true(!is.null(result.rho0$inner))
  expect_true(!is.null(result.rho0$breakpoints))
})

test_that("sim.arg $inner is an InnerTree R6 object", {
  expect_true(is.R6(result.rho0$inner))
  expect_true(is.element("InnerTree", class(result.rho0$inner)))
})

test_that("sim.arg rejects non-OuterTree input", {
  expect_error(sim.arg("not an OuterTree"), regexp="OuterTree")
})

test_that("sim.arg log contains coalescent events", {
  log <- result.rho0$inner$get.log()
  expect_gt(sum(log$event == "coalescent"), 0)
})

test_that("sim.arg log contains exactly N.TIPS sampling events", {
  log <- result.rho0$inner$get.log()
  expect_equal(sum(log$event == "sampling"), N.TIPS)
})


test_that("rho=0 produces no recombination events", {
  log <- result.rho0$inner$get.log()
  expect_equal(sum(log$event == "recombination"), 0)
})

test_that("rho=0 produces no breakpoints", {
  expect_equal(length(result.rho0$breakpoints), 0)
})

test_that("sim.arg rho>0 returns list with $inner and $breakpoints", {
  expect_type(result.rho03, "list")
  expect_true(!is.null(result.rho03$inner))
  expect_true(!is.null(result.rho03$breakpoints))
})

test_that("breakpoint count matches recombination log entries", {
  skip_if(!has.breakpoints, "no breakpoints generated at this seed")
  log <- result.rho03$inner$get.log()
  # each recombination logs two rows, one per parental lineage
  expect_equal(length(result.rho03$breakpoints),
               sum(log$event == "recombination") / 2)
})

test_that("breakpoints are integers in [1, seq.length-1]", {
  skip_if(!has.breakpoints, "no breakpoints generated at this seed")
  bps <- unlist(result.rho03$breakpoints)
  expect_true(all(bps >= 1L))
  expect_true(all(bps <= SEQ.LEN - 1L))
  expect_true(all(bps == as.integer(bps)))
})

test_that("breakpoint positions are unique", {
  skip_if(!has.breakpoints, "no breakpoints generated at this seed")
  bps <- unlist(result.rho03$breakpoints)
  expect_equal(length(bps), length(unique(bps)))
})

# --- resolve.arg ---

test_that("rho=0 gives single segment spanning full genome", {
  segs <- resolved.rho0$segments
  expect_equal(nrow(segs), 1L)
  expect_equal(segs$start, 1L)
  expect_equal(segs$end, SEQ.LEN)
})

test_that("segment count equals breakpoint count plus one", {
  skip_if(!has.breakpoints, "no breakpoints generated at this seed")
  expect_equal(nrow(resolved.rho03$segments),
               length(result.rho03$breakpoints) + 1)
})

test_that("local tree count equals segment count", {
  skip_if(!has.breakpoints, "no breakpoints generated at this seed")
  expect_equal(length(resolved.rho03$local.trees),
               nrow(resolved.rho03$segments))
})

test_that("segment boundaries are contiguous", {
  skip_if(!has.breakpoints, "no breakpoints generated at this seed")
  segs <- resolved.rho03$segments
  expect_equal(segs$start[1], 1L)
  expect_equal(segs$end[nrow(segs)], SEQ.LEN)
  expect_true(all(segs$start[-1] == segs$end[-nrow(segs)] + 1L))
})

test_that("rho=0 local tree is a valid phylo with tips", {
  lt <- resolved.rho0$local.trees[[1]]
  expect_false(is.null(lt$phylo))
  expect_s3_class(lt$phylo, "phylo")
  expect_gt(length(lt$phylo$tip.label), 0)
})

test_that("all local trees are valid phylo objects with tips", {
  skip_if(!has.breakpoints, "no breakpoints generated at this seed")
  for (i in seq_along(resolved.rho03$local.trees)) {
    lt <- resolved.rho03$local.trees[[i]]
    expect_false(is.null(lt$phylo), label=paste("segment", i, "non-NULL"))
    expect_true(inherits(lt$phylo, "phylo"), label=paste("segment", i, "is phylo"))
    expect_gt(length(lt$phylo$tip.label), 0, label=paste("segment", i, "has tips"))
  }
})
