## tests for sim.arg and resolve.arg

N.TIPS  <- 10L
SEQ.LEN <- 9000L
RHO     <- 0.3

# seed 31 verified to produce a valid outer tree for sim.arg
.make.outer <- function() {
  mod <- Model$new(read_yaml("test_Superinfection.yaml"))
  set.seed(31)
  dyn <- sim.dynamics(mod)
  suppressWarnings(sim.outer.tree(dyn))
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


test_that("recombinant pathogens have exactly two parents", {
  skip_if(!has.breakpoints, "no breakpoints generated at this seed")
  log <- result.rho03$inner$get.log()
  recomb.children <- unique(log$pathogen1[log$event == "recombination"])
  all.paths <- result.rho03$inner$get.all.pathogens()
  n.parents <- sapply(recomb.children, function(n) {
    p <- all.paths[[n]]
    if (is.null(p)) return(NA)
    length(p$get.parents())
  })
  expect_true(all(n.parents == 2, na.rm=TRUE))
})

test_that("non-recombinant pathogens have at most one parent", {
  log <- result.rho03$inner$get.log()
  recomb.children <- unique(log$pathogen1[log$event == "recombination"])
  all.paths <- result.rho03$inner$get.all.pathogens()
  other.names <- setdiff(names(all.paths), recomb.children)
  n.parents <- sapply(other.names, function(n) length(all.paths[[n]]$get.parents()))
  expect_true(all(n.parents <= 1))
})

test_that("each breakpoint is named by its recombinant child pathogen", {
  skip_if(!has.breakpoints, "no breakpoints generated at this seed")
  log <- result.rho03$inner$get.log()
  recomb.children <- unique(log$pathogen1[log$event == "recombination"])
  expect_true(all(names(result.rho03$breakpoints) %in% recomb.children))
})

test_that("every recombination event has a corresponding breakpoint", {
  skip_if(!has.breakpoints, "no breakpoints generated at this seed")
  log <- result.rho03$inner$get.log()
  recomb.children <- unique(log$pathogen1[log$event == "recombination"])
  expect_true(all(recomb.children %in% names(result.rho03$breakpoints)))
})

test_that("recombinant pathogens have breakpoint position stored", {
  skip_if(!has.breakpoints, "no breakpoints generated at this seed")
  log <- result.rho03$inner$get.log()
  recomb.children <- unique(log$pathogen1[log$event == "recombination"])
  all.paths <- result.rho03$inner$get.all.pathogens()
  bps <- sapply(recomb.children, function(n) {
    p <- all.paths[[n]]
    if (is.null(p)) return(NA)
    p$get.breakpoint()
  })
  expect_true(all(!is.na(bps)))
  expect_true(all(bps >= 1L & bps <= SEQ.LEN - 1L, na.rm=TRUE))
})

test_that("non-recombinant pathogens have NA breakpoint", {
  log <- result.rho03$inner$get.log()
  recomb.children <- unique(log$pathogen1[log$event == "recombination"])
  all.paths <- result.rho03$inner$get.all.pathogens()
  other.names <- setdiff(names(all.paths), recomb.children)
  bps <- sapply(other.names, function(n) all.paths[[n]]$get.breakpoint())
  expect_true(all(is.na(bps)))
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
  segs <- resolved.rho03$local.trees
  valid <- sapply(segs, function(lt) !is.null(lt$phylo) && inherits(lt$phylo, "phylo") && length(lt$phylo$tip.label) > 0)
  expect_true(all(valid))
})
test_that("resolve.arg produces genuinely divergent topology (hand-built positive control)", {
  times <- c(A=10, B=10, C=10, L=8, R=8, BL=5, CR=5, ROOT=1)
  paths <- lapply(names(times), function(n) Pathogen$new(name=n, end.time=times[[n]]))
  names(paths) <- names(times)
  log <- data.frame(
    time      = c(10,10,10,  8,8,  5,5,  5,5,  1,1),
    event     = c(rep("sampling", 3),
                  rep("recombination", 2),
                  rep("coalescent", 2),
                  rep("coalescent", 2),
                  rep("coalescent", 2)),
    pathogen1 = c("A","B","C", "A","A", "BL","BL", "CR","CR", "ROOT","ROOT"),
    pathogen2 = c(NA,NA,NA, "L","R", "L","B", "R","C", "BL","CR"),
    stringsAsFactors=FALSE
  )
  fake.inner <- list(get.log=function() log, get.all.pathogens=function() paths)
  arg.result <- list(inner=fake.inner, breakpoints=list(A=500L))
  res <- resolve.arg(arg.result, seq.length=1000L)

  phy1 <- res$local.trees[[1]]$phylo
  phy2 <- res$local.trees[[2]]$phylo
  expect_true(is.monophyletic(phy1, c("A","B")))
  expect_false(is.monophyletic(phy1, c("A","C")))
  expect_true(is.monophyletic(phy2, c("A","C")))
  expect_false(is.monophyletic(phy2, c("A","B")))
})
