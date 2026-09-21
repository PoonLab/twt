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
test_that("Host lineage pool: basic activate/deactivate/query", {
  h <- Host$new(name="H1", compartment="I")
  expect_false(h$is.pool.initialized())
  h$init.pool(5)
  expect_true(h$is.pool.initialized())
  expect_equal(h$get.pool.size(), 5)
  expect_equal(h$count.inactive.slots(), 5)
  expect_equal(h$count.active.slots(), 0)

  p <- Pathogen$new(name="P1", end.time=1)
  new.id <- h$activate.new.slot(p)
  expect_equal(p$get.slot.id(), new.id)
  expect_true(h$is.slot.active(new.id))
  expect_equal(h$get.slot.occupant(new.id)$get.name(), "P1")
  expect_equal(h$count.active.slots(), 1)
  expect_equal(h$count.inactive.slots(), 4)

  h$deactivate.slot(new.id)
  expect_false(h$is.slot.active(new.id))
  expect_equal(h$count.active.slots(), 0)
  # pool size must not shrink after deactivation
  expect_equal(h$get.pool.size(), 5)
})

test_that("Host lineage pool: init.pool is idempotent", {
  h <- Host$new(name="H1", compartment="I")
  h$init.pool(10)
  h$init.pool(999)  # should be ignored, pool already initialized
  expect_equal(h$get.pool.size(), 10)
})

test_that("sim.arg does not exceed pop.size at high rho (fixed lineage pool)", {
  settings <- read_yaml("test_Superinfection.yaml")
  settings$Parameters$sigma <- 0.05
  mod <- Model$new(settings)
  set.seed(33)
  dyn <- tryCatch(sim.dynamics(mod, max.attempts=10), error=function(e) NULL)
  if (is.null(dyn)) skip("could not build dynamics for this seed")
  outer <- tryCatch(
    withCallingHandlers(sim.outer.tree(dyn), warning=function(w) invokeRestart("muffleWarning")),
    error=function(e) NULL)
  if (is.null(outer)) skip("could not build outer tree for this seed")

  # rho=5,7,10 previously exceeded pop.size and crashed/errored before
  # the fixed lineage pool (Art's design: sample recombination parents
  # from a fixed pool of p.size lineages instead of always creating a
  # new lineage de novo). Now they should all succeed.
  for (rho in c(5, 7, 10)) {
    arg <- tryCatch(sim.arg(outer, rho=rho, seq.length=9000), error=function(e) e)
    expect_false(inherits(arg, "error"),
                 info=paste("rho =", rho, "should not error with fixed lineage pool"))
  }
})

test_that("resolve.arg produces valid trees on top of the fixed lineage pool", {
  settings <- read_yaml("test_Superinfection.yaml")
  settings$Parameters$sigma <- 0.05
  mod <- Model$new(settings)
  set.seed(33)
  dyn <- tryCatch(sim.dynamics(mod, max.attempts=10), error=function(e) NULL)
  if (is.null(dyn)) skip("could not build dynamics for this seed")
  outer <- tryCatch(
    withCallingHandlers(sim.outer.tree(dyn), warning=function(w) invokeRestart("muffleWarning")),
    error=function(e) NULL)
  if (is.null(outer)) skip("could not build outer tree for this seed")

  n.sampled <- outer$get.sampled()$count.type()
  arg <- sim.arg(outer, rho=2, seq.length=9000)
  res <- resolve.arg(arg, seq.length=9000)

  valid <- sapply(res$local.trees, function(lt) {
    phy <- lt$phylo
    inherits(phy, "phylo") &&
      length(phy$tip.label) == n.sampled &&
      sum(duplicated(phy$tip.label)) == 0 &&
      !any(is.na(phy$edge.length)) &&
      !any(phy$edge.length < 0, na.rm=TRUE)
  })
  expect_true(all(valid))
})
test_that("resolve.arg: per-child breakpoint lookup is correct with multiple distinct breakpoints", {
  # A and C recombine independently at different breakpoints (300, 700)
  times <- c(A=10,B=10,C=10,D=10, L1=8,R1=8, L2=7,R2=7,
             AB=5,CD=5, RR=4, ABCD=2, ROOT=1)
  paths <- lapply(names(times), function(n) Pathogen$new(name=n, end.time=times[[n]]))
  names(paths) <- names(times)
  log <- data.frame(
    time      = c(10,10,10,10, 8,8, 7,7, 5,5, 5,5, 4,4, 2,2, 1,1),
    event     = c(rep("sampling",4), rep("recombination",2), rep("recombination",2),
                  rep("coalescent",2), rep("coalescent",2), rep("coalescent",2),
                  rep("coalescent",2), rep("coalescent",2)),
    pathogen1 = c("A","B","C","D", "A","A", "C","C", "AB","AB", "CD","CD",
                  "RR","RR", "ABCD","ABCD", "ROOT","ROOT"),
    pathogen2 = c(NA,NA,NA,NA, "L1","R1", "L2","R2", "L1","B", "L2","D",
                  "R1","R2", "AB","CD", "ABCD","RR"),
    stringsAsFactors=FALSE
  )
  fake.inner <- list(get.log=function() log, get.all.pathogens=function() paths)
  arg.result <- list(inner=fake.inner, breakpoints=list(A=300L, C=700L))
  res <- resolve.arg(arg.result, seq.length=1000L)

  expect_equal(length(res$local.trees), 3)

  phy1 <- collapse.singles(res$local.trees[[1]]$phylo)
  phy2 <- collapse.singles(res$local.trees[[2]]$phylo)
  phy3 <- collapse.singles(res$local.trees[[3]]$phylo)

  expect_true(is.monophyletic(phy1, c("C","D")))
  expect_true(is.monophyletic(phy2, c("C","D")))
  expect_false(is.monophyletic(phy2, c("A","C")))
  expect_true(is.monophyletic(phy3, c("A","C")))
  expect_false(is.monophyletic(phy3, c("C","D")))
})


# --- skip.recomb.free ---

test_that(".compute.recomb.free.hosts flags hosts correctly under bottleneck/SI/tip-count conditions", {
  # H1: bottleneck=1, no SI, 1 sampled tip -> recomb.free TRUE
  # H2: bottleneck=1, no SI, 2 sampled tips -> recomb.free FALSE (tip count)
  # H3: bottleneck=1, 1 SI event, 1 sampled tip -> recomb.free FALSE (superinfected)
  # H4: bottleneck=2, no SI, 1 sampled tip -> recomb.free FALSE (loose bottleneck)
  events <- data.frame(
    time      = c(5,5,5,5, 3, 1,1,1,1,1),
    event     = c(rep("transmission", 4), "transmission",
                  rep("migration", 5)),
    from.comp = c("S","S","S","S", "I",
                  "I","I","I","I","I"),
    to.comp   = c("I","I","I","J", "I",
                  "T","T","T","T","T"),
    from.host = c("Src1","Src2","Src3","Src4", "OtherHost",
                  "H1","H2","H2","H3","H4"),
    to.host   = c("H1","H2","H3","H4", "H3",
                  NA,NA,NA,NA,NA),
    stringsAsFactors = FALSE
  )

  fake.mod <- list(
    get.infected        = function() c(S = FALSE, I = TRUE, J = FALSE),
    get.bottleneck.size = function(comp) if (comp == "J") "2" else "1"
  )
  fake.inner <- list(has.target = function(comp) comp == "T")

  profiles <- .compute.recomb.free.hosts(events, fake.mod, fake.inner, envir = new.env())

  expect_true(profiles$H1$recomb.free)
  expect_equal(profiles$H1$bottleneck.size, 1)

  expect_false(profiles$H2$recomb.free)  # two sampled tips
  expect_false(profiles$H3$recomb.free)  # superinfected
  expect_false(profiles$H4$recomb.free)  # bottleneck size 2
})

test_that(".draw.next.event skips recombination for hosts flagged recomb.free", {
  fake.host <- list(
    count.pathogens = function() 1L,
    get.compartment = function() "I",
    get.name        = function() "H1",
    get.pathogens   = function() list("P1")
  )
  fake.active <- list(
    get.hosts        = function() list(fake.host),
    get.host.by.name = function(name) fake.host
  )
  fake.mod <- list(get.coalescent.rate = function(comp) "0")

  set.seed(1)
  ev.free <- .draw.next.event(fake.active, fake.mod, rho = 0.5, envir = new.env(),
                               host.profiles = list(H1 = list(recomb.free = TRUE,
                                                               bottleneck.size = 1)))
  expect_null(ev.free)

  set.seed(1)
  ev.default <- .draw.next.event(fake.active, fake.mod, rho = 0.5, envir = new.env())
  expect_false(is.null(ev.default))
  expect_equal(ev.default$type, "recombination")
  expect_equal(ev.default$host, "H1")
})

test_that("sim.arg with skip.recomb.free=TRUE runs and produces valid trees", {
  result <- suppressWarnings(
    sim.arg(outer.tree, rho = RHO, seq.length = SEQ.LEN, skip.recomb.free = TRUE)
  )
  expect_type(result, "list")
  expect_true(is.R6(result$inner))

  resolved <- resolve.arg(result, seq.length = SEQ.LEN)
  valid <- sapply(resolved$local.trees, function(lt) {
    inherits(lt$phylo, "phylo") &&
      length(lt$phylo$tip.label) == N.TIPS &&
      sum(duplicated(lt$phylo$tip.label)) == 0
  })
  expect_true(all(valid))
})

test_that("skip.recomb.free eliminates recombination events for hosts flagged recomb-free", {
  inner.for.profile <- InnerTree$new(outer.tree)
  mod <- inner.for.profile$get.model()

  events <- outer.tree$get.log()
  events$time <- as.numeric(events$time)
  events <- events[order(events$time, decreasing = TRUE), ]

  profiles <- .compute.recomb.free.hosts(events, mod, inner.for.profile, envir = new.env())
  free.hosts <- names(profiles)[sapply(profiles, function(p) isTRUE(p$recomb.free))]
  skip_if(length(free.hosts) == 0, "no recomb-free host at this seed")

  result.skip <- suppressWarnings(
    sim.arg(outer.tree, rho = RHO, seq.length = SEQ.LEN, skip.recomb.free = TRUE)
  )
  log <- result.skip$inner$get.log()
  recomb.hosts <- unique(log$from.host[log$event == "recombination"])
  expect_false(any(free.hosts %in% recomb.hosts))
})


# --- pool-overflow bug: Host$activate.new.slot has no bound against pool.size ---

test_that("Host pool: activate.new.slot has no bound and can exceed pool.size", {
  # raw primitive has no guard -- the check lives at the
  # .do.recombination call site, tested next
  h <- Host$new(name = "H1", compartment = "I")
  h$init.pool(2)

  make.fake.pathogen <- function() {
    slot <- NA_integer_
    list(
      set.slot.id = function(id) slot <<- id,
      get.slot.id = function() slot
    )
  }
  p1 <- make.fake.pathogen()
  p2 <- make.fake.pathogen()
  p3 <- make.fake.pathogen()

  h$activate.new.slot(p1)
  h$activate.new.slot(p2)
  expect_equal(h$count.active.slots(), 2)
  expect_equal(h$get.pool.size(), 2)  # pool is now at capacity

  h$activate.new.slot(p3)
  expect_gt(h$count.active.slots(), h$get.pool.size())  # invariant violated
})

test_that(".do.recombination errors instead of silently exceeding pool.size when every slot is already claimed", {
  # pool's already full (2/2), so a fresh lineage needing a slot
  # should error instead of pushing count past pool.size
  fake.host <- list(
    is.pool.initialized = function() TRUE,
    count.active.slots   = function() 2L,
    get.pool.size        = function() 2L,
    get.name              = function() "H1",
    activate.new.slot     = function(pathogen) stop("should not be reached")
  )
  fake.active <- list(get.host.by.name = function(name) fake.host)
  fake.inner  <- list(get.active = function() fake.active)
  fake.pathogen <- list(get.slot.id = function() NA_integer_)

  expect_error(
    .do.recombination("H1", fake.pathogen, fake.inner, time = 0.5,
                       seq.length = 9000L, p.size = 2L),
    regexp = "lineage pool is full"
  )
})
