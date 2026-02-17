library(twt)

test_that("Initialize OuterTree object", {
  settings <- yaml.load_file("test_SIR.yaml")
  mod <- Model$new(settings)
  outer <- OuterTree$new(mod)
  
  result <- outer$get.targets()
  expected <- list('I_samp'=10)
  expect_equal(result, expected)
  
  outer.log <- outer$get.log()
  expect_true(is.data.frame(outer.log))
  expect_true(nrow(outer.log)==0)
  expect_equal(names(outer.log), c('time', 'event', 'from.comp', 'to.comp', 
                                   'from.host', 'to.host'))
  
  event <- list(time=1.0, event='transmission', from.comp='S', to.comp='I', 
                from.host='S', to.host='I')
  outer$add.event(event)
  expect_true(outer$get.nrow()==1)
})

test_that("Convert OuterTree to phylo", {
  # modify settings for smaller event log
  settings <- yaml.load_file("test_SIR.yaml")
  settings$Parameters$beta <- 0.02
  settings$Parameters$psi <- 2.0
  settings$Sampling$targets$I_samp <- 3
  mod <- Model$new(settings)
  
  # FIXME: This test is fragile - if the following simulation code is changed,
  #        then the expected values may be different.
  set.seed(13)
  event.log <- sim.dynamics(mod)
  
  outer <- sim.outer.tree(mod, event.log)
  phy <- as.phylo(outer)
  expect_equal(Ntip(phy), 3)
  expect_true(is.rooted(phy))
  expect_false(is.binary(phy))
  expect_true(has.singles(phy))
  
  phy2 <- collapse.singles(phy)
  expect.phy <- read.tree(
    text="((I_2_sample:0.05211251,I_1_sample:0.07238021):0.1641128,I_4_sample:0.2075166);")
  result <- comparePhylo(phy2, expect.phy)
  expect_true(is.element("Both trees have the same tip labels.", 
              result$messages))
  
  result <- sort(phy2$edge.length)
  expected <- sort(expect.phy$edge.length)
  expect_equal(result, expected, tolerance=1e-6)
  
  m1 <- cophenetic(phy2)
  m2 <- cophenetic(expect.phy)
  labels <- c("I_1_sample", "I_2_sample", "I_4_sample")
  result <- m1[labels, labels]
  expected <- m2[labels, labels]
  expect_equal(result, expected, tolerance=1e-6)
})


test_that("Reordering an eventlog", {
  events <- data.frame(
    time=c(0.4, 0.1, 0.5, 0.2, 0.7),
    from.host=c('C', 'A', 'B', 'A', 'D'),
    to.host=c('D', 'B', 'B_samp', 'C', 'D_samp')
    )
  result <- .reorder.events(events, 'A', order='postorder', decreasing=TRUE)
  expected <- c('D_samp', 'D', 'C', 'B_samp', 'B', 'A')
  expect_equal(result, expected)
  
  result <- .reorder.events(events, 'A', order='postorder', decreasing=FALSE)
  expected <- c('B_samp', 'B', 'D_samp', 'D', 'C', 'A')
  expect_equal(result, expected)
  
  result <- .reorder.events(events, 'A', order='preorder', decreasing=TRUE)
  expected <- c('A', 'C', 'D', 'D_samp', 'B', 'B_samp')
  expect_equal(result, expected)
  
  result <- .reorder.events(events, 'A', order='preorder', decreasing=FALSE)
  expected <- c('A', 'B', 'B_samp', 'C', 'D', 'D_samp')
  expect_equal(result, expected)
})


test_that("Filter out superinfection events", {
  events <- data.frame(
    time=c(0.4, 0.1, 0.5, 0.2, 0.7),
    event=rep('transmission', 5),
    from.host=c('E', 'A', 'B', 'A', 'B'),
    to.host=c('D', 'B', 'E', 'E', 'E')
  )
  result <- .filter.firsts(events)
  expected <- data.frame(
    time=c(0.4, 0.1, 0.2),
    event=rep('transmission', 3),
    from.host=c('E', 'A', 'A'),
    to.host=c('D', 'B', 'E'),
    row.names=as.integer(c(1, 2, 4))
  )
  expect_equal(result, expected)
})


test_that("Migration events relabel nodes", {
  events <- data.frame(
    time=seq(0.1, 0.4, 0.1),
    event=c("transmission", "migration", "transmission", "migration"),
    from.comp=c("S", "I1", "S", "I2"),
    to.comp=c("I1", "I2", "I1", "I2_samp"),
    from.host=c("A", "B", "B", "B"),
    to.host=c("B", NA, "C", NA)
  )
  targets <- list("I1_samp"=1, "I2_samp"=1)
  result <- .relabel.nodes(events, targets)
  expected <- data.frame(
    time=seq(0.1, 0.4, 0.1),
    event=c("transmission", "migration", "transmission", "migration"),
    from.comp=c("S", "I1", "S", "I2"),
    to.comp=c("I1", "I2", "I1", "I2_samp"),
    from.host=c("A", "B", "B_1", "B_1"),
    to.host=c("B", "B_1", "C", "B_sample")
  )
  expect_equal(result, expected)
  
  events <- data.frame(
    time=seq(0.1, 0.4, 0.1),
    event=c("transmission", "transmission", "transmission", "migration"),
    from.comp=c("S", "I1", "S", "I2"),
    to.comp=c("I1", "I1", "I1", "I2_samp"),
    from.host=c("A", "A", "B", "B"),
    to.host=c("B", "B", "C", NA)
  )
  expect_error(.relabel.nodes(events, targets))
})
