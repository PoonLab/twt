test_that("simulate.sequences errors clearly when phangorn isn't installed", {
  skip_if(requireNamespace("phangorn", quietly = TRUE),
          "phangorn is installed -- this test only covers the missing-package path")
  resolved <- list(local.trees = list())
  expect_error(simulate.sequences(resolved), regexp = "requires the 'phangorn' package")
})

test_that("simulate.sequences stitches per-segment sequences into full-length, tip-matching sequences", {
  skip_if_not_installed("phangorn")

  # same tip set in both segments -- matches resolve.arg's own invariant
  # (every segment's tree includes the full sample set, just possibly
  # different topology/branch lengths)
  tree1 <- ape::read.tree(text = "((A:1,B:1):1,C:2);")
  tree2 <- ape::read.tree(text = "(A:2,(B:1,C:1):1);")

  resolved <- list(
    segments = data.frame(start = c(1L, 101L), end = c(100L, 200L)),
    local.trees = list(
      list(start = 1L,   end = 100L, phylo = tree1),
      list(start = 101L, end = 200L, phylo = tree2)
    )
  )

  set.seed(42)
  seqs <- simulate.sequences(resolved, mu = 0.1)

  expect_setequal(names(seqs), c("A", "B", "C"))
  expect_true(all(nchar(seqs) == 200))       # 100 + 100 across the two segments
  expect_true(all(grepl("^[acgt]+$", seqs))) # valid DNA alphabet
})

test_that("simulate.sequences errors on a segment with no valid tree", {
  skip_if_not_installed("phangorn")

  resolved <- list(
    local.trees = list(list(start = 1L, end = 100L, phylo = NULL))
  )
  expect_error(simulate.sequences(resolved), regexp = "no valid tree")
})

test_that("simulate.sequences is reproducible with a seed and returns full-length sequences", {
  skip_if_not_installed("phangorn")
  set.seed(42)
  t1 <- ape::rtree(6, tip.label = paste0("t", 1:6))
  t2 <- ape::rtree(6, tip.label = paste0("t", 1:6))
  resolved <- list(local.trees = list(
    list(start = 1L, end = 40L, phylo = t1),
    list(start = 41L, end = 100L, phylo = t2)))
  a <- simulate.sequences(resolved, mu = 0.05, seed = 1)
  b <- simulate.sequences(resolved, mu = 0.05, seed = 1)
  expect_identical(a, b)
  expect_true(all(nchar(a) == 100))
  expect_setequal(names(a), paste0("t", 1:6))
})

test_that("simulate.sequences errors when a tip is missing from a segment tree", {
  skip_if_not_installed("phangorn")
  set.seed(42)
  t1 <- ape::rtree(6, tip.label = paste0("t", 1:6))
  resolved <- list(local.trees = list(
    list(start = 1L, end = 50L, phylo = t1),
    list(start = 51L, end = 100L, phylo = ape::drop.tip(t1, "t6"))))
  expect_error(simulate.sequences(resolved, mu = 0.05, seed = 1), "not every tip")
})

test_that("simulate.sequences errors when segments do not tile the genome", {
  skip_if_not_installed("phangorn")
  set.seed(42)
  t1 <- ape::rtree(6, tip.label = paste0("t", 1:6))
  resolved <- list(local.trees = list(
    list(start = 1L, end = 40L, phylo = t1),
    list(start = 51L, end = 100L, phylo = t1)))
  expect_error(simulate.sequences(resolved, mu = 0.05, seed = 1), "tile")
})
