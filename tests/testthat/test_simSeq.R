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
