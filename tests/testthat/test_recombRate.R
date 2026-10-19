# Tests for the per-nucleotide recombination rate (.resolve.rho and sim.arg).

mock.model <- function(rho.per.site = NULL) {
  params <- list(simTime = 10)
  if (!is.null(rho.per.site)) params$rho.per.site <- rho.per.site
  list(get.parameters = function() params)
}

test_that("per-nucleotide rate is multiplied by sequence length", {
  rho <- .resolve.rho(1e-4, FALSE, 1e-8, 9000L, mock.model())
  expect_equal(rho, 9e-5)
})

test_that("same per-nucleotide rate scales with genome length", {
  r1 <- .resolve.rho(1e-4, FALSE, 1e-8, 9000L, mock.model())
  r2 <- .resolve.rho(1e-4, FALSE, 1e-8, 18000L, mock.model())
  expect_equal(r2, 2 * r1)
  expect_equal(r1 / 9000, r2 / 18000)   # per-nucleotide rate unchanged
})

test_that("rho given alone is genome-level and not rescaled", {
  expect_equal(.resolve.rho(0.3, TRUE, NULL, 9000L, mock.model()), 0.3)
  expect_equal(.resolve.rho(0.3, TRUE, NULL, 18000L, mock.model()), 0.3)
})

test_that("rho and rho.per.site together are rejected", {
  expect_error(.resolve.rho(0.3, TRUE, 1e-8, 9000L, mock.model()),
               regexp = "not both")
})

test_that("rho.per.site is read from the model parameters when no argument is given", {
  rho <- .resolve.rho(1e-4, FALSE, NULL, 9000L, mock.model(rho.per.site = 2e-8))
  expect_equal(rho, 2e-8 * 9000)
})

test_that("an explicit rho overrides the YAML rate", {
  rho <- .resolve.rho(0.3, TRUE, NULL, 9000L, mock.model(rho.per.site = 2e-8))
  expect_equal(rho, 0.3)
})

test_that("an explicit rho.per.site overrides the YAML rate", {
  rho <- .resolve.rho(1e-4, FALSE, 5e-9, 9000L, mock.model(rho.per.site = 2e-8))
  expect_equal(rho, 5e-9 * 9000)
})

test_that("the default rho is used when nothing else is specified", {
  expect_equal(.resolve.rho(1e-4, FALSE, NULL, 9000L, mock.model()), 1e-4)
})

test_that("invalid per-nucleotide rates are rejected", {
  expect_error(.resolve.rho(1e-4, FALSE, -1e-8, 9000L, mock.model()),
               regexp = "non-negative")
  expect_error(.resolve.rho(1e-4, FALSE, c(1e-8, 2e-8), 9000L, mock.model()),
               regexp = "single")
  expect_error(.resolve.rho(1e-4, FALSE, NULL, 9000L,
                            mock.model(rho.per.site = "fast")),
               regexp = "rho.per.site")
})

test_that("a zero per-nucleotide rate gives no recombination", {
  expect_equal(.resolve.rho(1e-4, FALSE, 0, 9000L, mock.model()), 0)
})
