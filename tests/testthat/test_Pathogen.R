# --- Pathogen ---------------------------------------------------------------

test_that("Pathogen: accessors return what the constructor stored", {
  p <- Pathogen$new(name = "a", start.time = 1, end.time = 2)
  expect_equal(p$get.name(), "a")
  expect_equal(p$get.start.time(), 1)
  expect_equal(p$get.end.time(), 2)
  expect_length(p$get.parents(), 0)
  expect_length(p$get.children(), 0)
})

test_that("Pathogen: set.start.time changes the start time", {
  p <- Pathogen$new(name = "a", start.time = 1, end.time = 2)
  p$set.start.time(0.5)
  expect_equal(p$get.start.time(), 0.5)
})

test_that("Pathogen: print shows the object summary", {
  p <- Pathogen$new(name = "a", start.time = 1, end.time = 2)
  expect_output(print(p), "twt Pathogen")
  expect_output(print(p), "Start time")
})

test_that("Pathogen: set.parent accepts a single Pathogen or a list", {
  child <- Pathogen$new(name = "c")
  p1 <- Pathogen$new(name = "p1")
  p2 <- Pathogen$new(name = "p2")

  child$set.parent(p1)
  expect_length(child$get.parents(), 1)
  expect_identical(child$get.parent(), p1)

  child$set.parent(list(p1, p2))
  expect_length(child$get.parents(), 2)
  expect_identical(child$get.parent(), p1)   # first parent
})

test_that("Pathogen: set.parent rejects things that are not Pathogens", {
  child <- Pathogen$new(name = "c")
  expect_error(child$set.parent(list(Pathogen$new(name = "p"), "x")),
               "all elements of `parent` list must be Pathogen objects",
               fixed = TRUE)
  expect_error(child$set.parent("x"),
               "`parent` must be a Pathogen object or list of Pathogen objects",
               fixed = TRUE)
})

test_that("Pathogen: breakpoint can be set and read back", {
  p <- Pathogen$new(name = "a")
  p$set.breakpoint(42L)
  expect_equal(p$get.breakpoint(), 42L)
})
