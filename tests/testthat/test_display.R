# Smoke tests for print, summary and plot methods. They check that the output is
# produced without errors, not what the plot looks like. Plots are drawn to a
# null device so nothing is written to disk.

fresh.settings <- function() yaml.load_file("test_SIR.yaml")

with.null.device <- function(code) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  force(code)
}

make.objects <- function() {
  mod <- Model$new(fresh.settings())
  set.seed(1)
  dyn <- sim.dynamics(mod)
  outer <- sim.outer.tree(dyn)
  list(mod = mod, dyn = dyn, outer = outer)
}

test_that("Model print and summary show the model", {
  mod <- Model$new(fresh.settings())
  expect_output(print(mod), "twt Model")
  expect_output(summary(mod), "twt Model")
})

test_that("plot.Model draws without error", {
  mod <- Model$new(fresh.settings())
  with.null.device(expect_error(plot(mod), NA))
})

test_that("dynamics objects print and plot", {
  o <- make.objects()
  expect_output(print(o$dyn), "twt Dynamics object")
  with.null.device(expect_error(plot(o$dyn), NA))
})

test_that("plotting uncounted dynamics warns and counts first", {
  mod <- Model$new(fresh.settings())
  set.seed(1)
  dyn <- sim.dynamics(mod, counted = FALSE)
  with.null.device(
    expect_warning(plot(dyn), "plot.dynamics requires counts")
  )
})

test_that("OuterTree prints and plots", {
  o <- make.objects()
  expect_output(print(o$outer), "twt OuterTree")
  with.null.device(expect_error(plot(o$outer), NA))
  with.null.device(expect_error(plot(o$outer, dynamics = o$dyn), NA))
})

test_that("HostSet and Host print", {
  o <- make.objects()
  hosts <- o$outer$get.sampled()
  expect_output(print(hosts), "Host objects")
  expect_output(print(hosts$get.hosts()[[1]]), "twt Host")
})

test_that("InnerTree prints", {
  o <- make.objects()
  inner <- sim.inner.tree(o$outer)
  expect_output(print(inner), "twt InnerTree")
})
