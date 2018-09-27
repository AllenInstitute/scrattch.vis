library(scrattch.vis)

context("sankey.R")

test_that("erf() calculates error margins required to generate a sigmoidal curve", {
  expect_equal(erf(0.0000000000000000001),0)
  expect_equal(erf(0.0000000000001), 1.127987e-13)
  expect_equal(erf(10), 1)
  x <- 0.5
  expect_is(x, "numeric")
  expect_error(erf("string"))
  expect_error(erf(tasic_2016_anno))
})




test_that("sigline() Generates x-y coordinates for a sigmoidal line", {
  helper <- readRDS("inst/testdata/helper_sigline_plot.rds")
  plot <- plot(sigline())
  expect_is(sigline(), "data.frame")
  expect_known_output(plot, "inst/testdata/helper_heatmap_legend_plot.rds", print = TRUE)
})
