library(scrattch.vis)

context("Plot main")

test_that("Barplots of gene expression of individual samples", {
  skip("Function needs bug fix :)")
})

test_that("Build a heatmap legend plot",{
  helper <- readRDS("inst/testdata/helper_heatmapLegendPlot.rds")
  plot <- heatmap_legend_plot()
  expect_is(plot, "gg")
  expect_is(helper, "gg")
  expect_output(str(plot), "list")
  expect_known_output(plot, "inst/testdata/helper_heatmapLegendPlot.rds", print = TRUE)
  })