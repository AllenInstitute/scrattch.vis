library(scrattch.vis)

context("plot_main.R")

test_that("Barplots of gene expression of individual samples", {
  skip("Skip for now. Fix me!")
})

test_that("Build a heatmap legend plot", {
  #skip("Skip for now. Fix me!")
  helper <-
    readRDS(
      system.file("testdata", "helper_heatmap_legend_plot.rds", package = "scrattch.vis")
    )
  plot <- heatmap_legend_plot()
  expect_is(plot, "gg")
  expect_is(helper, "gg")
  expect_output(str(plot), "list")
  expect_equal(plot, helper)
})
