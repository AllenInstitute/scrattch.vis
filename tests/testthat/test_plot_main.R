library(scrattch.vis)

context("plot_main.R")

test_that("Barplots of gene expression of individual samples", {
  skip("Skip for now. Fix me!")
})

test_that("heatmap_legend_plot() Builds a heatmap legend plot", {
  #skip("Skip for now. Fix me!")
  plot <- heatmap_legend_plot()
  expect_is(plot, "gg")
  expect_output(str(plot), "list")
})
