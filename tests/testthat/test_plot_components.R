library(scrattch.vis)

context("Plot components")

test_that("Convert integers to scientific notation labels", {
  my_numbers <- c(100,15.359,32687,.000468)
  expect_type(my_numbers, "double")
  expect_type(sci_label(my_numbers), "character")
  expect_length(sci_label(my_numbers, sig_figs = 3, type = "text"), 4)
  expect_error(expect_output(sci_label(), NULL))
})



test_that("Remove the X-axis (and most other margins)", {
  skip("Function doesn't work at the moment. Skip.")
})


test_that("Build polygons from plot data for fancy headers built into the plot area", {
  df <- system.file("testdata","helper_plot_data.csv", package = "scrattch.vis")
  df <- read.csv(df)
  grouping <- "primary_type"
  ymin <- 5
  expect_output(str(df), "xpos")
  expect_type(df, "list")
  expect_length(grouping, 1)
  expect_type(grouping, "character")
  expect_type(ymin, "double")
  expect_length(build_header_polygons(df, grouping, ymin), 4)
  expect_output(str(build_header_polygons(df, grouping, ymin)), "data.frame")
  expect_error(expect_output(build_header_polygons(), NULL))
})


test_that("Build colorful, rectangular labels for plot headers in plot space", {
  df <- system.file("testdata","helper_plot_data.csv", package = "scrattch.vis")
  df <- read.csv(df)
  grouping <- "primary_type"
  expect_output(str(build_header_labels(df, grouping, ymin = 5)), "data.frame")
  expect_output(str(build_header_labels(df, grouping, ymin = 5)), "label")
  expect_length(build_header_labels(df, grouping, ymin = 5), 6)
  expect_error(expect_output(build_header_labels(), NULL))
})


test_that("Covert hclust objects to segments for use in ggplots", {
  library(tidyverse)
  df <- mtcars %>%
    dist(.) %>%
    hclust(.)
  out <- hclust_to_seg(df, tree.dir = "down", dir.lims = c(0,1))
  expect_type(df, "list")
  expect_output(str(df), "hclust")
  expect_length(out, 4)
  expect_length(hclust_to_seg(df, tree.dir = "up", dir.lims = c(0,1)), 4)
  expect_output(str(out), "data.frame")
})


test_that("Jitter x-y coordinates in a spiral pattern",{
  out <- mtcars %>%
    dist(.) %>%
    hclust(.) %>%
    hclust_to_seg(.)
  expect_length(out, 4)
  expect_output(str(spiral_jitter(out$x,out$y)), "data.frame")
  expect_output(str(spiral_jitter(out$x,out$y)), "y")
  expect_length(spiral_jitter(out$x,out$y), 2)
})