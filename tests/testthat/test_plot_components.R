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
  expect_type(df, "data.frame")
  expect_length(grouping, 1)
  expect_type(grouping, "character")
  expect_type(ymin, "numeric")
  expect_length(build_header_polygons(df, grouping, ymin), 4)
  expect_type(str(build_header_polygons(df, grouping, ymin)), "data.frame")
  expect_error(expect_output(build_header_polygons(), NULL))
})

