library(scrattch.vis)

context("data formatting")

test_that("Matrix is converted to a df", {
  
  library(tasic2016data)
  counts <- tasic2016data::tasic_2016_counts
  df <- mat_to_data_df(counts, cols_are = "sample_names")
  expect_true(class(counts) == "matrix")
  expect_true(class(df) == "data.frame")
  expect_output(str(df), "1809 obs")
  expect_error(expect_output(mat_to_data_df(), NULL))
})


test_that("Data frame is correctly melted", {
  df <- system.file("testdata","helper_dataframe.csv", package = "scrattch.vis")
  df <- read.csv(df)
  melt <- melt_data_df(df, grouping = "sample_name", value_cols = NULL)
  expect_true(length(colnames(melt)) == 3)
  expect_true(class(df) == "data.frame")
  expect_true(class(melt) == "data.frame")
  expect_true(is.factor(melt$sample_name))
  expect_error(expect_output(melt_data_df(), NULL))
})


test_that("Compute stats for samples grouped by one or more annotations", {
  
  library(tasic2016data)
  anno <- tasic2016data::tasic_2016_anno
  df <- system.file("testdata", "helper_melt_df.csv", package = "scrattch.vis")
  df <- read.csv(df)
  stats <- group_stats(df, value_cols = NULL, anno, grouping = c("mouse_line", "cre_reporter"), stat = "mean")
  expect_false(length(colnames(stats)) == 4)
  expect_true(class(df) == "data.frame")
  expect_output(str(df), "sample_name")
  expect_true(class(stats) == "data.frame")
  expect_true(length(stats) > 0)  
  expect_error(expect_output(group_stats(), NULL))
})


test_that("Convert expression data to heatmap colors for plotting", {
  df <- system.file("testdata","helper_dataframe.csv", package = "scrattch.vis")
  df <- read.csv(df)
  output <- data_df_to_colors(df)
  expect_true(class(df) == "data.frame")
  expect_output(str(output), "sample_name")
  expect_output(str(output), "#00008B")
  expect_true(class(output) == "data.frame")
  expect_error(expect_output(data_df_to_colors(), NULL))
})

test_that("Build plot positions for a character vector", {
  df <- system.file("testdata","helper_dataframe.csv", package = "scrattch.vis")
  df <- read.csv(df)
  vector <- as.vector(colnames(df))
  output <- build_vec_pos(vector, vec_name = "genes")
  expect_true(class(vector) == "character")
  expect_true(class(output) == "data.frame")
  expect_output(str(output), "genes")
  expect_true(length(output) == 2)
  expect_error(expect_output(build_vec_pos(), NULL))
})

test_that("Format data provided in list format for scrattch plots", {
 skip('skip this test for now') 
})

