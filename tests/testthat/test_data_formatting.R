library(scrattch.vis)

context("data formatting")

test_that("Matrix is converted to a df", {

  library(tasic2016data)
  counts <- tasic2016data::tasic_2016_counts
  df <- mat_to_data_df(counts, cols_are = "sample_names")
  expect_true(class(df) == "data.frame")
  expect_output(str(df), "1809 obs")            
})


test_that("Data frame is correctly melted", {
  df <- system.file("inst","testdata","helper_dataframe.csv", package = "scrattch.vis")
  df <- read.csv(df)
  melt <- melt_data_df(df, grouping = "sample_name", value_cols = NULL)
  expect_true(length(colnames(melt)) == 3)
  expect_true(class(melt) == "data.frame")
  expect_true(is.factor(melt$sample_name))
})


test_that("Compute stats for samples grouped by one or more annotations", {
  library(tasic2016data)
  anno <- tasic2016data::tasic_2016_anno
  df <- system.file("inst", "testdata", "helper_melt_df.csv", package = "scrattch.vis")
  df <- read.csv(df)
  stats <- group_stats(df, value_cols = NULL, anno, grouping = c("mouse_line", "cre_reporter"), stat = "mean")
  expect_false(length(colnames(stats)) == 4)
  expect_true(class(stats) == "data.frame")
  expect_true(length(test) > 0)
})

