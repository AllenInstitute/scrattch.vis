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
  out <- sigline()
  helper <- read.csv(system.file("testdata","helper_sigline.csv",  package = "scrattch.vis"))
  plot <- plot(sigline())
  expect_is(sigline(), "data.frame")
  expect_error(sigline(sigfun != "erf"))
  expect_equal(out, helper)
})



test_that("sigribbon() Expands a sigmoidal line to a ribbon by adding ymin", {
  expect_warning(sigribbon(sigline(), height = NA))
  expect_error(sigribbon(sigline(), height = "NA"))
  expect_length(sigribbon(sigline(), height = 5, from = "top"), 3)
  expect_length(sigribbon(sigline(), height = 5, from = "bot"), 3)
  expect_length(sigribbon(sigline(), height = 5, from = "mid"), 3)
  expect_warning(sigribbon(sigline(), height = 1/0, from = "top"))
  expect_warning(sigribbon(sigline(), height = 1/0, from = "top"))
  
})


test_that("make_group_nodes() Makes group nodes for categorical data", {
  library(tasic2016data)
  expect_length(tasic_2016_anno, 17)
  expect_error(make_group_nodes(tasic_2016_anno, "cre_driver"))
  expect_length(make_group_nodes(tasic_2016_anno, "primary_type"), 6)
  out <- make_group_nodes(tasic_2016_anno, "secondary_type")
  expect_output(str(out), "xpos")
})


test_that("make_plot_nodes() generates the plot nodes from make_group_nodes() output", {
  library(tasic2016data)
  input <- make_group_nodes(tasic_2016_anno, "primary_type")
  expect_type(input, "list")
  expect_length(make_plot_nodes(input),13)
  expect_equal(colnames(make_plot_nodes(input)), c("id", "name", "color", "n", "group", "xpos", "xmin", 
                                                   "xmax", "n_groups", "group_pad", "n_cum", "ymin", "ymax"))
})
