library(scrattch.vis)

context("sankey.R")

test_that("erf() calculates error margins required to generate a sigmoidal curve",
          {
            expect_equal(erf(0.0000000000000000001), 0)
            expect_equal(erf(0.0000000000001), 1.127987e-13)
            expect_equal(erf(10), 1)
            x <- 0.5
            expect_is(x, "numeric")
            expect_error(erf("string"))
            expect_error(erf(tasic_2016_anno))
          })



test_that("sigline() Generates x-y coordinates for a sigmoidal line", {
  out <- sigline()
  helper <-
    readRDS(system.file("testdata", "helper_sigline.RData",  package = "scrattch.vis"))
  plot <- plot(sigline())
  expect_is(sigline(), "data.frame")
  expect_error(sigline(sigfun != "erf"))
  expect_equal(out, helper)
})



test_that("sigribbon() Expands a sigmoidal line to a ribbon by adding ymin", {
  expect_warning(sigribbon(sigline(), height = NA))
  expect_error(sigribbon(sigline(), height = "NA"))
  expect_error(sigribbon())
  expect_length(sigribbon(sigline(), height = 5, from = "top"), 3)
  expect_length(sigribbon(sigline(), height = 5, from = "bot"), 3)
  expect_length(sigribbon(sigline(), height = 5, from = "mid"), 3)
  expect_warning(sigribbon(sigline(), height = 1 / 0, from = "top"))
  expect_warning(sigribbon(sigline(), height = 1 / 0, from = "top"))
  
})



test_that("make_group_nodes() Makes group nodes for categorical data", {
  library(tasic2016data)
  expect_length(tasic_2016_anno, 17)
  expect_error(make_group_nodes())
  expect_error(make_group_nodes(tasic_2016_anno, "cre_driver"))
  expect_error(make_group_nodes(tasic_2016_anno, "primary_type"))
  out <-
    make_group_nodes(tasic_2016_anno, c("primary_type", "secondary_type"))
  expect_output(str(out), "xpos")
})


test_that("make_plot_nodes() generates the plot nodes from make_group_nodes() output",
          {
            library(tasic2016data)
            input <-
              make_group_nodes(tasic_2016_anno, c("primary_type", "secondary_type"))
            expect_type(input, "list")
            expect_error(make_plot_nodes())
            expect_length(make_plot_nodes(input), 13)
            expect_equal(
              colnames(make_plot_nodes(input)),
              c(
                "id",
                "name",
                "color",
                "n",
                "group",
                "xpos",
                "xmin",
                "xmax",
                "n_groups",
                "group_pad",
                "n_cum",
                "ymin",
                "ymax"
              )
            )
          })


test_that("make_group_links() generates the group links from the plot notes",
          {
            library(tasic2016data)
            library(tidyverse)
            inputnodes <-
              make_group_nodes(tasic_2016_anno, c("primary_type", "secondary_type")) %>%
              make_plot_nodes(.)
            
            expect_length(make_group_links(
              tasic_2016_anno,
              c("primary_type", "secondary_type"),
              inputnodes
            ), 16)
            expect_error(make_group_links())
            expect_error(make_group_links(tasic_2016_anno, "primary_type", inputnodes))
            expect_is(make_group_links(
              tasic_2016_anno,
              c("primary_type", "secondary_type"),
              inputnodes
            ), "tbl_df")
            expect_equal(
              colnames(make_group_links(
                tasic_2016_anno,
                c("primary_type", "secondary_type"),
                inputnodes
              )),
              c(
                "group1_id",
                "group2_id",
                "group1_label",
                "group2_label",
                "group1_color",
                "group2_color",
                "n",
                "group1",
                "group2",
                "x",
                "group1_min",
                "xend",
                "group2_min",
                "y",
                "yend",
                "link_id"
              )
            )
          })



test_that("make_plot_links() generates the link ids",
          {
            library(tasic2016data)
            library(tidyverse)
            input <-
              make_group_nodes(tasic_2016_anno, c("primary_type", "secondary_type")) %>%
              make_plot_nodes(.) %>%
              make_group_links(tasic_2016_anno, c("primary_type", "secondary_type"), .)
            
            expect_length(input, 16)
            expect_length(make_plot_links(input), 5)
            expect_equal(colnames(make_plot_links(input)),
                         c("x", "y", "ymin", "fill", "link_id"))
            expect_type(make_plot_links(input), "list")
            expect_type(make_plot_links(input)$link_id, "integer")
            expect_error(make_plot_links())
          })


test_that("build_river_plot() generates a river plot",
          {
            library(tasic2016data)
            expect_error(build_river_plot())
            expect_error(build_river_plot(tasic_2016_anno, "primary_type"))
            expect_warning(build_river_plot(tasic_2016_anno, c("primary_type", "secondary_type")))
            group <- as.factor(c("primary_type", "secondary_type"))
            expect_error(build_river_plot(tasic_2016_anno, group))
            expect_length(build_river_plot(tasic_2016_anno, c("primary_type", "secondary_type")), 9)
            helper <-
              readRDS(
                system.file("testdata", "helper_riverplot_tasicdata.RData", package = "scrattch.vis")
              )
            out <-
              build_river_plot(tasic_2016_anno, c("primary_type", "secondary_type"))
            expect_equal(helper, out)
          })
