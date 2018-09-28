library(scrattch.vis)

context("plot_components.R")

test_that("sci_labels() converts values to scientific notation for ggplot2",
          {
            skip("Skip for now. Fix me!")
            test_data <-
              readRDS(system.file("testdata", "test_sci_label_data.RData", package = "scrattch.vis"))
            
            input <- test_data$input
            expected <- test_data$output_3_ggplot2
            
            output <- sci_label(input,
                                sig_figs = 3,
                                type = "ggplot2")
            
            expect_equal(output, expected)
          })

test_that("sci_labels() converts values to scientific notation for DT",
          {
            skip("Skip for now. Fix me!")
            test_data <-
              readRDS(system.file("testdata", "test_sci_label_data.RData", package = "scrattch.vis"))
            
            input <- test_data$input
            expected <- test_data$output_3_dt
            
            output <- sci_label(input,
                                sig_figs = 3,
                                type = "DT")
            
            expect_equal(output, expected)
          })

test_that("sci_labels() converts values to scientific notation for text",
          {
            skip("Skip for now. Fix me!")
            test_data <-
              readRDS(system.file("testdata", "test_sci_label_data.RData", package = "scrattch.vis"))
            
            input <- test_data$input
            expected <- test_data$output_3_text
            
            output <- sci_label(input,
                                sig_figs = 3,
                                type = "text")
            
            expect_equal(output, expected)
          })

test_that(
  "build_header_polygons() constructs a data.frame of polygon coordinates for each group",
  {
    skip("Skip for now. Fix me!")
    anno <-
      readRDS(system.file("testdata", "test_anno_incomplete.RData", package = "scrattch.vis"))
    data <-
      readRDS(system.file("testdata", "test_data_df.RData", package = "scrattch.vis"))
    
    test_groups <- 1:5
    
    test_anno <- anno[anno$primary_type_id %in% test_groups, ]
    
    output_angle <- build_header_polygons(
      data,
      test_anno,
      grouping = "primary_type",
      group_order = NULL,
      ymin = 10,
      label_height = 25,
      fraction_of_label = 0.1,
      poly_type = "angle"
    )
    
    expect_is(output_angle, "data.frame")
    expect_equal(ncol(output_angle), 4)
    expect_equal(nrow(output_angle), length(test_groups) * 4)
    expect_equal(max(output_angle$poly.x), nrow(test_anno))
    expect_equal(min(output_angle$poly.y), 10)
    
    order_groups <- output_angle$id[1:length(test_groups) * 4 - 3]
    order_xmin <- output_angle$poly.x[1:length(test_groups) * 4 - 3]
    order_groups <- order_groups[order(order_xmin)]
    
    expect_equal(order_groups, test_groups)
    
    output_square <- build_header_polygons(
      data,
      test_anno,
      grouping = "primary_type",
      group_order = NULL,
      ymin = 10,
      label_height = 25,
      fraction_of_label = 0.1,
      poly_type = "square"
    )
    
    expect_is(output_square, "data.frame")
    expect_equal(ncol(output_square), 4)
    expect_equal(nrow(output_square), length(test_groups) * 4)
    expect_equal(max(output_square$poly.x), nrow(test_anno))
    expect_equal(min(output_square$poly.y), 10)
    
    
    order_groups <- output_square$id[1:length(test_groups) * 4 - 3]
    order_xmin <-
      output_square$poly.x[1:length(test_groups) * 4 - 3]
    order_groups <- order_groups[order(order_xmin)]
    
    expect_equal(order_groups, test_groups)
    
    
    # Testing with alternative order
    test_groups <- sample(1:5, 5, replace = FALSE)
    
    test_anno <- anno[anno$primary_type_id %in% test_groups, ]
    
    output_angle <- build_header_polygons(
      data,
      test_anno,
      grouping = "primary_type",
      group_order = test_groups,
      ymin = 10,
      label_height = 25,
      fraction_of_label = 0.1,
      poly_type = "angle"
    )
    
    expect_is(output_angle, "data.frame")
    expect_equal(ncol(output_angle), 4)
    expect_equal(nrow(output_angle), length(test_groups) * 4)
    expect_equal(max(output_angle$poly.x), nrow(test_anno))
    expect_equal(min(output_angle$poly.y), 10)
    
    order_groups <- output_angle$id[1:length(test_groups) * 4 - 3]
    order_xmin <- output_angle$poly.x[1:length(test_groups) * 4 - 3]
    order_groups <- order_groups[order(order_xmin)]
    
    expect_equal(order_groups, test_groups)
    
    output_square <- build_header_polygons(
      data,
      test_anno,
      grouping = "primary_type",
      group_order = test_groups,
      ymin = 10,
      label_height = 25,
      fraction_of_label = 0.1,
      poly_type = "square"
    )
    
    expect_is(output_square, "data.frame")
    expect_equal(ncol(output_square), 4)
    expect_equal(nrow(output_square), length(test_groups) * 4)
    expect_equal(max(output_square$poly.x), nrow(test_anno))
    expect_equal(min(output_square$poly.y), 10)
    
    
    order_groups <- output_square$id[1:length(test_groups) * 4 - 3]
    order_xmin <-
      output_square$poly.x[1:length(test_groups) * 4 - 3]
    order_groups <- order_groups[order(order_xmin)]
    
    expect_equal(order_groups, test_groups)
  }
)


test_that("Remove the X-axis (and most other margins)", {
  skip("Skip for now. Fix me!")
})

test_that(
  "build_header_labels() builds colorful, rectangular labels for plot headers in plot space",
  {
    skip("Skip for now. Fix me!")
    df <-
      system.file("testdata", "helper_plot_data.csv", package = "scrattch.vis")
    df <- read.csv(df)
    grouping <- "primary_type"
    expect_output(str(build_header_labels(df, grouping, ymin = 5)), "data.frame")
    expect_output(str(build_header_labels(df, grouping, ymin = 5)), "label")
    expect_length(build_header_labels(df, grouping, ymin = 5), 6)
    expect_error(expect_output(build_header_labels(), NULL))
  }
)


test_that("hclust_to_seg() Converts hclust objects to segments for use in ggplots",
          {
            skip("Skip for now. Fix me!")
            library(tidyverse)
            df <- mtcars %>%
              dist(.) %>%
              hclust(.)
            out <- hclust_to_seg(df, tree.dir = "down", dir.lims = c(0, 1))
            expect_type(df, "list")
            expect_output(str(df), "hclust")
            expect_length(out, 4)
            expect_length(hclust_to_seg(df, tree.dir = "up", dir.lims = c(0, 1)), 4)
            expect_output(str(out), "data.frame")
          })


test_that("spiral_jitter() Generates Jitter x-y coordinates in a spiral pattern",
          {
            out <- mtcars %>%
              dist(.) %>%
              hclust(.) %>%
              hclust_to_seg(.)
            expect_length(out, 4)
            expect_is(out, "data.frame")
            expect_warning(spiral_jitter(out$x, out$y, n = 1))
            expect_error(spiral_jitter(out$x, out$y, n = NA))
            expect_error(spiral_jitter(out$x, out$y, n = "string"))
            expect_error(spiral_jitter(out$x, out$y, max_n = "string"))
            expect_error(spiral_jitter(out$x, out$y, ratio = "notgolden"))
          })
