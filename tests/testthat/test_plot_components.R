library(scrattch.vis)

context("plot_components.R")

test_that(
  "sci_labels() converts values to scientific notation for ggplot2",
  {
    test_data <- readRDS(system.file("testdata", "test_sci_label_data.RData", package = "scrattch.vis"))
    
    input <- test_data$input
    expected <- test_data$output_3_ggplot2
    
    output <- sci_label(input,
                        sig_figs = 3,
                        type = "ggplot2")
    
    expect_equal(output, expected)
  }
)

test_that(
  "sci_labels() converts values to scientific notation for DT",
  {
    test_data <- readRDS(system.file("testdata", "test_sci_label_data.RData", package = "scrattch.vis"))
    
    input <- test_data$input
    expected <- test_data$output_3_dt
    
    output <- sci_label(input,
                        sig_figs = 3,
                        type = "DT")
    
    expect_equal(output, expected)
  }
)

test_that(
  "sci_labels() converts values to scientific notation for text",
  {
    test_data <- readRDS(system.file("testdata", "test_sci_label_data.RData", package = "scrattch.vis"))
    
    input <- test_data$input
    expected <- test_data$output_3_text
    
    output <- sci_label(input,
                        sig_figs = 3,
                        type = "text")
    
    expect_equal(output, expected)
  }
)

test_that(
  "build_header_polygons() constructs a data.frame of polygon coordinates for each group",
  {
    anno <- readRDS(system.file("testdata", "test_anno_incomplete.RData", package = "scrattch.vis"))
    data <- readRDS(system.file("testdata", "test_data_df.RData", package = "scrattch.vis"))
    
    test_groups <- 1:5
    
    test_anno <- anno[anno$primary_type_id %in% test_groups,]
    
    output_angle <- build_header_polygons(data,
                                          test_anno,
                                          grouping = "primary_type",
                                          group_order = NULL,
                                          ymin = 10,
                                          label_height = 25,
                                          fraction_of_label = 0.1,
                                          poly_type = "angle")
    
    expect_is(output_angle, "data.frame")
    expect_equal(ncol(output_angle), 4)
    expect_equal(nrow(output_angle), length(test_groups) * 4)
    expect_equal(max(output_angle$poly.x), nrow(test_anno))
    expect_equal(min(output_angle$poly.y), 10)
    
    order_groups <- output_angle$id[1:length(test_groups) * 4 - 3]
    order_xmin <- output_angle$poly.x[1:length(test_groups) * 4 - 3]
    order_groups <- order_groups[order(order_xmin)]
    
    expect_equal(order_groups, test_groups)
    
    output_square <- build_header_polygons(data,
                                           test_anno,
                                           grouping = "primary_type",
                                           group_order = NULL,
                                           ymin = 10,
                                           label_height = 25,
                                           fraction_of_label = 0.1,
                                           poly_type = "square")
    
    expect_is(output_square, "data.frame")
    expect_equal(ncol(output_square), 4)
    expect_equal(nrow(output_square), length(test_groups) * 4)
    expect_equal(max(output_square$poly.x), nrow(test_anno))
    expect_equal(min(output_square$poly.y), 10)
    
    
    order_groups <- output_square$id[1:length(test_groups) * 4 - 3]
    order_xmin <- output_square$poly.x[1:length(test_groups) * 4 - 3]
    order_groups <- order_groups[order(order_xmin)]
    
    expect_equal(order_groups, test_groups)
    
    
    # Testing with alternative order
    test_groups <- sample(1:5, 5, replace = FALSE)
    
    test_anno <- anno[anno$primary_type_id %in% test_groups,]
    
    output_angle <- build_header_polygons(data,
                                          test_anno,
                                          grouping = "primary_type",
                                          group_order = test_groups,
                                          ymin = 10,
                                          label_height = 25,
                                          fraction_of_label = 0.1,
                                          poly_type = "angle")
    
    expect_is(output_angle, "data.frame")
    expect_equal(ncol(output_angle), 4)
    expect_equal(nrow(output_angle), length(test_groups) * 4)
    expect_equal(max(output_angle$poly.x), nrow(test_anno))
    expect_equal(min(output_angle$poly.y), 10)
    
    order_groups <- output_angle$id[1:length(test_groups) * 4 - 3]
    order_xmin <- output_angle$poly.x[1:length(test_groups) * 4 - 3]
    order_groups <- order_groups[order(order_xmin)]
    
    expect_equal(order_groups, test_groups)
    
    output_square <- build_header_polygons(data,
                                           test_anno,
                                           grouping = "primary_type",
                                           group_order = test_groups,
                                           ymin = 10,
                                           label_height = 25,
                                           fraction_of_label = 0.1,
                                           poly_type = "square")
    
    expect_is(output_square, "data.frame")
    expect_equal(ncol(output_square), 4)
    expect_equal(nrow(output_square), length(test_groups) * 4)
    expect_equal(max(output_square$poly.x), nrow(test_anno))
    expect_equal(min(output_square$poly.y), 10)
    
    
    order_groups <- output_square$id[1:length(test_groups) * 4 - 3]
    order_xmin <- output_square$poly.x[1:length(test_groups) * 4 - 3]
    order_groups <- order_groups[order(order_xmin)]
    
    expect_equal(order_groups, test_groups)
  }
)


