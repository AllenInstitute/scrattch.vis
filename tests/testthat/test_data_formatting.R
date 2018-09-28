library(scrattch.vis)

context("data_formatting.R")

test_that(
  "melt_data_df() correctly melts a data.frame to a 3-column format", 
  {
    df <- readRDS(system.file("testdata","test_data_df.RData", package = "scrattch.vis"))
    expect_is(df, "data.frame")
    
    melt_df <- melt_data_df(df, grouping = "sample_name", value_cols = NULL)
    expect_true(ncol(melt_df) == 3)
    
    expect_is(melt_df, "data.frame")
    expect_true(is.character(melt_df$sample_name))
  })


test_that(
  "group_stats() computes statistics for each value column x all combinations of grouping column values", 
  {
    
    library(tasic2016data)
    anno <- tasic2016data::tasic_2016_anno
    
    df <- readRDS(system.file("testdata", "test_data_df.RData", package = "scrattch.vis"))
    
    grouping <- c("mouse_line", "cre_reporter")
    
    all_stats <- group_stats(df, 
                             value_cols = NULL, 
                             anno, 
                             grouping = grouping, 
                             stat = "mean")
    
    n_all_value_cols <- ncol(df) - 1
    n_grouping <- length(grouping)
    n_groups <- nrow(unique(anno[,grouping]))
    
    expect_is(all_stats, "data.frame")
    expect_equal(ncol(all_stats), n_all_value_cols + n_grouping)
    expect_equal(nrow(all_stats), n_groups)
    
    value_cols <- c("Rorb","1110059E24Rik","Nkx2-1")
    select_stats <- group_stats(df, 
                                value_cols = value_cols, 
                                anno, 
                                grouping = grouping, 
                                stat = "mean")
    
    n_value_cols <- length(value_cols)
    n_grouping <- length(grouping)
    n_groups <- nrow(unique(anno[,grouping]))
    
    expect_is(select_stats, "data.frame")
    expect_equal(ncol(select_stats), n_value_cols + n_grouping)
    expect_equal(nrow(select_stats), n_groups)
  }
)


test_that(
  "Convert expression data to heatmap colors for plotting", 
  {
    df <- readRDS(system.file("testdata","test_data_df.RData", package = "scrattch.vis"))
    
    output <- data_df_to_colors(df)
    
    expect_is(output, "data.frame")
    expect_equal(ncol(df), ncol(output))
  }
)

test_that(
  "build_vec_pos() builds a data.frame of plot positions for a vector", 
  {
    genes <- readRDS(system.file("testdata","test_genes.RData", package = "scrattch.vis"))
    
    output <- build_vec_pos(genes, 
                            vec_name = "gene",
                            sort = "none",
                            axis_name = "y")
    
    expect_is(output, "data.frame")
    expect_equal(nrow(output), length(genes))
    expect_equal(names(output), c("gene","y"))
    expect_equal(output$gene, genes)
    
    output_sorted <- build_vec_pos(genes,
                                   vec_name = "gene",
                                   sort = "alpha",
                                   axis_name = "y")
    
    expect_is(output_sorted, "data.frame")
    expect_equal(nrow(output_sorted), length(genes))
    expect_equal(names(output_sorted), c("gene","y"))
    expect_equal(output_sorted$gene, genes[order(genes)])
    
    
  })

test_that("Format data provided in list format for scrattch plots", {
  skip('skip this test for now') 
})
