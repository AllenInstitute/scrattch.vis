# Note on order of operations:
# For plotting, data should first be filtered based on annotations and group_order (if provided)
# Then, max values can be obtained for max_val
# Next, data can be rescaled based on log_scale
# Then, global or within-row normalization can occur (during color selection for heatmaps)
# 1. Filter
# 2. Compute Max Vals
# 3. Scale
# 4. Normalize
# FilMScaN

#' Barplots of gene expression of individual samples
#' 
#' This function will generate plots similar to those in Figure 3a-c of Tasic, et al. (2016).
#' 
#' @param data A data frame containing gene expression values. The first column should be sample_name
#' @param anno Sample annotations. The first column should be sample_name, and each annotation should have \_id, \_label, and \_color columns
#' @param genes A character vector containing gene symbols to be plotted. 
#' @param grouping A character string specifying the desc base (column) that should be used to group cells. 
#' @param group_order Optional: Explicit specification of group order by supplying a vector of group_ids.
#' @param log_scale Logical , determines if data is log scaled before plotting. Default = FALSE.
#' @param font_size numeric object, the font size (in pts) used to make the plot.
#' @param label_height numeric object, Percent of the plot height that should be used for the labels (0 to 100). Default is 25.
#' @param label_type Label shape, "angle" or "square"
#' @param max_width numeric object, percent of plot width that should be used for maximum expression values (0 to 100). Default is 10.
#' @param bg_color plot background color. Default is a light blue ("#ADCFE0")
#' @param return_type What values to return - can be "plot", "data", or "both". Default is "plot".
#' 
#' @return a ggplot2 plot object
#'
sample_bar_plot <- function(data,
                            anno,
                            genes,
                            grouping,
                            group_order = NULL,
                            log_scale = FALSE,
                            font_size = 7, 
                            label_height = 25, 
                            label_type = "angle",
                            max_width = 10,
                            bg_color = "#ADCFE0",
                            return_type = "plot") {
  
  # Reverse order so genes are plotted from the top down
  genes <- rev(genes)
  
  group_cols <- group_columns(grouping)
  
  # Filter data to genes and samples in anno
  gene_data <- filter_gene_data(data, 
                                genes, 
                                anno, 
                                group_cols,
                                group_order, 
                                "sample_name")
  
  # Filter annotations if group_order is provided
  if(!is.null(group_order)) {
    anno <- anno[anno[[group_cols$id]] %in% group_order,]
  }
  
  # Get maximum values for each gene before rescaling to plot space.
  max_vals_unscaled <- max_gene_vals(gene_data, genes)
  
  # Left-join data to anno. This will ensure that data is filtered for the cells provided in anno
  plot_data <- dplyr::left_join(anno, gene_data, by = "sample_name")
  
  # Rescale values if needed
  if(log_scale) {
    plot_data <- scale_gene_data(plot_data, genes, scale_type = "log10")
  }
  
  # Add x-positions for each sample
  plot_data <- add_sample_xpos(plot_data,
                               group_cols = group_cols,
                               group_order = group_order)
  
  # Compute basic count stats that are used downstream
  # n_stats$genes, n_stats$samples, and n_stats$groups
  n_stats <- get_n_stats(plot_data, group_cols, genes)
  
  # Calculate segments for cluster sepration lines
  parsed_id <- rlang::parse_expr(group_cols$id)
  segment_lines <- plot_data %>%
    dplyr::group_by(!!parsed_id) %>%
    dplyr::summarise(x = max(xpos) )
  
  # The background of the plot is a rectangular object.
  background_data <- data.frame(xmin = 0, 
                                xmax = n_stats$samples, 
                                ymin = 1, 
                                ymax = n_stats$genes + 1, 
                                fill = bg_color)
  
  ### Plot setup
  p <-ggplot2::ggplot() +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0), 
                                breaks = 1:n_stats$genes + 0.45, 
                                labels = genes) +
    ggplot2::theme_classic(base_size = font_size) +
    theme(axis.text = ggplot2::element_text(size=ggplot2::rel(1), face = "italic"),
          axis.ticks = ggplot2::element_blank(),
          axis.line = ggplot2::element_blank(),
          axis.title = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_blank()) +
    ggplot2::geom_rect(data = background_data,
                       ggplot2::aes(xmin = xmin,
                                    xmax = xmax,
                                    ymin = ymin,
                                    ymax = ymax,
                                    fill = fill)) +
    ggplot2::geom_segment(data = segment_lines,
                          ggplot2::aes(x = x, 
                                       xend = x, 
                                       y = 1, 
                                       yend = n_stats$genes + 1),
                          linewidth = 0.2, 
                          color = "gray60", 
                          linetype = "dashed")
  
  ### plot the bars for each gene
  for(i in 1:n_stats$genes) {
    
    gene <- genes[i]
    
    plot_data[[paste0(gene, "_min")]] <- i
    plot_data[[paste0(gene, "_max")]] <- scale_values_plot_space(plot_data[[gene]],
                                                 min_ps = i)
    
    # plot the rectangles for the barplots
    parsed_min <- rlang::parse_expr(paste0(gene, "_min"))
    parsed_max <- rlang::parse_expr(paste0(gene, "_max"))
    parsed_color <- rlang::parse_expr(group_cols$color)
    
    p <- p + 
      ggplot2::geom_rect(data = plot_data,
                         ggplot2::aes(xmin = xpos - 1,
                                      xmax = xpos,
                                      ymin = !!parsed_min,
                                      ymax = !!parsed_max,
                                      fill = !!parsed_color))
    
  }
  
  ### Cluster labels at the top of the plot
  
  # build_header_polygons from plot_components.R
  header_polygons <- build_header_polygons(data = gene_data, 
                                           anno = anno,
                                           grouping = grouping,
                                           group_order = group_order,
                                           ymin = n_stats$genes + 1, 
                                           label_height = label_height, 
                                           poly_type = label_type)
  
  # Build the cell type label rectangles from plot_components.R
  header_labels <- build_header_labels(data = plot_data, 
                                       grouping = grouping,
                                       ymin = n_stats$genes + 1, 
                                       label_height = label_height, 
                                       label_type = label_type)
  
  p <- ggplot_header_labels(p,
                            header_labels = header_labels,
                            header_polygons = header_polygons,
                            font_size = font_size)
  
  ### Scale bar elements
  p <- ggplot_scale_bars(p, 
                         n_genes = n_stats$genes, 
                         n_samples = n_stats$samples,
                         extent = 0.9)
  
  ### Maximum value labels at the right edge of the plot
  max_val_dfs <- build_max_dfs(n_stats, width_stat = "samples", max_vals_unscaled, max_width)
  
  p <- ggplot_max_vals(p,
                       n_stats = n_stats,
                       max_val_dfs = max_val_dfs,
                       font_size = font_size)
  
  if(return_type == "plot") {
    return(p)
  } else if(return_type == "data") {
    return(list(plot_data = plot_data,
                header_labels = header_labels,
                header_polygons = header_polygons,
                max_val_dfs = max_val_dfs,
                n_stats = n_stats))
  } else if(return_type == "both") {
    return(list(plot = p,
                plot_data = plot_data,
                header_labels = header_labels,
                header_polygons = header_polygons,
                max_val_dfs = max_val_dfs,
                n_stats = n_stats))
  }
}

#' Heatmaps of gene expression for individual samples
#' 
#' @param data A data frame containing gene expression values. The first column should be sample_name
#' @param anno Sample annotations. The first column should be sample_name, and each annotation should have \_id, \_label, and \_color columns
#' @param genes A character vector containing gene symbols to be plotted. 
#' @param grouping A character vector specifying the desc base that should be used to group cells
#' @param group_order Optional: Explicit specification of group order by supplying a vector of group_ids.
#' @param log_scale Logical , determines if data is log scaled before plotting. Default = FALSE.
#' @param normalize_rows Logical, determines if heatmaps will be normalized for each gene. Default = FALSE.
#' @param font_size numeric object, the font size (in pts) used to make the plot.
#' @param label_height numeric object, Percent of the plot height that should be used for the labels (0 to 100). Default is 25.
#' @param label_type Label shape, "angle" or "square"
#' @param max_width numeric object, percent of plot width that should be used for maximum expression values (0 to 100). Default is 10.
#' @param return_type What values to return - can be "plot", "data", or "both". Default is "plot".
#' 
#' @return a ggplot2 plot object
#'
sample_heatmap_plot <- function(data,
                                anno,
                                genes,
                                grouping,
                                group_order = NULL,
                                log_scale = TRUE, 
                                normalize_rows = FALSE,
                                colorset = NULL,
                                font_size = 7, 
                                label_height = 25,
                                label_type = "angle",
                                max_width = 10,
                                return_type = "plot") {
  
  genes <- rev(genes)
  
  group_cols <- group_columns(grouping)
  
  # Filter data to genes and samples in anno
  gene_data <- filter_gene_data(data, 
                                genes, 
                                anno, 
                                group_cols,
                                group_order, 
                                "sample_name")
  
  # Filter annotations if group_order is provided
  if(!is.null(group_order)) {
    anno <- anno[anno[[group_cols$id]] %in% group_order,]
  }
  
  # Get maximum values for each gene before rescaling to plot space.
  max_vals_unscaled <- max_gene_vals(gene_data, genes)
  
  # Left-join data to anno. This will ensure that data is filtered for the cells provided in anno
  plot_data <- dplyr::left_join(anno, gene_data, by = "sample_name")
  
  # Rescale values if needed
  if(log_scale) {
    plot_data <- scale_gene_data(plot_data, genes, scale_type = "log10")
  }
  
  # Convert the data values to heatmap colors
  plot_data <- data_df_to_colors(plot_data,
                                 value_cols = genes,
                                 per_col = normalize_rows,
                                 colorset = colorset)
  
  # Add x-positions for each sample
  plot_data <- add_sample_xpos(plot_data,
                               group_cols = group_cols,
                               group_order = group_order)
  
  # Compute basic count stats that are used downstream
  # n_stats$genes, n_stats$samples, and n_stats$groups
  n_stats <- get_n_stats(plot_data, group_cols, genes)
  
  # Plot setup
  p <-ggplot2::ggplot(plot_data) +
    ggplot2::scale_fill_identity() +
    ggplot2::theme_classic(base_size = font_size) +
    theme(axis.text = ggplot2::element_text(size=ggplot2::rel(1), face = "italic"),
          axis.ticks = ggplot2::element_blank(),
          axis.line = ggplot2::element_blank(),
          axis.title = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_blank()) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0), 
                                breaks = 1:n_stats$genes + 0.45, 
                                labels = genes)
  
  # plot the rectangles for each gene
  for(i in seq_along(genes)) {
    
    # Explicit specification of rectangle components
    # prevents overplotting and makes the function faster.
    rect_data <- data.frame(xmin = plot_data$xpos - 1,
                            xmax = plot_data$xpos,
                            ymin = i,
                            ymax = i + 1,
                            fill = plot_data[[genes[i]]])
    
    # plot the rectangles for the heatmap
    p <- p + ggplot2::geom_rect(data = rect_data,
                                ggplot2::aes(xmin = xmin, 
                                             xmax = xmax, 
                                             ymin = ymin, 
                                             ymax = ymax + 1, 
                                             fill = fill))
    
  }
  
  ### Cluster labels at the top of the plot
  
  # build_header_polygons from plot_components.R
  header_polygons <- build_header_polygons(data = gene_data, 
                                           anno = anno,
                                           grouping = grouping,
                                           group_order = group_order,
                                           ymin = n_stats$genes + 1, 
                                           label_height = label_height, 
                                           poly_type = label_type)
  
  # Build the cell type label rectangles from plot_components.R
  header_labels <- build_header_labels(data = plot_data, 
                                       grouping = grouping,
                                       ymin = n_stats$genes + 1, 
                                       label_height = label_height, 
                                       label_type = label_type)
  
  p <- ggplot_header_labels(p,
                            header_labels = header_labels,
                            header_polygons = header_polygons,
                            font_size = font_size)
  
  ### Scale bar elements
  p <- ggplot_scale_bars(p, 
                         n_genes = n_stats$genes, 
                         n_samples = n_stats$samples,
                         extent = 0.9)
  
  ### Maximum value labels at the right edge of the plot
  max_val_dfs <- build_max_dfs(n_stats, width_stat = "samples", max_vals_unscaled, max_width)
  
  p <- ggplot_max_vals(p,
                       n_stats = n_stats,
                       max_val_dfs = max_val_dfs,
                       font_size = font_size)
  
  if(return_type == "plot") {
    return(p)
  } else if(return_type == "data") {
    return(list(plot_data = plot_data,
                header_labels = header_labels,
                header_polygons = header_polygons,
                max_val_dfs = max_val_dfs,
                n_stats = n_stats))
  } else if(return_type == "both") {
    return(list(plot = p,
                plot_data = plot_data,
                header_labels = header_labels,
                header_polygons = header_polygons,
                max_val_dfs = max_val_dfs,
                n_stats = n_stats))
  }
  
}


#' Fire Heatmaps of gene expression for individual samples
#' 
#' @param data A data frame containing gene expression values. The first column should be sample_name
#' @param anno Sample annotations. The first column should be sample_name, and each annotation should have \_id, \_label, and \_color columns
#' @param genes A character vector containing gene symbols to be plotted. 
#' @param grouping A character vector specifying the desc base that should be used to group cells
#' @param group_order Optional: Explicit specification of group order by supplying a vector of group_ids.
#' @param log_scale Logical , determines if data is log scaled before plotting. Default = FALSE.
#' @param normalize_rows Logical, determines if heatmaps will be normalized for each gene. Default = FALSE.
#' @param font_size numeric object, the font size (in pts) used to make the plot.
#' @param label_height numeric object, Percent of the plot height that should be used for the labels (0 to 100). Default is 25.
#' @param max_width numeric object, percent of plot width that should be used for maximum expression values (0 to 100). Default is 10.
#' @param return_type What values to return - can be "plot", "data", or "both". Default is "plot".
#' 
#' @return a ggplot2 plot object
#'
sample_fire_plot <- function(data,
                             anno,
                             genes,
                             grouping,
                             group_order = NULL,
                             log_scale = TRUE, 
                             normalize_rows = FALSE,
                             colorset = NULL,
                             top_values = "lowest",
                             font_size = 7, 
                             label_height = 25,
                             max_width = 10,
                             return_type = "plot") {
  
  genes <- rev(genes)
  
  group_cols <- group_columns(grouping)
  
  # Filter data to genes and samples in anno
  gene_data <- filter_gene_data(data, 
                                genes, 
                                anno, 
                                group_cols,
                                group_order, 
                                "sample_name")
  
  # Filter annotations if group_order is provided
  if(!is.null(group_order)) {
    anno <- anno[anno[[group_cols$id]] %in% group_order,]
  }
  
  # Get maximum values for each gene before rescaling to plot space.
  max_vals_unscaled <- max_gene_vals(gene_data, genes)
  
  # Left-join data to anno. This will ensure that data is filtered for the cells provided in anno
  plot_data <- dplyr::left_join(anno, gene_data, by = "sample_name")
  
  # Rescale values if needed
  if(log_scale) {
    plot_data <- scale_gene_data(plot_data, genes, scale_type = "log10")
  }
  
  # Add x-positions for each group
  plot_data <- add_group_xpos(plot_data,
                              group_cols = group_cols,
                              group_order = group_order)
  
  # Compute basic count stats that are used downstream
  # n_stats$genes, n_stats$samples, and n_stats$groups
  n_stats <- get_n_stats(plot_data, group_cols, genes)
  
  # Build the cell type label rectangles from plot_components.R
  header_labels <- build_header_labels(data = plot_data, 
                                       grouping = grouping,
                                       group_order = group_order,
                                       ymin = n_stats$genes + 1, 
                                       label_height = label_height, 
                                       label_type = "simple")
  
  # Plot setup
  p <-ggplot2::ggplot(data) +
    ggplot2::scale_fill_identity() +
    ggplot2::theme_classic(base_size = font_size) +
    theme(axis.text = ggplot2::element_text(size = ggplot2::rel(1), face = "italic"),
          axis.ticks = ggplot2::element_blank(),
          axis.line = ggplot2::element_blank(),
          axis.title = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_blank()) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0), 
                                breaks = 1:n_stats$genes + 0.45, 
                                labels = genes)
  
  # plot the rectangles for each gene
  for (i in seq_along(genes)) {
    
    gene <- genes[i]
    
    gene_values <- plot_data %>%
      dplyr::select(dplyr::one_of("sample_name", gene))
    
    gene_colors <- data_df_to_colors(gene_values,
                                     value_cols = gene,
                                     per_col = normalize_rows,
                                     colorset = colorset)
    
    names(gene_colors)[names(gene_colors) == gene] <- "plot_fill"
    
    if (top_values == "highest") {
      parsed_gene <- rlang::parse_expr(gene)
      rect_data <- plot_data %>%
        dplyr::left_join(gene_colors, by = "sample_name") %>%
        dplyr::arrange(!!parsed_gene) %>%
        dplyr::group_by(xpos) %>%
        dplyr::mutate(group_n = dplyr::n()) %>%
        dplyr::mutate(xmin = xpos - 0.5,
                      xmax = xpos + 0.5,
                      ymin = seq(i, i + 1, length.out = group_n[1] + 1)[-(group_n[1] + 1)],
                      ymax = seq(i, i + 1, length.out = group_n[1] + 1)[-1])
    } else {
      arr_gene <- paste0("-",gene)
      parsed_arr_gene <- rlang::parse_expr(arr_gene)
      rect_data <- plot_data %>%
        dplyr::left_join(gene_colors, by = "sample_name") %>%
        dplyr::arrange(!!parsed_arr_gene) %>%
        dplyr::group_by(xpos) %>%
        dplyr::mutate(group_n = dplyr::n()) %>%
        dplyr::mutate(xmin = xpos - 0.5,
                      xmax = xpos + 0.5,
                      ymin = seq(i, i + 1, length.out = group_n[1] + 1)[-(group_n[1] + 1)],
                      ymax = seq(i, i + 1, length.out = group_n[1] + 1)[-1])
    }
    
    
    # plot the rectangles for the heatmap
    p <- p + ggplot2::geom_rect(data = rect_data,
                                ggplot2::aes(xmin = xmin, 
                                             xmax = xmax, 
                                             ymin = ymin, 
                                             ymax = ymax, 
                                             fill = plot_fill))
    
  }
  
  ### Cluster labels at the top of the plot
  
  # Build the cell type label rectangles from plot_components.R
  header_labels <- build_header_labels(data = plot_data, 
                                       grouping = grouping,
                                       ymin = n_stats$genes + 1, 
                                       label_height = label_height, 
                                       label_type = "simple")
  
  p <- ggplot_header_labels(p,
                            header_labels = header_labels,
                            header_polygons = NULL,
                            font_size = font_size)
  
  ### Maximum value labels at the right edge of the plot
  max_val_dfs <- build_max_dfs(n_stats, 
                               width_stat = "groups", 
                               max_vals_unscaled, 
                               max_width)
  
  p <- ggplot_max_vals(p,
                       n_stats = n_stats,
                       width_stat = "groups",
                       max_val_dfs = max_val_dfs,
                       font_size = font_size)
  
  if(return_type == "plot") {
    return(p)
  } else if(return_type == "data") {
    return(list(plot_data = plot_data,
                header_labels = header_labels,
                header_polygons = header_polygons,
                max_val_dfs = max_val_dfs,
                n_stats = n_stats))
  } else if(return_type == "both") {
    return(list(plot = p,
                plot_data = plot_data,
                header_labels = header_labels,
                max_val_dfs = max_val_dfs,
                n_stats = n_stats))
  }
  
}


#' Violin plots of gene expression for grouped samples
#' 
#' @param data A data frame containing gene expression values. The first column should be sample_name
#' @param anno Sample annotations. The first column should be sample_name, and each annotation should have \_id, \_label, and \_color columns
#' @param genes A character vector containing gene symbols to be plotted. 
#' @param grouping A character vector specifying the desc base that should be used to group cells
#' @param group_order Optional: Explicit specification of group order by supplying a vector of group_ids.
#' @param log_scale Logical , determines if data is log scaled before plotting. Default = FALSE.
#' @param font_size numeric object, the font size (in pts) used to make the plot.
#' @param label_height numeric object, Percent of the plot height that should be used for the labels (0 to 100). Default is 25.
#' @param show_counts Logical, whether or not to display sample counts at the top of labels. Default = TRUE.
#' @param rotate_counts Logical, whether or not to rotate sample counts by 90 degrees. Default = FALSE.
#' @param max_width numeric object, percent of plot width that should be used for maximum expression values (0 to 100). Default is 10.
#' @param return_type What values to return - can be "plot", "data", or "both". Default is "plot".
#' 
#' @return a ggplot2 plot object
#'
group_violin_plot <- function(data,
                              anno,
                              genes,
                              grouping,
                              group_order = NULL,
                              log_scale = TRUE,
                              font_size = 7, 
                              label_height = 25,
                              show_counts = TRUE, 
                              rotate_counts = FALSE,
                              max_width = 10,
                              return_type = "plot") {
  
  # Reverse so that genes go from top to bottom
  # instead of bottom to top.
  genes <- rev(genes)
  
  group_cols <- group_columns(grouping)
  
  # Filter data to genes and samples in anno
  gene_data <- filter_gene_data(data, 
                                genes, 
                                anno, 
                                group_cols,
                                group_order, 
                                "sample_name")
  
  # Filter annotations if group_order is provided
  if(!is.null(group_order)) {
    anno <- anno[anno[[group_cols$id]] %in% group_order,]
  }
  
  # Get maximum values for each gene before rescaling to plot space.
  max_vals_unscaled <- max_gene_vals(gene_data, genes)
  
  # Left-join data to anno. This will ensure that data is filtered for the cells provided in anno
  plot_data <- dplyr::left_join(anno, gene_data, by = "sample_name")
  
  max_vals_scaled <- max_vals_unscaled
  
  # Rescale values if needed
  if(log_scale) {
    plot_data <- scale_gene_data(plot_data, genes, scale_type = "log10")
    max_vals_scaled <- log10(max_vals_unscaled + 1)
  }
  
  ## Arrange
  
  # Add x-positions for each group
  plot_data <- add_group_xpos(plot_data,
                              group_cols = group_cols,
                              group_order = group_order)
  
  # Compute basic count stats that are used downstream
  # n_stats$genes, n_stats$samples, and n_stats$groups
  n_stats <- get_n_stats(plot_data, group_cols, genes)
  
  # Scale the data between i and i + 0.9
  for (i in 1:length(genes)) {
    gene <- genes[i]
    gene_max <- max_vals_scaled[gene]
    
    plot_data[[gene]] <- scale_values_plot_space(plot_data[[gene]],
                                                 min_ps = i)
    
    
  }
  
  header_labels <- build_header_labels(data = plot_data, 
                                       grouping = grouping,
                                       group_order = group_order,
                                       ymin = n_stats$genes + 1, 
                                       label_height = label_height, 
                                       label_type = "simple")
  
  label_y_size <- max(header_labels$ymax) - min(header_labels$ymin)
  
  group_data <- plot_data %>%
    dplyr::group_by(xpos) %>%
    dplyr::summarise(group_n = dplyr::n()) %>%
    dplyr::mutate(label_y = n_stats$genes + label_y_size * 0.05,
                  group_n_y = max(header_labels$ymax) - 0.1 * label_y_size)
  
  # Plot setup
  p <- ggplot2::ggplot() +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_y_continuous("", 
                                breaks = 1:n_stats$genes + 0.45, 
                                labels = genes, 
                                expand = c(0, 0)) +
    ggplot2::scale_x_continuous("", 
                                expand = c(0, 0)) +
    ggplot2::theme_classic(font_size) +
    ggplot2::theme(axis.text = ggplot2::element_text(size = ggplot2::rel(1), face = "italic"),
                   axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   legend.position = "none") +
    ggplot2::geom_hline(ggplot2::aes(yintercept = 1:(n_stats$genes)), size = 0.2)
  
  # plot the violins for each gene
  for (i in 1:length(genes)) {
    gene <- genes[[i]]
    
    # Check for lack of variance. If no variance, we only plot the median value
    # instead of a violin.
    
    if(sum(plot_data[[gene]] == 0) == nrow(plot_data)) {
      has_variance <- FALSE
    } else if(var(plot_data[[gene]]) == 0) {
      has_variance <- FALSE
    } else {
      has_variance <- TRUE
    }
    
    parsed_gene <- rlang::parse_expr(genes[i])
    parsed_color <- rlang::parse_expr(group_cols$color)
    
    if(has_variance) {
      p <- p + 
        ggplot2::geom_violin(data = plot_data,
                             ggplot2::aes(x = xpos, 
                                          y = !!parsed_gene, 
                                          fill = !!parsed_color),
                             scale = "width",
                             adjust = 2,
                             linewidth = 0.1)
    }
    p <- p +
      ggplot2::stat_summary(data = plot_data,
                            ggplot2::aes(x = xpos, 
                                         y = !!parsed_gene),
                            fun = "median", 
                            fun.min = "median", 
                            fun.max = "median", 
                            geom = "point", 
                            size = 0.7)
  }
  
  # Cluster labels
  p <- ggplot_header_labels(p,
                            header_labels = header_labels,
                            header_polygons = NULL,
                            font_size = font_size)
  
  ### Maximum value labels at the right edge of the plot
  max_val_dfs <- build_max_dfs(n_stats, 
                               width_stat = "groups", 
                               max_vals_unscaled, 
                               max_width)
  
  p <- ggplot_max_vals(p,
                       n_stats = n_stats,
                       width_stat = "groups",
                       max_val_dfs = max_val_dfs,
                       font_size = font_size)
  
  # Cluster counts
  if(show_counts) {
    if(rotate_counts) {
      p <- p + 
        ggplot2::geom_text(data = group_data,
                           ggplot2::aes(x = xpos,
                                        y = group_n_y, 
                                        label = group_n),
                           angle = 90,
                           hjust = 1, vjust = 0.35, 
                           size = pt2mm(font_size))
      
    } else {
      p <- p + 
        ggplot2::geom_text(data = group_data,
                           ggplot2::aes(x = xpos,
                                        y = group_n_y, 
                                        label = group_n),
                           size = pt2mm(font_size))
    }
  }
  
  if(return_type == "plot") {
    return(p)
  } else if(return_type == "data") {
    return(list(plot_data = plot_data,
                header_labels = header_labels,
                header_polygons = header_polygons,
                max_val_dfs = max_val_dfs,
                n_stats = n_stats))
  } else if(return_type == "both") {
    return(list(plot = p,
                plot_data = plot_data,
                header_labels = header_labels,
                group_counts = group_data,
                max_val_dfs = max_val_dfs,
                n_stats = n_stats))
  }
}


#' Quasirandom jittered plots of gene expression for grouped samples
#' 
#' @param data A data frame containing gene expression values. The first column should be sample_name
#' @param anno Sample annotations. The first column should be sample_name, and each annotation should have \_id, \_label, and \_color columns
#' @param genes A character vector containing gene symbols to be plotted. 
#' @param grouping A character vector specifying the desc base that should be used to group cells
#' @param group_order Optional: Explicit specification of group order by supplying a vector of group_ids.
#' @param log_scale Logical , determines if data is log scaled before plotting. Default = FALSE.
#' @param font_size numeric object, the font size (in pts) used to make the plot.
#' @param label_height numeric object, Percent of the plot height that should be used for the labels (0 to 100). Default is 25.
#' @param show_counts Logical, whether or not to display sample counts at the top of labels. Default = TRUE.
#' @param rotate_counts Logical, whether or not to rotate sample counts by 90 degrees. Default = FALSE.
#' @param max_width numeric object, percent of plot width that should be used for maximum expression values (0 to 100). Default is 10.
#' @param return_type What values to return - can be "plot", "data", or "both". Default is "plot".
#' 
#' @return a ggplot2 plot object
#'
group_quasirandom_plot <- function(data,
                                   anno,
                                   genes,
                                   grouping,
                                   group_order = NULL,
                                   log_scale = TRUE,
                                   font_size = 7, 
                                   label_height = 25,
                                   show_counts = TRUE, 
                                   rotate_counts = FALSE,
                                   max_width = 10,
                                   return_type = "plot") {
  
  # Reverse so that genes go from top to bottom
  # instead of bottom to top.
  genes <- rev(genes)
  
  group_cols <- group_columns(grouping)
  
  # Filter data to genes and samples in anno
  gene_data <- filter_gene_data(data, 
                                genes, 
                                anno, 
                                group_cols,
                                group_order, 
                                "sample_name")
  
  # Filter annotations if group_order is provided
  if(!is.null(group_order)) {
    anno <- anno[anno[[group_cols$id]] %in% group_order,]
  }
  
  # Get maximum values for each gene before rescaling to plot space.
  max_vals_unscaled <- max_gene_vals(gene_data, genes)
  
  # Left-join data to anno. This will ensure that data is filtered for the cells provided in anno
  plot_data <- dplyr::left_join(anno, gene_data, by = "sample_name")
  
  # Rescale values if needed
  if(log_scale) {
    plot_data <- scale_gene_data(plot_data, genes, scale_type = "log10")
  }
  
  # Add x-positions for each group
  plot_data <- add_group_xpos(plot_data,
                              group_cols = group_cols,
                              group_order = group_order)
  
  # Compute basic count stats that are used downstream
  # n_stats$genes, n_stats$samples, and n_stats$groups
  n_stats <- get_n_stats(plot_data, group_cols, genes)
  
  # Scale the data between i and i + 0.9
  for (i in 1:length(genes)) {
    gene <- genes[i]
    
    plot_data[[gene]] <- scale_values_plot_space(plot_data[[gene]],
                                                 min_ps = i)
  }
  
  header_labels <- build_header_labels(data = plot_data, 
                                       grouping = grouping,
                                       group_order = group_order,
                                       ymin = n_stats$genes + 1, 
                                       label_height = label_height, 
                                       label_type = "simple")
  
  label_y_size <- max(header_labels$ymax) - min(header_labels$ymin)
  
  group_data <- plot_data %>%
    dplyr::group_by(xpos) %>%
    dplyr::summarise(group_n = dplyr::n()) %>%
    dplyr::mutate(label_y = n_stats$genes + label_y_size * 0.05,
                  group_n_y = max(header_labels$ymax) - 0.1 * label_y_size)
  
  # Plot setup
  p <-ggplot2::ggplot() +
    ggplot2::scale_color_identity() +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_y_continuous("", 
                                breaks = 1:length(genes) + 0.45, 
                                labels = genes, 
                                expand = c(0, 0)) +
    ggplot2::scale_x_continuous("", 
                                expand = c(0, 0)) +
    ggplot2::theme_classic(font_size) +
    theme(axis.text = ggplot2::element_text(size = ggplot2::rel(1), face = "italic"),
          axis.text.x = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank(),
          legend.position = "none") +
    ggplot2::geom_hline(aes(yintercept = 1:(n_stats$genes)), size = 0.2)
  
  # plot the swarms for each gene
  for (i in 1:n_stats$genes) {
    gene <- genes[[i]]
    if (var(plot_data[[gene]]) != 0) {
      parsed_gene <- rlang::parse_expr(genes[i])
      parsed_color <- rlang::parse_expr(group_cols$color)
      p <- p + 
        ggbeeswarm::geom_quasirandom(data = plot_data,
                                     ggplot2::aes(x = xpos, 
                                                  y = !!parsed_gene, 
                                                  color = !!parsed_color),
                                     size = 0.2,
                                     groupOnX = TRUE)
    }
    
  }
  
  # Cluster labels
  p <- ggplot_header_labels(p,
                            header_labels = header_labels,
                            header_polygons = NULL,
                            font_size = font_size)
  
  ### Maximum value labels at the right edge of the plot
  max_val_dfs <- build_max_dfs(n_stats, 
                               width_stat = "groups", 
                               max_vals_unscaled, 
                               max_width)
  
  p <- ggplot_max_vals(p,
                       n_stats = n_stats,
                       width_stat = "groups",
                       max_val_dfs = max_val_dfs,
                       font_size = font_size)
  
  # Cluster counts
  if (show_counts) {
    if (rotate_counts) {
      p <- p + ggplot2::geom_text(data = group_data,
                                  ggplot2::aes(x = xpos,
                                               y = group_n_y, 
                                               label = group_n),
                                  angle = 90,
                                  hjust = 1, vjust = 0.35, 
                                  size = pt2mm(font_size))
    } else {
      p <- p + ggplot2::geom_text(data = group_data,
                                  ggplot2::aes(x = xpos,
                                               y = group_n_y, 
                                               label = group_n),
                                  size = pt2mm(font_size))
    }
  }
  
  if(return_type == "plot") {
    return(p)
  } else if(return_type == "data") {
    return(list(plot_data = plot_data,
                header_labels = header_labels,
                header_polygons = header_polygons,
                max_val_dfs = max_val_dfs,
                n_stats = n_stats))
  } else if(return_type == "both") {
    return(list(plot = p,
                plot_data = plot_data,
                header_labels = header_labels,
                group_counts = group_data,
                max_val_dfs = max_val_dfs,
                n_stats = n_stats))
  }
}

#' Box Plots plots of gene expression for grouped samples
#' 
#' @param data A data frame containing gene expression values. The first column should be sample_name
#' @param anno Sample annotations. The first column should be sample_name, and each annotation should have \_id, \_label, and \_color columns
#' @param genes A character vector containing gene symbols to be plotted. 
#' @param grouping A character vector specifying the desc base that should be used to group cells
#' @param group_order Optional: Explicit specification of group order by supplying a vector of group_ids.
#' @param log_scale Logical , determines if data is log scaled before plotting. Default = FALSE.
#' @param font_size numeric object, the font size (in pts) used to make the plot.
#' @param label_height numeric object, Percent of the plot height that should be used for the labels (0 to 100). Default is 25.
#' @param show_counts Logical, whether or not to display sample counts at the top of labels. Default = TRUE.
#' @param rotate_counts Logical, whether or not to rotate sample counts by 90 degrees. Default = FALSE.
#' @param max_width numeric object, percent of plot width that should be used for maximum expression values (0 to 100). Default is 10.
#' @param return_type What values to return - can be "plot", "data", or "both". Default is "plot".
#' 
#' @return a ggplot2 plot object
#'
group_box_plot <- function(data,
                           anno,
                           genes,
                           grouping,
                           group_order = NULL,
                           log_scale = TRUE,
                           font_size = 7, 
                           label_height = 25,
                           show_counts = TRUE, 
                           rotate_counts = FALSE,
                           max_width = 10,
                           return_type = "plot") {
  # Reverse so that genes go from top to bottom
  # instead of bottom to top.
  genes <- rev(genes)
  
  group_cols <- group_columns(grouping)
  
  # Filter data to genes and samples in anno
  gene_data <- filter_gene_data(data, 
                                genes, 
                                anno, 
                                group_cols,
                                group_order, 
                                "sample_name")
  
  # Filter annotations if group_order is provided
  if(!is.null(group_order)) {
    anno <- anno[anno[[group_cols$id]] %in% group_order,]
  }
  
  # Get maximum values for each gene before rescaling to plot space.
  max_vals_unscaled <- max_gene_vals(gene_data, genes)
  
  # Left-join data to anno. This will ensure that data is filtered for the cells provided in anno
  plot_data <- dplyr::left_join(anno, gene_data, by = "sample_name")
  
  # Rescale values if needed
  if(log_scale) {
    plot_data <- scale_gene_data(plot_data, genes, scale_type = "log10")
  }
  
  # Add x-positions for each group
  plot_data <- add_group_xpos(plot_data,
                              group_cols = group_cols,
                              group_order = group_order)
  
  # Compute basic count stats that are used downstream
  # n_stats$genes, n_stats$samples, and n_stats$groups
  n_stats <- get_n_stats(plot_data, group_cols, genes)
  
  # Scale the data between i and i + 0.9
  for (i in 1:length(genes)) {
    gene <- genes[i]
    
    plot_data[[gene]] <- scale_values_plot_space(plot_data[[gene]],
                                                 min_ps = i)
  }
  
  header_labels <- build_header_labels(data = plot_data, 
                                       grouping = grouping,
                                       group_order = group_order,
                                       ymin = n_stats$genes + 1, 
                                       label_height = label_height, 
                                       label_type = "simple")
  
  label_y_size <- max(header_labels$ymax) - min(header_labels$ymin)
  
  group_data <- plot_data %>%
    dplyr::group_by(xpos) %>%
    dplyr::summarise(group_n = n()) %>%
    dplyr::mutate(label_y = n_stats$genes + label_y_size * 0.05,
                  group_n_y = max(header_labels$ymax) - 0.1 * label_y_size)
  
  # Plot setup
  p <-ggplot2::ggplot() +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_y_continuous("", 
                                breaks = 1:length(genes) + 0.45, 
                                labels = genes, 
                                expand = c(0, 0)) +
    ggplot2::scale_x_continuous("", 
                                expand = c(0, 0)) +
    ggplot2::theme_classic(font_size) +
    theme(axis.text = ggplot2::element_text(size = ggplot2::rel(1), face = "italic"),
          axis.text.x = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank(),
          legend.position = "none") +
    ggplot2::geom_hline(aes(yintercept = 1:(n_stats$genes)), size = 0.2)
  
  # plot the boxplots for each gene
  for (i in 1:length(genes)) {
    gene <- genes[[i]]
    if (var(plot_data[[gene]]) != 0) {
      parsed_gene <- rlang::parse_expr(genes[i])
      parsed_color <- rlang::parse_expr(group_cols$color)
      p <- p + 
        ggplot2::geom_boxplot(data = plot_data,
                              ggplot2::aes(x = "xpos", 
                                           y = !!parsed_gene, 
                                           fill = !!parsed_color),
                              width = 0.8,
                              outlier.size = 0.2,
                              size = 0.1)
    }
    
  }
  
  # Cluster labels
  p <- ggplot_header_labels(p,
                            header_labels = header_labels,
                            header_polygons = NULL,
                            font_size = font_size)
  
  ### Maximum value labels at the right edge of the plot
  max_val_dfs <- build_max_dfs(n_stats, 
                               width_stat = "groups", 
                               max_vals_unscaled, 
                               max_width)
  
  p <- ggplot_max_vals(p,
                       n_stats = n_stats,
                       width_stat = "groups",
                       max_val_dfs = max_val_dfs,
                       font_size = font_size)
  
  # Cluster counts
  if (show_counts) {
    if (rotate_counts) {
      p <- p + ggplot2::geom_text(data = group_data,
                                  ggplot2::aes(x = xpos,
                                               y = group_n_y, 
                                               label = group_n),
                                  angle = 90,
                                  hjust = 1, vjust = 0.35, 
                                  size = pt2mm(font_size))
    } else {
      p <- p + ggplot2::geom_text(data = group_data,
                                  ggplot2::aes(x = xpos,
                                               y = group_n_y, 
                                               label = group_n),
                                  size = pt2mm(font_size))
    }
  }
  
  if(return_type == "plot") {
    return(p)
  } else if(return_type == "data") {
    return(list(plot_data = plot_data,
                header_labels = header_labels,
                header_polygons = header_polygons,
                max_val_dfs = max_val_dfs,
                n_stats = n_stats))
  } else if(return_type == "both") {
    return(list(plot = p,
                plot_data = plot_data,
                header_labels = header_labels,
                group_counts = group_data,
                max_val_dfs = max_val_dfs,
                n_stats = n_stats))
  }
}


#' Heatmap plots of group summary statistics
#' 
#' @param data A data frame containing gene expression values. The first column should be sample_name
#' @param anno Sample annotations. The first column should be sample_name, and each annotation should have \_id, \_label, and \_color columns
#' @param genes A character vector containing gene symbols to be plotted. 
#' @param grouping A character vector specifying the desc base that should be used to group cells
#' @param group_order Optional: Explicit specification of group order by supplying a vector of group_ids.
#' @param stat The statistic to apply to each group. Options are: 
#' \itemize{
#'   \item "median"
#'   \item "mean"
#'   \item "tmean" (25\% trimmed mean)
#'   \item "nzmean" (mean of non-zero values)
#'   \item "nzmedian" (median of non-zero values)
#'   \item "prop_gt0" (proportion of samples > 0)
#'   \item "prop_gt1" (proportion of samples > 1)
#'   \item "min"
#'   \item "max"
#'   }
#' @param log_scale Logical , determines if data is log scaled before plotting. Default = FALSE.
#' @param normalize_rows Logical, whether or not to rescale data within each row of the plot. Default = FALSE.
#' @param font_size numeric object, the font size (in pts) used to make the plot.
#' @param label_height numeric object, Percent of the plot height that should be used for the labels (0 to 100). Default is 25.
#' @param show_counts Logical, whether or not to display sample counts at the top of labels. Default = TRUE.
#' @param rotate_counts Logical, whether or not to rotate sample counts by 90 degrees. Default = FALSE.
#' @param max_width numeric object, percent of plot width that should be used for maximum expression values (0 to 100). Default is 10.
#' @param return_type What values to return - can be "plot", "data", or "both". Default is "plot".
#' 
#' @return a ggplot2 plot object
#'
group_heatmap_plot <- function(data,
                               anno,
                               genes,
                               grouping,
                               group_order = NULL,
                               stat = "median",
                               log_scale = TRUE,
                               normalize_rows = FALSE,
                               colorset = NULL,
                               font_size = 7, 
                               label_height = 25,
                               show_counts = TRUE, 
                               rotate_counts = FALSE,
                               max_width = 10,
                               return_type = "plot") {
  
  # Reverse so that genes go from top to bottom
  # instead of bottom to top.
  genes <- rev(genes)
  
  group_cols <- group_columns(grouping)
  
  # Filter data to genes and samples in anno
  gene_data <- filter_gene_data(data, 
                                genes, 
                                anno, 
                                group_cols,
                                group_order, 
                                "sample_name")
  
  # Filter annotations if group_order is provided
  if(!is.null(group_order)) {
    anno <- anno[anno[[group_cols$id]] %in% group_order,]
  }
  
  gene_stats <- group_stats(gene_data,
                            value_cols = genes,
                            anno = anno,
                            grouping = group_cols$label,
                            stat = stat)
  
  # Get maximum values for each gene before rescaling to plot space.
  max_vals_unscaled <- max_gene_vals(gene_stats, genes)
  
  if (log_scale) {
    gene_stats <- scale_gene_data(gene_stats, genes, scale_type = "log10")
  }
  
  # Convert the data values to heatmap colors
  gene_stats <- data_df_to_colors(gene_stats,
                                  value_cols = genes,
                                  per_col = normalize_rows,
                                  colorset = colorset)
  
  # Left-join data to anno. This will ensure that data is filtered for the cells provided in anno
  plot_anno <- anno %>%
    dplyr::select(dplyr::one_of(group_cols$id, group_cols$label, group_cols$color)) %>%
    unique()
  
  parsed_id <- rlang::parse_expr(group_cols$id)
  group_counts <- anno %>%
    dplyr::group_by(!!parsed_id) %>%
    dplyr::summarise(group_n = n())
  
  plot_data <- dplyr::left_join(plot_anno, gene_stats, by = group_cols$label)
  plot_data <- dplyr::left_join(plot_data, group_counts, by = group_cols$id)
  
  # Add x-positions for each group
  plot_data <- add_group_xpos(plot_data,
                              group_cols = group_cols,
                              group_order = group_order)
  
  # Compute basic count stats that are used downstream
  # n_stats$genes, n_stats$samples, and n_stats$groups
  n_stats <- get_n_stats(plot_data, group_cols, genes)
  
  header_labels <- build_header_labels(data = plot_data, 
                                       grouping = grouping,
                                       group_order = group_order,
                                       ymin = n_stats$genes + 1, 
                                       label_height = label_height, 
                                       label_type = "simple")
  
  label_y_size <- max(header_labels$ymax) - min(header_labels$ymin)
  
  group_data <- plot_data %>%
    dplyr::select(xpos, group_n) %>%
    dplyr::mutate(label_y = n_stats$genes + label_y_size * 0.05,
                  group_n_y = max(header_labels$ymax) - 0.1 * label_y_size)
  
  # Plot setup
  p <-ggplot2::ggplot() +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_y_continuous("", 
                                breaks = 1:length(genes) + 0.45, 
                                labels = genes, 
                                expand = c(0, 0)) +
    ggplot2::scale_x_continuous("", 
                                expand = c(0, 0)) +
    ggplot2::theme_classic(font_size) +
    theme(axis.text = ggplot2::element_text(size = ggplot2::rel(1), face = "italic"),
          axis.text.x = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank(),
          legend.position = "none") +
    ggplot2::geom_hline(aes(yintercept = 1:(n_stats$genes)), size = 0.2)
  
  # plot the heatmap for each gene
  for (i in 1:length(genes)) {
    gene <- genes[[i]]
    parsed_gene <- rlang::parse_expr(gene)
    parsed_y <- rlang::parse_expr(i)
    
    p <- p + 
      ggplot2::geom_rect(data = plot_data,
                         ggplot2::aes(xmin = xpos - 0.5, 
                                      xmax = xpos + 0.5,
                                      ymin = !!parsed_y, 
                                      ymax = !!parsed_ + 1, 
                                      fill = !!parsed_gene))
    
    
  }
  
  # Cluster labels
  p <- ggplot_header_labels(p,
                            header_labels = header_labels,
                            header_polygons = NULL,
                            font_size = font_size)
  
  ### Maximum value labels at the right edge of the plot
  max_val_dfs <- build_max_dfs(n_stats, 
                               width_stat = "groups", 
                               max_vals_unscaled, 
                               max_width)
  
  p <- ggplot_max_vals(p,
                       n_stats = n_stats,
                       width_stat = "groups",
                       max_val_dfs = max_val_dfs,
                       font_size = font_size)
  
  # Cluster counts
  if (show_counts) {
    if (rotate_counts) {
      p <- p + ggplot2::geom_text(data = group_data,
                                  ggplot2::aes(x = xpos,
                                               y = group_n_y, 
                                               label = group_n),
                                  angle = 90,
                                  hjust = 1, vjust = 0.35, 
                                  size = pt2mm(font_size))
    } else {
      p <- p + ggplot2::geom_text(data = group_data,
                                  ggplot2::aes(x = xpos,
                                               y = group_n_y, 
                                               label = group_n),
                                  size = pt2mm(font_size))
    }
  }
  
  if(return_type == "plot") {
    return(p)
  } else if(return_type == "data") {
    return(list(plot_data = plot_data,
                header_labels = header_labels,
                header_polygons = header_polygons,
                max_val_dfs = max_val_dfs,
                n_stats = n_stats))
  } else if(return_type == "both") {
    return(list(plot = p,
                plot_data = plot_data,
                header_labels = header_labels,
                group_counts = group_data,
                max_val_dfs = max_val_dfs,
                n_stats = n_stats))
  }
}

#' Dot-plot Heatmap plots of group summary statistics
#' 
#' @param data A data frame containing gene expression values. The first column should be sample_name
#' @param anno Sample annotations. The first column should be sample_name, and each annotation should have \_id, \_label, and \_color columns
#' @param genes A character vector containing gene symbols to be plotted. 
#' @param grouping A character vector specifying the desc base that should be used to group cells
#' @param group_order Optional: Explicit specification of group order by supplying a vector of group_ids.
#' @param fill_stat The statistic to apply to each group for use as dot fill color. Default = "median". Options are: 
#' \itemize{
#'   \item "median"
#'   \item "mean"
#'   \item "tmean" (25\% trimmed mean)
#'   \item "nzmean" (mean of non-zero values)
#'   \item "nzmedian" (median of non-zero values)
#'   \item "prop_gt0" (proportion of samples > 0)
#'   \item "prop_gt1" (proportion of samples > 1)
#'   \item "min"
#'   \item "max"
#'   }
#' @param size_stat The statistic to apply to each group for scaling dot size. Same options as fill_stat. Default = "prop_gt0".
#' @param max_size Maximum size of dots, in pts.
#' @param log_scale Logical , determines if data is log scaled before plotting. Default = FALSE.
#' @param normalize_rows Logical, whether or not to rescale data within each row of the plot. Default = FALSE.
#' @param font_size numeric object, the font size (in pts) used to make the plot.
#' @param label_height numeric object, Percent of the plot height that should be used for the labels (0 to 100). Default is 25.
#' @param show_counts Logical, whether or not to display sample counts at the top of labels. Default = TRUE.
#' @param rotate_counts Logical, whether or not to rotate sample counts by 90 degrees. Default = FALSE.
#' @param max_width numeric object, percent of plot width that should be used for maximum expression values (0 to 100). Default is 10.
#' @param return_type What values to return - can be "plot", "data", or "both". Default is "plot".
#' 
#' @return a ggplot2 plot object
#'
group_dot_plot <- function(data,
                           anno,
                           genes,
                           grouping,
                           group_order = NULL,
                           fill_stat = "median",
                           size_stat = "prop_gt0",
                           max_size = 10,
                           log_scale = TRUE,
                           normalize_rows = FALSE,
                           colorset = NULL,
                           font_size = 7, 
                           label_height = 25,
                           show_counts = TRUE, 
                           rotate_counts = FALSE,
                           max_width = 10,
                           return_type = "plot") {
  # Reverse so that genes go from top to bottom
  # instead of bottom to top.
  genes <- rev(genes)
  
  group_cols <- group_columns(grouping)
  
  # Filter data to genes and samples in anno
  gene_data <- filter_gene_data(data, 
                                genes, 
                                anno, 
                                group_cols,
                                group_order, 
                                "sample_name")
  
  # Filter annotations if group_order is provided
  if(!is.null(group_order)) {
    anno <- anno[anno[[group_cols$id]] %in% group_order,]
  }
  
  gene_fill_stats <- group_stats(gene_data,
                                 value_cols = genes,
                                 anno = anno,
                                 grouping = group_cols$label,
                                 stat = fill_stat)
  
  # Get maximum values for each gene before rescaling to plot space.
  max_vals_unscaled <- max_gene_vals(gene_fill_stats, genes)
  
  gene_size_stats <- group_stats(gene_data,
                                 value_cols = genes,
                                 anno = anno,
                                 grouping = group_cols$label,
                                 stat = size_stat)
  
  if(log_scale) {
    gene_fill_stats <- scale_gene_data(gene_fill_stats, genes, scale_type = "log10")
  }
  
  # Convert the data values to heatmap colors
  gene_fill_data <- data_df_to_colors(gene_fill_stats,
                                      value_cols = genes,
                                      per_col = normalize_rows,
                                      colorset = colorset)
  
  names(gene_fill_data)[match(genes, names(gene_fill_data))] <- paste0(genes, "_fill")
  names(gene_size_stats)[match(genes, names(gene_size_stats))] <- paste0(genes, "_size")
  
  # Left-join data to anno. This will ensure that data is filtered for the cells provided in anno
  plot_anno <- anno %>%
    dplyr::select(dplyr::one_of(group_cols$id, group_cols$label, group_cols$color)) %>%
    unique()
  
  parsed_id <- rlang::parse_expr(group_cols$id)
  group_counts <- anno %>%
    dplyr::group_by(!!parsed_id) %>%
    dplyr::summarise(group_n = dplyr::n())
  
  plot_data <- plot_anno %>%
    dplyr::left_join(gene_fill_data, by = group_cols$label) %>%
    dplyr::left_join(gene_size_stats, by = group_cols$label) %>%
    dplyr::left_join(group_counts, by = group_cols$id)
  
  # Add x-positions for each group
  plot_data <- add_group_xpos(plot_data,
                              group_cols = group_cols,
                              group_order = group_order)
  
  # Compute basic count stats that are used downstream
  # n_stats$genes, n_stats$samples, and n_stats$groups
  n_stats <- get_n_stats(plot_data, group_cols, genes)
  
  header_labels <-build_header_labels(data = plot_data, 
                                      grouping = grouping,
                                      group_order = group_order,
                                      ymin = n_stats$genes + 1, 
                                      label_height = label_height, 
                                      label_type = "simple")
  
  label_y_size <- max(header_labels$ymax) - min(header_labels$ymin)
  
  group_data <- plot_data %>%
    dplyr::select(xpos, group_n) %>%
    dplyr::mutate(label_y = n_stats$genes + label_y_size * 0.05,
                  group_n_y = max(header_labels$ymax) - 0.1 * label_y_size)
  
  # Plot setup
  p <-ggplot2::ggplot() +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_size_area(max_size = pt2mm(max_size)) +
    ggplot2::scale_y_continuous("", 
                                breaks = 1:length(genes) + 0.45, 
                                labels = genes, 
                                expand = c(0, 0)) +
    ggplot2::scale_x_continuous("", 
                                expand = c(0, 0)) +
    ggplot2::theme_classic(font_size) +
    theme(axis.text = ggplot2::element_text(size = ggplot2::rel(1), face = "italic"),
          axis.text.x = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank(),
          legend.position = "none") +
    ggplot2::geom_hline(aes(yintercept = 1:(n_stats$genes)), size = 0.2)
  
  # plot the dots for each gene
  for(i in 1:length(genes)) {
    gene <- genes[[i]]
    gene_fill <- paste0(gene, "_fill")
    parsed_gene_fill <- rlang::parse_expr(gene_fill)
    gene_size <- paste0(gene, "_size")
    parsed_gene_size <- rlang::parse_expr(gene_size)
    
    plot_data[[paste0(gene,"_y")]] <- i + 0.5
    parsed_y <- rlang::parse_expr(paste0(gene,"_y"))
    
    p <- p + 
      ggplot2::geom_point(data = plot_data,
                          ggplot2::aes(x = xpos,
                                       y = !!parsed_y, 
                                       fill = !!parsed_gene_fill,
                                       size = !!parsed_gene_size),
                          pch = 21)
    
    
  }
  
  # Cluster labels
  p <- ggplot_header_labels(p,
                            header_labels = header_labels,
                            header_polygons = NULL,
                            font_size = font_size)
  
  ### Maximum value labels at the right edge of the plot
  max_val_dfs <- build_max_dfs(n_stats, 
                               width_stat = "groups", 
                               max_vals_unscaled, 
                               max_width)
  
  p <- ggplot_max_vals(p,
                       n_stats = n_stats,
                       width_stat = "groups",
                       max_val_dfs = max_val_dfs,
                       font_size = font_size)
  
  # Cluster counts
  if (show_counts) {
    if (rotate_counts) {
      p <- p + ggplot2::geom_text(data = group_data,
                                  ggplot2::aes(x = xpos,
                                               y = group_n_y, 
                                               label = group_n),
                                  angle = 90,
                                  hjust = 1, vjust = 0.35, 
                                  size = pt2mm(font_size))
    } else {
      p <- p + ggplot2::geom_text(data = group_data,
                                  ggplot2::aes(x = xpos,
                                               y = group_n_y, 
                                               label = group_n),
                                  size = pt2mm(font_size))
    }
  }
  
  if(return_type == "plot") {
    return(p)
  } else if(return_type == "data") {
    return(list(plot_data = plot_data,
                header_labels = header_labels,
                header_polygons = header_polygons,
                max_val_dfs = max_val_dfs,
                n_stats = n_stats))
  } else if(return_type == "both") {
    return(list(plot = p,
                plot_data = plot_data,
                header_labels = header_labels,
                group_counts = group_data,
                max_val_dfs = max_val_dfs,
                n_stats = n_stats))
  }
}

#' Build a heatmap legend plot
#' 
#' @param min_val numeric, the minimum value in the plot scale (default = 0)
#' @param max_val numeric, the maximum value in the plot scale (default = 4)
#' @param ggplot2::scale_name character, the name for the values displayed (default = "FPKM")
#' @param colorset character vector, the colors to interpolate between using colorRampPalette.
#' 
#' @return a ggplot2 heatmap legend plot
heatmap_legend_plot <- function(min_val = 0, 
                                max_val = 4, 
                                scale_name = "FPKM", 
                                colorset = c("darkblue","dodgerblue","gray80","orange","orangered")) {
  
  colors <- colorRampPalette(colorset)(1001)
  
  ## Build ggplot2::geom_rect() compatible table
  legend_data <- data.frame(xmin = 1:1001,
                            xmax = 1:1001 + 1,
                            ymin = 0,
                            ymax = 1,
                            fill = colors)
  
  border_data <- data.frame(x = min(legend_data$xmin), 
                            xend = max(legend_data$xmax))
  
  legend_plot <-ggplot2::ggplot() +
    ggplot2::geom_rect(data = legend_data,
                       ggplot2::aes(xmin = xmin, xmax = xmax, 
                                    ymin = ymin, ymax = ymax, 
                                    fill = fill)) +
    ggplot2::geom_segment(
      data = border_data,
      ggplot2::aes(x = x, xend = xend, 
                   y = 0, yend = 0)) +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_y_continuous(expand = c(0,0)) +
    ggplot2::scale_x_continuous(scale_name, 
                                breaks = c(0, 250, 500, 750, 1000),
                                labels = round(seq(min_val, max_val, by = (max_val - min_val) / 4), 2)) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   axis.line.y = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank(),
                   axis.line.x = ggplot2::element_blank())
  
  return(legend_plot)
  
}

#' Build a point size legend plot
#' 
#' @param min_val numeric, the minimum value in the plot scale (default = 0)
#' @param max_val numeric, the maximum value in the plot scale (default = 1)
#' @param max_size numeric, the max size of the point scale (default = 6, default for ggplot2::scale_size_area() )
#' @param ggplot2::scale_name character, the name for the values displayed (default = "Fraction of cells")
#' @param n_sizes numeric, the number of points in the scale from min_val to max_val, inclusive (default = 6)
#' 
#' @return a ggplot2 point size legend plot
#' 
point_size_legend_plot <- function(min_val = 0,
                                   max_val = 1,
                                   max_size = 6,
                                   scale_name = "Fraction of cells",
                                   n_sizes = 6) {
  data <- data.frame(x = 1:n_sizes,
                     y = 1,
                     vals = seq(min_val, max_val, length.out = n_sizes))
  
  legend_plot <-ggplot2::ggplot() +
    ggplot2::geom_point(data = data,
                        ggplot2::aes(x = x,
                                     y = y,
                                     size = vals),
                        pch = 21,
                        fill = NA) +
    ggplot2::scale_size_area(max_size = max_size) +
    ggplot2::scale_x_continuous(scale_name,
                                breaks = data$x,
                                labels = data$vals) +
    ggplot2::scale_y_continuous("") +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none",
                   axis.text.y = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   axis.line = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank())
  
  return(legend_plot)
}