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
                            bg_color = "#ADCFE0") {
  
  library(dplyr)
  library(ggplot2)
  library(purrr)
  
  genes <- rev(genes)
  
  group_id <- paste0(grouping, "_id")
  group_label <- paste0(grouping, "_label")
  group_color <- paste0(grouping, "_color")
  
  gene_data <- data[, c("sample_name",genes)]
  gene_data <- gene_data[match(anno$sample_name, data$sample_name),]
  
  # Get maximum values for each gene before rescaling to plot space.
  max_vals <- map_dbl(genes, function(x) { max(gene_data[[x]]) })
  
  # Left-join data to anno. This will ensure that data is filtered for the cells provided in anno
  plot_data <- left_join(anno, gene_data, by = "sample_name")
  
  if(log_scale) {
    plot_data[,genes] <- log10(plot_data[,genes] + 1)
  }
  
  # Add an x position to each sample.
  if(!is.null(group_order)) {
    # Because we allow ranges, and groups may not necessarily be continuous integer sets
    # We have to filter out any that don't match first.
    group_order <- group_order[group_order %in% anno[[group_id]]]
    
    group_order_df <- data.frame(group = group_order) %>%
      mutate(.plot_order = 1:n())
    names(group_order_df)[1] <- group_id
    
    plot_data <- plot_data %>%
      filter_(paste0(group_id, " %in% group_order")) %>%
      left_join(group_order_df, by = group_id) %>%
      arrange(.plot_order) %>%
      mutate(xpos = 1:n()) %>%
      select(-.plot_order)
  } else {
    # Otherwise, arrange using the group_id for the group_by parameter, and use that order.
    plot_data <- plot_data %>%
      arrange_(group_id) %>%
      mutate(xpos = 1:n())
  }
  
  # Calculate the number of genes and samples for use as plot dimensions
  n_genes <- length(genes)
  n_groups <- length(unique(plot_data[[group_id]]))
  n_samples <- nrow(plot_data)
  

  
  # build_header_polygons from plot_components.R
  header_polygons <- build_header_polygons(data = plot_data, 
                                           grouping = grouping,
                                           #group_order = group_order,
                                           ymin = n_genes + 1, 
                                           label_height = label_height, 
                                           poly_type = label_type)
  
  # Build the cell type label rectangles from plot_components.R
  header_labels <-build_header_labels(data = plot_data, 
                                      grouping = grouping,
                                      ymin = n_genes + 1, 
                                      label_height = label_height, 
                                      label_type = label_type)
  
  # Calculate Plot Scale bars
  scale_bars <- data.frame(gene = genes,
                           ymin = 1:n_genes,
                           ymid = 1:n_genes + 0.45,
                           ymax = 1:n_genes + 0.9,
                           xmin = -nrow(plot_data) * 0.01,
                           xmax = 0)
  
  # Calculate segments for cluster sepration lines
  segment_lines <- plot_data %>%
    group_by_(group_id) %>%
    summarise(x = max(xpos) )
  
  # Build the maximum value labels for the right edge
  max_labels <- data.frame(x = n_samples * 1.01,
                           y = 1:n_genes + 0.5,
                           label = sci_label(max_vals))
  max_header <- data.frame(x = n_samples * 1.01,
                           y = n_genes + 1,
                           label = "Max value")
  max_width <- n_samples*(max_width/100)/(1-max_width/100)
  
  # The background of the plot is a rectangular object.
  background_data <- data.frame(xmin = 0, 
                                xmax = n_samples, 
                                ymin = 1, 
                                ymax = n_genes + 1, 
                                fill = bg_color)
  
  # pt2mm function is used for text labels, from plot_components.R
  
  # Plot setup
  p <- ggplot() +
    scale_fill_identity() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0), 
                       breaks = 1:n_genes + 0.45, 
                       labels = genes) +
    theme_classic(base_size = font_size) +
    theme(axis.text = element_text(size=rel(1)),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_blank()) +
    geom_rect(data = background_data,
              aes(xmin = xmin,
                  xmax = xmax,
                  ymin = ymin,
                  ymax = ymax,
                  fill = fill)) +
    geom_segment(data = segment_lines,
                 aes(x = x, 
                     xend = x, 
                     y = 1, 
                     yend = n_genes + 1),
                 size = 0.2, 
                 color = "gray60", 
                 linetype = "dashed")
  
  # plot the bars for each gene
  for(i in 1:length(genes)) {
    # if sort is true, arrange the values within each row from high to low.
    # if(sort) {
    #   plot_data <- plot_data %>% 
    #     arrange_(group_id, paste0("-",genes[i])) %>% 
    #     mutate(xpos = 1:n())
    # }
    # 
    plot_data[[genes[i]]] <- i + plot_data[[genes[i]]]/max(plot_data[[genes[i]]]) * 0.9
    
    # plot the rectangles for the barplots
    p <- p + 
      geom_rect(data = plot_data,
                aes_string(xmin = "xpos - 1",
                           xmax = "xpos",
                           ymin = i,
                           ymax = genes[i],
                           fill = group_color))
    
  }
  
  p <- p + 
    # Cluster labels at the top of the plot
    geom_rect(data = header_labels,
              aes(xmin = xmin , 
                  xmax = xmax, 
                  ymin = ymin, 
                  ymax = ymax, 
                  fill = color) ) +
    geom_text(data = header_labels, 
              aes(x = (xmin + xmax) / 2, 
                  y = ymin + 0.05, 
                  label = label),
              angle = 90, 
              vjust = 0.35, 
              hjust = 0, 
              size = pt2mm(font_size)) +
    geom_polygon(data = header_polygons,
                 aes(x = poly.x, 
                     y = poly.y, 
                     fill = color, 
                     group = id) ) +
    # Scale bar elements
    geom_hline(data = scale_bars,
               aes(yintercept = ymin), 
               size = 0.2) +
    geom_segment(data = scale_bars, 
                 aes(x = xmin,
                     xend = xmax,
                     y = ymid, 
                     yend = ymid),
                 size = 0.2) +
    geom_segment(data = scale_bars,
                 aes(x = xmin, 
                     xend = xmax, 
                     y = ymax, 
                     yend = ymax),
                 size = 0.2) +
    geom_segment(data = scale_bars,
                 aes(x = xmax, 
                     xend = xmax, 
                     y = ymin, 
                     yend = ymax),
                 size = 0.2) +
    # Maximum value labels at the right edge of the plot
    geom_rect(aes(xmin = n_samples + 1, 
                  xmax = n_samples + max_width, 
                  ymin = 1, 
                  ymax = max(header_labels$ymax)), 
              fill = "#FFFFFF") +
    geom_text(data = max_header,
              aes(x = x, 
                  y = y, 
                  label = label),
              angle = 90, 
              hjust = 0, 
              vjust = 0.5, 
              size = pt2mm(font_size) ) +
    geom_text(data = max_labels,
              aes(x = x, 
                  y = y, 
                  label = label),
              hjust = 0, 
              vjust = 0.5, 
              size = pt2mm(font_size) , 
              parse = TRUE)
  
  return(p)
  
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
                                max_width = 10) {
  
  library(dplyr)
  library(ggplot2)
  
  genes <- rev(genes)
  
  group_id <- paste0(grouping, "_id")
  group_label <- paste0(grouping, "_label")
  group_color <- paste0(grouping, "_color")
  
  gene_data <- data[,c("sample_name",genes)]
  gene_data <- gene_data[match(anno$sample_name, data$sample_name),]
  
  # Get maximum values for each gene before rescaling to plot space.
  max_vals <- map_dbl(genes, function(x) { max(gene_data[[x]]) })
  
  if(log_scale) {
    gene_data[,genes] <- log10(gene_data[,genes] + 1)
  }
  
  # Convert the data values to heatmap colors
  gene_data <- data_df_to_colors(gene_data,
                                 value_cols = genes,
                                 per_col = normalize_rows,
                                 colorset = colorset)

  # Left-join data to anno. This will ensure that data is filtered for the cells provided in anno
  plot_data <- left_join(anno, gene_data, by = "sample_name")
  
  # Add an x position to each sample.
  if(!is.null(group_order)) {
    # Because we allow ranges, and groups may not necessarily be continuous integer sets
    # We have to filter out any that don't match first.
    group_order <- group_order[group_order %in% anno[[group_id]]]
    
    group_order_df <- data.frame(group = group_order) %>%
      mutate(.plot_order = 1:n())
    names(group_order_df)[1] <- group_id
    
    plot_data <- plot_data %>%
      filter_(paste0(group_id, " %in% group_order")) %>%
      left_join(group_order_df, by = group_id) %>%
      arrange(.plot_order) %>%
      mutate(xpos = 1:n()) %>%
      select(-.plot_order)
  } else {
    # Otherwise, arrange using the group_id for the group_by parameter, and use that order.
    plot_data <- plot_data %>%
      arrange_(group_id) %>%
      mutate(xpos = 1:n())
  }
  
  # Calculate the number of genes and samples for use as plot dimensions
  n_genes <- length(genes)
  n_groups <- length(unique(plot_data[[group_id]]))
  n_samples <- nrow(plot_data)
  
  # build_header_polygons from plot_components.R
  header_polygons <- build_header_polygons(data = plot_data, 
                                           grouping = grouping,
                                           group_order = group_order,
                                           ymin = n_genes + 1, 
                                           label_height = label_height, 
                                           poly_type = label_type)
  
  # Build the cell type label rectangles from plot_components.R
  header_labels <-build_header_labels(data = plot_data, 
                                      grouping = grouping,
                                      ymin = n_genes + 1, 
                                      label_height = label_height, 
                                      label_type = label_type)
  
  # Build the maximum value labels for the right edge
  max_labels <- data.frame(x = n_samples * 1.01,
                           y = 1:n_genes + 0.5,
                           label = sci_label(max_vals) )
  max_header <- data.frame(x = n_samples * 1.01,
                           y = n_genes + 1,
                           label = "Max value")
  max_width <- n_samples*(max_width/100)/(1-max_width/100)
  
  # Plot setup
  p <- ggplot(data) +
    scale_fill_identity() +
    theme_classic(base_size = font_size) +
    theme(axis.text = element_text(size=rel(1)),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_blank()) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0), breaks = 1:n_genes + 0.45, labels = genes)
  
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
    p <- p + geom_rect(data = rect_data,
                       aes(xmin = xmin, 
                           xmax = xmax, 
                           ymin = ymin, 
                           ymax = ymax + 1, 
                           fill = fill))
    
  }
  
  # Label elements
  # pt2mm() is in plot_components.R
  p <- p + 
    # Cluster labels at the top of the plot
    geom_rect(data = header_labels,
              aes(xmin = xmin, xmax = xmax, 
                  ymin = ymin, ymax = ymax, 
                  fill = color)) +
    geom_text(data = header_labels,
              aes(x = (xmin + xmax) / 2, 
                  y = ymin + 0.05, 
                  label = label),
              angle = 90, 
              vjust = 0.35, hjust = 0, 
              size = pt2mm(font_size)) +
    geom_polygon(data = header_polygons,
                 aes(x = poly.x, 
                     y = poly.y, 
                     fill = color, 
                     group = id)) +
    # Maximum value labels at the right edge of the plot
    geom_rect(aes(xmin = n_samples + 1, xmax = n_samples + max_width, 
                  ymin = 1, ymax = max(header_labels$ymax)), 
              fill = "#FFFFFF") +
    geom_text(data = max_header,
              aes(x = x, y = y, 
                  label = label),
              angle = 90, 
              hjust = 0, vjust = 0.5, 
              size = pt2mm(font_size)) +
    geom_text(data = max_labels,
              aes(x = x, y = y, 
                  label = label),
              hjust = 0, vjust = 0.5, 
              size = pt2mm(font_size), parse = TRUE)
  
  p
  
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
                             max_width = 10) {
  
  library(dplyr)
  library(ggplot2)
  
  genes <- rev(genes)
  
  group_id <- paste0(grouping, "_id")
  group_label <- paste0(grouping, "_label")
  group_color <- paste0(grouping, "_color")
  
  gene_data <- data[,c("sample_name",genes)]
  gene_data <- gene_data[match(anno$sample_name, data$sample_name),]
  
  # Get maximum values for each gene before rescaling to plot space.
  max_vals <- map_dbl(genes, function(x) { max(gene_data[[x]]) })
  
  if (log_scale) {
    gene_data[,genes] <- log10(gene_data[,genes] + 1)
  }
  
  # Left-join data to anno. This will ensure that data is filtered for the cells provided in anno
  plot_data <- left_join(anno, gene_data, by = "sample_name")
  
  # Add an x position to each group
  if (!is.null(group_order)) {
    # Because we allow ranges, and groups may not necessarily be continuous integer sets
    # We have to filter out any that don't match first.
    group_order <- group_order[group_order %in% anno[[group_id]]]
    
    group_order_df <- data.frame(group = group_order) %>%
      mutate(xpos = 1:n())
    names(group_order_df)[1] <- group_id
    
    plot_data <- plot_data %>%
      filter_(paste0(group_id, " %in% group_order")) %>%
      left_join(group_order_df, by = group_id)
    
  } else {
    # Otherwise, arrange using the group_id for the group_by parameter, and use that order.
    group_order_df <- plot_data %>%
      select(one_of(group_id)) %>%
      unique() %>%
      arrange_(group_id) %>%
      mutate(xpos = 1:n())
    
    plot_data <- plot_data %>%
      left_join(group_order_df, by = group_id)
  }
  
  # Calculate the number of genes and samples for use as plot dimensions
  n_genes <- length(genes)
  n_groups <- length(unique(plot_data[[group_id]]))
  n_samples <- nrow(plot_data)
  
  # Build the cell type label rectangles from plot_components.R
  header_labels <- build_header_labels(data = plot_data, 
                                      grouping = grouping,
                                      group_order = group_order,
                                      ymin = n_genes + 1, 
                                      label_height = label_height, 
                                      label_type = "simple")
  
  # Build the maximum value labels for the right edge
  max_labels <- data.frame(x = n_groups * 1.01,
                           y = 1:n_genes + 0.5,
                           label = sci_label(max_vals) )
  max_header <- data.frame(x = n_groups * 1.01,
                           y = n_genes + 1,
                           label = "Max value")
  max_width <- n_groups*(max_width/100)/(1 - max_width/100)
  
  # Plot setup
  p <- ggplot(data) +
    scale_fill_identity() +
    theme_classic(base_size = font_size) +
    theme(axis.text = element_text(size = rel(1)),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_blank()) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0), breaks = 1:n_genes + 0.45, labels = genes)
  
  # plot the rectangles for each gene
  for (i in seq_along(genes)) {
    
    gene <- genes[i]
    
    gene_values <- plot_data %>%
      select(one_of("sample_name", gene))
    
    gene_colors <- data_df_to_colors(gene_values,
                                     value_cols = gene,
                                     per_col = normalize_rows,
                                     colorset = colorset)
    
    names(gene_colors)[names(gene_colors) == gene] <- "plot_fill"
    
    if (top_values == "highest") {
      rect_data <- plot_data %>%
        left_join(gene_colors, by = "sample_name") %>%
        arrange_(gene) %>%
        group_by(xpos) %>%
        mutate(group_n = n()) %>%
        mutate(xmin = xpos - 0.5,
               xmax = xpos + 0.5,
               ymin = seq(i, i + 1*((group_n[1] + 1)/group_n[1]), length.out = group_n[1] + 1)[-group_n[1]],
               ymax = seq(i, i + 1*((group_n[1] + 1)/group_n[1]), length.out = group_n[1] + 1)[-1])
    } else {
      arr_gene <- paste0("-",gene)
      rect_data <- plot_data %>%
        left_join(gene_colors, by = "sample_name") %>%
        arrange_(arr_gene) %>%
        group_by(xpos) %>%
        mutate(group_n = n()) %>%
        mutate(xmin = xpos - 0.5,
               xmax = xpos + 0.5,
               ymin = seq(i, i + 1*((group_n[1] + 1)/group_n[1]), length.out = group_n[1] + 1)[-group_n[1]],
               ymax = seq(i, i + 1*((group_n[1] + 1)/group_n[1]), length.out = group_n[1] + 1)[-1])
    }

             
    # plot the rectangles for the heatmap
    p <- p + geom_rect(data = rect_data,
                       aes(xmin = xmin, 
                           xmax = xmax, 
                           ymin = ymin, 
                           ymax = ymax, 
                           fill = plot_fill))
    
  }
  
  # Label elements
  # pt2mm() is in plot_components.R
  p <- p + 
    # Cluster labels at the top of the plot
    geom_rect(data = header_labels,
              aes(xmin = xmin, xmax = xmax, 
                  ymin = ymin, ymax = ymax, 
                  fill = color)) +
    geom_text(data = header_labels,
              aes(x = (xmin + xmax) / 2, 
                  y = ymin + 0.05, 
                  label = label),
              angle = 90, 
              vjust = 0.35, hjust = 0, 
              size = pt2mm(font_size)) +
    # Maximum value labels at the right edge of the plot
    geom_rect(aes(xmin = n_groups + 1, xmax = n_groups + max_width, 
                  ymin = 1, ymax = max(header_labels$ymax)), 
              fill = "#FFFFFF") +
    geom_text(data = max_header,
              aes(x = x, y = y, 
                  label = label),
              angle = 90, 
              hjust = 0, vjust = 0.5, 
              size = pt2mm(font_size)) +
    geom_text(data = max_labels,
              aes(x = x, y = y, 
                  label = label),
              hjust = 0, vjust = 0.5, 
              size = pt2mm(font_size), parse = TRUE)
  
  p
  
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
                              max_width = 10) {

  # Reverse so that genes go from top to bottom
  # instead of bottom to top.
  genes <- rev(genes)
  
  group_id <- paste0(grouping, "_id")
  group_label <- paste0(grouping, "_label")
  group_color <- paste0(grouping, "_color")
  
  ## Filter
  
  gene_data <- data[,c("sample_name",genes)]
  gene_data <- gene_data[match(anno$sample_name, data$sample_name),]
  
  # Left-join data to anno. This will ensure that data is filtered for the cells provided in anno
  plot_data <- dplyr::left_join(anno, gene_data, by = "sample_name")
  
  # If group_order is provided, filter for samples that match
  if(!is.null(group_order)) {
    plot_data <- plot_data[plot_data[[group_id]] %in% group_order,]
  }
  
  ## Max
  # Get maximum values for each gene before rescaling to plot space.
  max_vals <- purrr::map_dbl(genes, function(x) { max(plot_data[[x]]) })
  names(max_vals) <- genes
  scaled_max_vals <- max_vals
  
  ## Scale
  if(log_scale) {
    plot_data[,genes] <- log10(plot_data[,genes] + 1)
    scaled_max_vals <- log10(max_vals + 1)
  }
  
  ## Arrange
  
  # Add an x position to each group
  if (!is.null(group_order)) {
    # Because we allow ranges, and groups may not necessarily be continuous integer sets
    # We have to filter out any that don't match first.
    group_order <- group_order[group_order %in% anno[[group_id]]]
    
    group_order_df <- build_vec_pos(group_order,
                                    vec_name = group_id,
                                    axis_name = "xpos")
    
    plot_data <- plot_data %>%
      dplyr::left_join(group_order_df, by = group_id)
    
  } else {
    # Otherwise, arrange using the group_id for the group_by parameter, and use that order.
    group_order_df <- plot_data %>%
      dplyr::select(dplyr::one_of(group_id)) %>%
      unique() %>%
      dplyr::arrange_(group_id) %>%
      dplyr::mutate(xpos = 1:dplyr::n())
    
    plot_data <- plot_data %>%
      dplyr::left_join(group_order_df, by = group_id)
  }
  
  n_genes <- length(genes)
  n_groups <- length(unique(plot_data[[group_id]]))
  n_samples <- nrow(plot_data)
  
  # Scale the data between i and i + 0.9
  for (i in 1:length(genes)) {
    gene <- genes[i]
    gene_max <- scaled_max_vals[gene]
    
    plot_data[[gene]] <- i + plot_data[[gene]] / gene_max * 0.9
  }
  
  header_labels <- build_header_labels(data = plot_data, 
                                      grouping = grouping,
                                      group_order = group_order,
                                      ymin = n_genes + 1, 
                                      label_height = label_height, 
                                      label_type = "simple")
  
  # Build the maximum value labels for the right edge
  max_labels <- data.frame(x = (n_groups + 0.5) * 1.01,
                           y = 1:n_genes + 0.5,
                           label = sci_label(max_vals) )
  max_header <- data.frame(x = (n_groups + 0.5) * 1.01,
                           y = n_genes + 1,
                           label = "Max value")
  max_width <- n_groups*(max_width/100)/(1 - max_width/100)
  
  label_y_size <- max(header_labels$ymax) - min(header_labels$ymin)
  
  group_data <- plot_data %>%
    dplyr::group_by(xpos) %>%
    dplyr::summarise(group_n = dplyr::n()) %>%
    dplyr::mutate(label_y = n_genes + label_y_size * 0.05,
                  group_n_y = max(header_labels$ymax) - 0.1 * label_y_size)
  
  # Plot setup
  p <- ggplot2::ggplot() +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_y_continuous("", 
                                breaks = 1:length(genes) + 0.45, 
                                labels = genes, 
                                expand = c(0, 0)) +
    ggplot2::scale_x_continuous("", 
                                expand = c(0, 0)) +
    ggplot2::theme_classic(font_size) +
    ggplot2::theme(axis.text = ggplot2::element_text(size = ggplot2::rel(1)),
                   axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   legend.position = "none") +
    ggplot2::geom_hline(ggplot2::aes(yintercept = 1:(n_genes)), size = 0.2)
  
  # plot the violins for each gene
  for (i in 1:length(genes)) {
    gene <- genes[[i]]

    # Check for lack of variance. If no variance, we only plot the median value
    # instead of a violin.
    # Note that values of 0 will now be i due to row position scaling, above.
    if(sum(plot_data[,gene] == i) == nrow(plot_data)) {
      has_variance <- FALSE
    } else if(var(plot_data[[gene]]) == 0) {
      has_variance <- FALSE
    } else {
      has_variance <- TRUE
    }
    
    if(has_variance) {
      p <- p + 
        ggplot2::geom_violin(data = plot_data,
                             ggplot2::aes_string(x = "xpos", 
                                                 y = genes[i], 
                                                 fill = group_color),
                             scale = "width",
                             adjust = 2)
    }
    
    p <- p +
      ggplot2::stat_summary(data = plot_data,
                            ggplot2::aes_string(x = "xpos", 
                                                y = genes[i]),
                            fun.y = "median", 
                            fun.ymin = "median", 
                            fun.ymax = "median", 
                            geom = "point", 
                            size = 0.7)
  }
  
  # Cluster labels
  p <- p +
    ggplot2::geom_rect(data = header_labels, 
                       ggplot2::aes(xmin = xmin, ymin = ymin, 
                                    xmax = xmax, ymax = ymax, 
                                    fill = color)) +
    ggplot2::geom_text(data = header_labels,
                       ggplot2::aes(x = (xmin + xmax) / 2, 
                                    y = ymin + 0.05, 
                                    label = label),
                       angle = 90, 
                       vjust = 0.35, hjust = 0, 
                       size = pt2mm(font_size)) +
    # Maximum value labels on right side of plot
    ggplot2::geom_rect(ggplot2::aes(xmin = n_groups + 0.5, xmax = n_groups + 0.5 + max_width, 
                                    ymin = 1, ymax = max(header_labels$ymax)),
                       fill = "white") +
    ggplot2::geom_text(data = max_header,
                       ggplot2::aes(x = x, y = y, 
                                    label = label),
                       angle = 90, 
                       hjust = 0, vjust = 0.35, 
                       size = pt2mm(font_size)) +
    ggplot2::geom_text(data = max_labels,
                       ggplot2::aes(x = x, y = y, 
                                    label = label),
                       hjust = 0, vjust = 0.35, 
                       size = pt2mm(font_size), 
                       parse = TRUE)
  
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
  
  return(p)
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
                                   max_width = 10) {
  library(dplyr)
  library(ggplot2)
  library(ggbeeswarm)
  
  genes <- rev(genes)
  
  group_id <- paste0(grouping, "_id")
  group_label <- paste0(grouping, "_label")
  group_color <- paste0(grouping, "_color")
  
  gene_data <- data[,c("sample_name",genes)]
  gene_data <- gene_data[match(anno$sample_name, data$sample_name),]
  
  # Get maximum values for each gene before rescaling to plot space.
  max_vals <- map_dbl(genes, function(x) { max(gene_data[[x]]) })
  names(max_vals) <- genes
  
  if (log_scale) {
    gene_data[,genes] <- log10(gene_data[,genes] + 1)
  }
  
  # Left-join data to anno. This will ensure that data is filtered for the cells provided in anno
  plot_data <- left_join(anno, gene_data, by = "sample_name")
  
  # Add an x position to each group
  if (!is.null(group_order)) {
    # Because we allow ranges, and groups may not necessarily be continuous integer sets
    # We have to filter out any that don't match first.
    group_order <- group_order[group_order %in% anno[[group_id]]]
    
    group_order_df <- data.frame(group = group_order) %>%
      mutate(xpos = 1:n())
    names(group_order_df)[1] <- group_id
    
    plot_data <- plot_data %>%
      filter_(paste0(group_id, " %in% group_order")) %>%
      left_join(group_order_df, by = group_id)
    
  } else {
    # Otherwise, arrange using the group_id for the group_by parameter, and use that order.
    group_order_df <- plot_data %>%
      select(one_of(group_id)) %>%
      unique() %>%
      arrange_(group_id) %>%
      mutate(xpos = 1:n())
    
    plot_data <- plot_data %>%
      left_join(group_order_df, by = group_id)
  }
  
  n_genes <- length(genes)
  n_groups <- length(unique(plot_data[[group_id]]))
  n_samples <- nrow(data)
  
  # Scale the data between i and i + 0.9
  for (i in 1:length(genes)) {
    gene <- genes[i]
    gene_max <- max_vals[gene]
    if (log_scale) {
      plot_data[[gene]] <- i + plot_data[[gene]] / log10(gene_max + 1) * 0.9
    } else {
      plot_data[[gene]] <- i + plot_data[[gene]] / gene_max * 0.9
    }
  }
  
  header_labels <- build_header_labels(data = plot_data, 
                                      grouping = grouping,
                                      group_order = group_order,
                                      ymin = n_genes + 1, 
                                      label_height = label_height, 
                                      label_type = "simple")
  
  # Build the maximum value labels for the right edge
  max_labels <- data.frame(x = (n_groups + 0.5) * 1.01,
                           y = 1:n_genes + 0.5,
                           label = sci_label(max_vals) )
  max_header <- data.frame(x = (n_groups + 0.5) * 1.01,
                           y = n_genes + 1,
                           label = "Max value")
  max_width <- n_groups*(max_width/100)/(1 - max_width/100)
  
  label_y_size <- max(header_labels$ymax) - min(header_labels$ymin)
  
  group_data <- plot_data %>%
    group_by(xpos) %>%
    summarise(group_n = n()) %>%
    mutate(label_y = n_genes + label_y_size * 0.05,
           group_n_y = max(header_labels$ymax) - 0.1 * label_y_size)
  
  # Plot setup
  p <- ggplot() +
    scale_color_identity() +
    scale_fill_identity() +
    scale_y_continuous("", 
                       breaks = 1:length(genes) + 0.45, 
                       labels = genes, 
                       expand = c(0, 0)) +
    scale_x_continuous("", 
                       expand = c(0, 0)) +
    theme_classic(font_size) +
    theme(axis.text = element_text(size = rel(1)),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none") +
    geom_hline(aes(yintercept = 1:(n_genes)), size = 0.2)
  
  # plot the swarms for each gene
  for (i in 1:length(genes)) {
    gene <- genes[[i]]
    if (var(plot_data[[gene]]) != 0) {
      p <- p + 
        geom_quasirandom(data = plot_data,
                         aes_string(x = "xpos", 
                                    y = genes[i], 
                                    color = group_color),
                         size = 0.5,
                         groupOnX = TRUE)
    }
    
  }
  
  # Cluster labels
  p <- p +
    geom_rect(data = header_labels, 
              aes(xmin = xmin, ymin = ymin, 
                  xmax = xmax, ymax = ymax, 
                  fill = color)) +
    geom_text(data = header_labels,
              aes(x = (xmin + xmax) / 2, 
                  y = ymin + 0.05, 
                  label = label),
              angle = 90, 
              vjust = 0.35, hjust = 0, 
              size = pt2mm(font_size)) +
    # Maximum value labels on right side of plot
    geom_rect(aes(xmin = n_groups + 0.5, xmax = n_groups + 0.5 + max_width, 
                  ymin = 1, ymax = max(header_labels$ymax)),
              fill = "white") +
    geom_text(data = max_header,
              aes(x = x, y = y, 
                  label = label),
              angle = 90, 
              hjust = 0, vjust = 0.35, 
              size = pt2mm(font_size)) +
    geom_text(data = max_labels,
              aes(x = x, y = y, 
                  label = label),
              hjust = 0, vjust = 0.35, 
              size = pt2mm(font_size), 
              parse = TRUE)
  
  # Cluster counts
  if (show_counts) {
    if (rotate_counts) {
      p <- p + geom_text(data = group_data,
                         aes(x = xpos,
                             y = group_n_y, 
                             label = group_n),
                         angle = 90,
                         hjust = 1, vjust = 0.35, 
                         size = pt2mm(font_size))
    } else {
      p <- p + geom_text(data = group_data,
                         aes(x = xpos,
                             y = group_n_y, 
                             label = group_n),
                         size = pt2mm(font_size))
    }
  }
  
  return(p)
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
                           max_width = 10) {
  library(dplyr)
  library(ggplot2)
  
  genes <- rev(genes)
  
  group_id <- paste0(grouping, "_id")
  group_label <- paste0(grouping, "_label")
  group_color <- paste0(grouping, "_color")
  
  gene_data <- data[,c("sample_name",genes)]
  gene_data <- gene_data[match(anno$sample_name, data$sample_name),]
  
  # Get maximum values for each gene before rescaling to plot space.
  max_vals <- map_dbl(genes, function(x) { max(gene_data[[x]]) })
  names(max_vals) <- genes
  
  if (log_scale) {
    gene_data[,genes] <- log10(gene_data[,genes] + 1)
  }
  
  # Left-join data to anno. This will ensure that data is filtered for the cells provided in anno
  plot_data <- left_join(anno, gene_data, by = "sample_name")
  
  # Add an x position to each group
  if (!is.null(group_order)) {
    # Because we allow ranges, and groups may not necessarily be continuous integer sets
    # We have to filter out any that don't match first.
    group_order <- group_order[group_order %in% anno[[group_id]]]
    
    group_order_df <- data.frame(group = group_order) %>%
      mutate(xpos = 1:n())
    names(group_order_df)[1] <- group_id
    
    plot_data <- plot_data %>%
      filter_(paste0(group_id, " %in% group_order")) %>%
      left_join(group_order_df, by = group_id)
    
  } else {
    # Otherwise, arrange using the group_id for the group_by parameter, and use that order.
    group_order_df <- plot_data %>%
      select(one_of(group_id)) %>%
      unique() %>%
      arrange_(group_id) %>%
      mutate(xpos = 1:n())
    
    plot_data <- plot_data %>%
      left_join(group_order_df, by = group_id)
  }
  
  n_genes <- length(genes)
  n_groups <- length(unique(plot_data[[group_id]]))
  n_samples <- nrow(data)
  
  # Scale the data between i and i + 0.9
  for (i in 1:length(genes)) {
    gene <- genes[i]
    gene_max <- max_vals[gene]
    if (log_scale) {
      plot_data[[gene]] <- i + plot_data[[gene]] / log10(gene_max + 1) * 0.9
    } else {
      plot_data[[gene]] <- i + plot_data[[gene]] / gene_max * 0.9
    }
  }
  
  header_labels <- build_header_labels(data = plot_data, 
                                      grouping = grouping,
                                      group_order = group_order,
                                      ymin = n_genes + 1, 
                                      label_height = label_height, 
                                      label_type = "simple")
  
  # Build the maximum value labels for the right edge
  max_labels <- data.frame(x = (n_groups + 0.5) * 1.01,
                           y = 1:n_genes + 0.5,
                           label = sci_label(max_vals) )
  max_header <- data.frame(x = (n_groups + 0.5) * 1.01,
                           y = n_genes + 1,
                           label = "Max value")
  max_width <- n_groups*(max_width/100)/(1 - max_width/100)
  
  label_y_size <- max(header_labels$ymax) - min(header_labels$ymin)
  
  group_data <- plot_data %>%
    group_by(xpos) %>%
    summarise(group_n = n()) %>%
    mutate(label_y = n_genes + label_y_size * 0.05,
           group_n_y = max(header_labels$ymax) - 0.1 * label_y_size)
  
  # Plot setup
  p <- ggplot() +
    scale_fill_identity() +
    scale_y_continuous("", 
                       breaks = 1:length(genes) + 0.45, 
                       labels = genes, 
                       expand = c(0, 0)) +
    scale_x_continuous("", 
                       expand = c(0, 0)) +
    theme_classic(font_size) +
    theme(axis.text = element_text(size = rel(1)),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none") +
    geom_hline(aes(yintercept = 1:(n_genes)), size = 0.2)
  
  # plot the swarms for each gene
  for (i in 1:length(genes)) {
    gene <- genes[[i]]
    if (var(plot_data[[gene]]) != 0) {
      p <- p + 
        geom_boxplot(data = plot_data,
                     aes_string(x = "xpos", 
                                y = genes[i], 
                                fill = group_color),
                     width = 0.8,
                     outlier.size = 0.6)
    }
    
  }
  
  # Cluster labels
  p <- p +
    geom_rect(data = header_labels, 
              aes(xmin = xmin, ymin = ymin, 
                  xmax = xmax, ymax = ymax, 
                  fill = color)) +
    geom_text(data = header_labels,
              aes(x = (xmin + xmax) / 2, 
                  y = ymin + 0.05, 
                  label = label),
              angle = 90, 
              vjust = 0.35, hjust = 0, 
              size = pt2mm(font_size)) +
    # Maximum value labels on right side of plot
    geom_rect(aes(xmin = n_groups + 0.5, xmax = n_groups + 0.5 + max_width, 
                  ymin = 1, ymax = max(header_labels$ymax)),
              fill = "white") +
    geom_text(data = max_header,
              aes(x = x, y = y, 
                  label = label),
              angle = 90, 
              hjust = 0, vjust = 0.35, 
              size = pt2mm(font_size)) +
    geom_text(data = max_labels,
              aes(x = x, y = y, 
                  label = label),
              hjust = 0, vjust = 0.35, 
              size = pt2mm(font_size), 
              parse = TRUE)
  
  # Cluster counts
  if (show_counts) {
    if (rotate_counts) {
      p <- p + geom_text(data = group_data,
                         aes(x = xpos,
                             y = group_n_y, 
                             label = group_n),
                         angle = 90,
                         hjust = 1, vjust = 0.35, 
                         size = pt2mm(font_size))
    } else {
      p <- p + geom_text(data = group_data,
                         aes(x = xpos,
                             y = group_n_y, 
                             label = group_n),
                         size = pt2mm(font_size))
    }
  }
  
  return(p)
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
                               max_width = 10) {
  library(dplyr)
  library(ggplot2)
  
  genes <- rev(genes)
  
  group_id <- paste0(grouping, "_id")
  group_label <- paste0(grouping, "_label")
  group_color <- paste0(grouping, "_color")
  
  gene_data <- data[match(anno$sample_name, data$sample_name), c("sample_name",genes)]
  gene_data <- gene_data[match(anno$sample_name, data$sample_name),]
  
  
  gene_stats <- group_stats(gene_data,
                            value_cols = genes,
                            anno = anno,
                            grouping = group_label,
                            stat = stat)
  
  # Get maximum values for each gene before rescaling to plot space.
  max_vals <- map_dbl(genes, function(x) { max(gene_stats[[x]]) })
  names(max_vals) <- genes
  
  if (log_scale) {
    gene_stats[,genes] <- log10(gene_stats[,genes] + 1)
  }
  
  # Convert the data values to heatmap colors
  gene_stats <- data_df_to_colors(gene_stats,
                                  value_cols = genes,
                                  per_col = normalize_rows,
                                  colorset = colorset)
  
  # Left-join data to anno. This will ensure that data is filtered for the cells provided in anno
  plot_anno <- anno %>%
    select(one_of(group_id, group_label, group_color)) %>%
    unique()
  
  group_n <- anno %>%
    group_by_(group_id) %>%
    summarise(group_n = n())
  
  plot_data <- left_join(plot_anno, gene_stats, by = group_label)
  
  # Add an x position to each group
  if (!is.null(group_order)) {
    # Because we allow ranges, and groups may not necessarily be continuous integer sets
    # We have to filter out any that don't match first.
    group_order <- group_order[group_order %in% anno[[group_id]]]
    
    group_order_df <- data.frame(group = group_order) %>%
      mutate(xpos = 1:n())
    names(group_order_df)[1] <- group_id
    
    plot_data <- plot_data %>%
      filter_(paste0(group_id, " %in% group_order")) %>%
      left_join(group_order_df, by = group_id) %>%
      left_join(group_n, by = group_id)
    
  } else {
    # Otherwise, arrange using the group_id for the group_by parameter, and use that order.
    group_order_df <- plot_data %>%
      select(one_of(group_id)) %>%
      arrange_(group_id) %>%
      mutate(xpos = 1:n())
    
    plot_data <- plot_data %>%
      left_join(group_order_df, by = group_id) %>%
      left_join(group_n, by = group_id)
  }
  
  n_genes <- length(genes)
  n_groups <- length(unique(plot_data[[group_id]]))
  n_samples <- nrow(data)
  
  header_labels <- build_header_labels(data = plot_data, 
                                      grouping = grouping,
                                      group_order = group_order,
                                      ymin = n_genes + 1, 
                                      label_height = label_height, 
                                      label_type = "simple")
  
  # Build the maximum value labels for the right edge
  max_labels <- data.frame(x = (n_groups + 0.5) * 1.01,
                           y = 1:n_genes + 0.5,
                           label = sci_label(max_vals) )
  max_header <- data.frame(x = (n_groups + 0.5) * 1.01,
                           y = n_genes + 1,
                           label = "Max value")
  max_width <- n_groups*(max_width/100)/(1-max_width/100)
  
  label_y_size <- max(header_labels$ymax) - min(header_labels$ymin)
  
  group_data <- plot_data %>%
    select(xpos, group_n) %>%
    mutate(label_y = n_genes + label_y_size * 0.05,
           group_n_y = max(header_labels$ymax) - 0.1 * label_y_size)
  
  # Plot setup
  p <- ggplot() +
    scale_fill_identity() +
    scale_y_continuous("", 
                       breaks = 1:length(genes) + 0.45, 
                       labels = genes, 
                       expand = c(0, 0)) +
    scale_x_continuous("", 
                       expand = c(0, 0)) +
    theme_classic(font_size) +
    theme(axis.text = element_text(size = rel(1)),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none") +
    geom_hline(aes(yintercept = 1:(n_genes)), size = 0.2)
  
  # plot the heatmap for each gene
  for (i in 1:length(genes)) {
    gene <- genes[[i]]
    
    p <- p + 
      geom_rect(data = plot_data,
                aes_string(xmin = "xpos - 0.5", xmax = "xpos + 0.5",
                           ymin = i, ymax = i + 1, 
                           fill = gene))
    
    
  }
  
  # Cluster labels
  p <- p +
    geom_rect(data = header_labels, 
              aes(xmin = xmin, ymin = ymin, 
                  xmax = xmax, ymax = ymax, 
                  fill = color)) +
    geom_text(data = header_labels,
              aes(x = (xmin + xmax) / 2, 
                  y = ymin + 0.05, 
                  label = label),
              angle = 90, 
              vjust = 0.35, hjust = 0, 
              size = pt2mm(font_size)) +
    # Maximum value labels on right side of plot
    geom_rect(aes(xmin = n_groups + 0.5, xmax = n_groups + 0.5 + max_width, 
                  ymin = 1, ymax = max(header_labels$ymax)),
              fill = "white") +
    geom_text(data = max_header,
              aes(x = x, y = y, 
                  label = label),
              angle = 90, 
              hjust = 0, vjust = 0.35, 
              size = pt2mm(font_size)) +
    geom_text(data = max_labels,
              aes(x = x, y = y, 
                  label = label),
              hjust = 0, vjust = 0.35, 
              size = pt2mm(font_size), 
              parse = TRUE)
  
  # Cluster counts
  if (show_counts) {
    if (rotate_counts) {
      p <- p + geom_text(data = group_data,
                         aes(x = xpos,
                             y = group_n_y, 
                             label = group_n),
                         angle = 90,
                         hjust = 1, vjust = 0.35, 
                         size = pt2mm(font_size))
    } else {
      p <- p + geom_text(data = group_data,
                         aes(x = xpos,
                             y = group_n_y, 
                             label = group_n),
                         size = pt2mm(font_size))
    }
  }
  
  return(p)
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
#'   \item "prop_gt0" (proportion of samples > 0)
#'   \item "prop_gt1" (proportion of samples > 1)
#'   \item "min"
#'   \item "max"
#'   }
#' @param size_stat The statistic to apply to each group for scaling dot size. Same options as fill_stat. Default = "prop_gt0".
#' @param log_scale Logical , determines if data is log scaled before plotting. Default = FALSE.
#' @param normalize_rows Logical, whether or not to rescale data within each row of the plot. Default = FALSE.
#' @param font_size numeric object, the font size (in pts) used to make the plot.
#' @param label_height numeric object, Percent of the plot height that should be used for the labels (0 to 100). Default is 25.
#' @param show_counts Logical, whether or not to display sample counts at the top of labels. Default = TRUE.
#' @param rotate_counts Logical, whether or not to rotate sample counts by 90 degrees. Default = FALSE.
#' @param max_width numeric object, percent of plot width that should be used for maximum expression values (0 to 100). Default is 10.
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
                           log_scale = TRUE,
                           normalize_rows = FALSE,
                           colorset = NULL,
                           font_size = 7, 
                           label_height = 25,
                           show_counts = TRUE, 
                           rotate_counts = FALSE,
                           max_width = 10) {
  library(dplyr)
  library(ggplot2)
  
  genes <- rev(genes)
  
  group_id <- paste0(grouping, "_id")
  group_label <- paste0(grouping, "_label")
  group_color <- paste0(grouping, "_color")
  
  gene_data <- data[,c("sample_name",genes)]
  gene_data <- gene_data[match(anno$sample_name, data$sample_name),]
  
  # Get maximum values for each gene before rescaling to plot space.
  max_vals <- map_dbl(genes, function(x) { max(gene_data[[x]]) })
  names(max_vals) <- genes
  
  gene_fill_stats <- group_stats(gene_data,
                                 value_cols = genes,
                                 anno = anno,
                                 grouping = group_label,
                                 stat = fill_stat)
  
  gene_size_stats <- group_stats(gene_data,
                                 value_cols = genes,
                                 anno = anno,
                                 grouping = group_label,
                                 stat = size_stat)
  
  if(log_scale) {
    gene_fill_stats[,genes] <- log10(gene_fill_stats[,genes] + 1)
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
    select(one_of(group_id, group_label, group_color)) %>%
    unique()
  
  group_n <- anno %>%
    group_by_(group_id) %>%
    summarise(group_n = n())
  
  plot_data <- plot_anno %>%
    left_join(gene_fill_data, by = group_label) %>%
    left_join(gene_size_stats, by = group_label)
  
  # Add an x position to each group
  if(!is.null(group_order)) {
    # Because we allow ranges, and groups may not necessarily be continuous integer sets
    # We have to filter out any that don't match first.
    group_order <- group_order[group_order %in% anno[[group_id]]]
    
    group_order_df <- data.frame(group = group_order) %>%
      mutate(xpos = 1:n())
    names(group_order_df)[1] <- group_id
    
    plot_data <- plot_data %>%
      filter_(paste0(group_id, " %in% group_order")) %>%
      left_join(group_order_df, by = group_id) %>%
      left_join(group_n, by = group_id)
    
  } else {
    # Otherwise, arrange using the group_id for the group_by parameter, and use that order.
    group_order_df <- plot_data %>%
      select(one_of(group_id)) %>%
      arrange_(group_id) %>%
      mutate(xpos = 1:n())
    
    plot_data <- plot_data %>%
      left_join(group_order_df, by = group_id) %>%
      left_join(group_n, by = group_id)
  }
  
  n_genes <- length(genes)
  n_groups <- length(unique(plot_data[[group_id]]))
  n_samples <- nrow(data)
  
  header_labels <-build_header_labels(data = plot_data, 
                                      grouping = grouping,
                                      group_order = group_order,
                                      ymin = n_genes + 1, 
                                      label_height = label_height, 
                                      label_type = "simple")
  
  # Build the maximum value labels for the right edge
  max_labels <- data.frame(x = (n_groups + 0.5) * 1.01,
                           y = 1:n_genes + 0.5,
                           label = sci_label(max_vals) )
  max_header <- data.frame(x = (n_groups + 0.5) * 1.01,
                           y = n_genes + 1,
                           label = "Max value")
  max_width <- n_groups*(max_width/100)/(1-max_width/100)
  
  label_y_size <- max(header_labels$ymax) - min(header_labels$ymin)
  
  group_data <- plot_data %>%
    select(xpos, group_n) %>%
    mutate(label_y = n_genes + label_y_size * 0.05,
           group_n_y = max(header_labels$ymax) - 0.1 * label_y_size)
  
  # Plot setup
  p <- ggplot() +
    scale_fill_identity() +
    scale_size_area() +
    scale_y_continuous("", 
                       breaks = 1:length(genes) + 0.45, 
                       labels = genes, 
                       expand = c(0, 0)) +
    scale_x_continuous("", 
                       expand = c(0, 0)) +
    theme_classic(font_size) +
    theme(axis.text = element_text(size = rel(1)),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none") +
    geom_hline(aes(yintercept = 1:(n_genes)), size = 0.2)
  
  # plot the swarms for each gene
  for(i in 1:length(genes)) {
    gene <- genes[[i]]
    gene_fill <- paste0(gene, "_fill")
    gene_size <- paste0(gene, "_size")
    p <- p + 
      geom_point(data = plot_data,
                 aes_string(x = "xpos",
                            y = i + 0.5, 
                            fill = gene_fill,
                            size = gene_size),
                 pch = 21)
    
    
  }
  
  # Cluster labels
  p <- p +
    geom_rect(data = header_labels, 
              aes(xmin = xmin, ymin = ymin, 
                  xmax = xmax, ymax = ymax, 
                  fill = color)) +
    geom_text(data = header_labels,
              aes(x = (xmin + xmax) / 2, 
                  y = ymin + 0.05, 
                  label = label),
              angle = 90, 
              vjust = 0.35, hjust = 0, 
              size = pt2mm(font_size)) +
    # Maximum value labels on right side of plot
    geom_rect(aes(xmin = n_groups + 0.5, xmax = n_groups + 0.5 + max_width, 
                  ymin = 1, ymax = max(header_labels$ymax)),
              fill = "white") +
    geom_text(data = max_header,
              aes(x = x, y = y, 
                  label = label),
              angle = 90, 
              hjust = 0, vjust = 0.35, 
              size = pt2mm(font_size)) +
    geom_text(data = max_labels,
              aes(x = x, y = y, 
                  label = label),
              hjust = 0, vjust = 0.35, 
              size = pt2mm(font_size), 
              parse = TRUE)
  
  # Cluster counts
  if(show_counts) {
    if(rotate_counts) {
      p <- p + geom_text(data = group_data,
                         aes(x = xpos,
                             y = group_n_y, 
                             label = group_n),
                         angle = 90,
                         hjust = 1, vjust = 0.35, 
                         size = pt2mm(font_size))
    } else {
      p <- p + geom_text(data = group_data,
                         aes(x = xpos,
                             y = group_n_y, 
                             label = group_n),
                         size = pt2mm(font_size))
    }
  }
  
  return(p)
}

#' Build a heatmap legend plot
#' 
#' @param minval numeric, the minimum value in the plot scale (default = 0)
#' @param maxval numeric, the maximum value in the plot scale (default = 4)
#' @param scale_name character, the name for the values displayed (default = "FPKM")
#' @param colorset character vector, the colors to interpolate between using colorRampPalette.
#' 
#' @return a ggplot2 heatmap legend plot
heatmap_legend_plot <- function(minval = 0, maxval = 4, scale_name = "FPKM", colorset=c("darkblue","dodgerblue","gray80","orange","orangered")) {
  
  library(ggplot2)
  
  colors <- colorRampPalette(colorset)(1001)
  
  ## Build geom_rect() compatible table
  legend_data <- data.frame(xmin = 1:1001,
                            xmax = 1:1001+1,
                            ymin = 0,
                            ymax = 1,
                            fill = colors)
  
  legend_plot <- ggplot(legend_data) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill)) +
    geom_segment(aes(x = min(xmin), xend = max(xmax), y = 0, yend = 0)) +
    scale_fill_identity() +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(scale_name, breaks=c(0,250,500,750,1000),
                       labels=round(seq(minval, maxval, by = (maxval-minval)/4),2)) +
    theme_classic() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          axis.title.y = element_blank(),
          axis.line.x = element_blank())
  
  return(legend_plot)
  
}