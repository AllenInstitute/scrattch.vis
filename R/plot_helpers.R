get_n_stats <- function(plot_data, group_cols, genes) {
  list(genes = length(genes),
       groups = length(unique(plot_data[[group_cols$id]])),
       samples = nrow(plot_data))
}

group_columns <- function(grouping) {
  list(id = paste0(grouping, "_id"),
       label = paste0(grouping, "_label"),
       color = paste0(grouping, "_color"))
}

filter_gene_data <- function(data, genes, anno, group_cols, group_order = NULL, sample_key_column = "sample_name") {
  
  if(!is.null(group_order)) {
    group_order <- group_order[group_order %in% anno[[group_cols$id]]]
    
    anno_filter <- anno[[group_cols$id]] %in% group_order
    anno <- anno[anno_filter, ]
  }
  
  row_filter <- match(anno[[sample_key_column]], data[[sample_key_column]])
  col_filter <- c(sample_key_column, genes)
  
  data[row_filter, col_filter]
  
}

max_gene_vals <- function(data, genes) {
  purrr::map_dbl(genes, function(x) { max(data[[x]], na.rm = TRUE) })
}

scale_gene_data <- function(data, genes, scale_type = "log10") {
  if(scale_type == "linear") {
    data
  } else if(scale_type == "log10") {
    data[, genes] <- log10(data[, genes] + 1)
    data
  }
}


add_sample_xpos <- function(data, group_cols, group_order = NULL) {
  if(is.null(group_order)) {
    data %>%
      dplyr::arrange_(group_cols$id) %>%
      dplyr::mutate(xpos = 1:n())
  } else {
    group_order_df <- data.frame(group = group_order) %>%
      dplyr::mutate(.plot_order = 1:n())
    names(group_order_df)[1] <- group_id
    
    group_filter <- data[[group_cols$id]] %in% group_order
    
    data[group_filter, ] %>%
      dplyr::left_join(group_order_df, 
                       by = group_cols$id) %>%
      dplyr::arrange(.plot_order) %>%
      dplyr::mutate(xpos = 1:n()) %>%
      dplyr::select(-.plot_order)
  }
}

build_max_dfs <- function(n_stats, max_vals, max_width) {
  labels <- data.frame(x = n_stats$samples * 1.01,
                       y = 1:n_stats$genes + 0.5,
                       label = sci_label(max_vals))
  header <- data.frame(x = n_stats$samples * 1.01,
                       y = n_stats$genes + 1,
                       label = "Max value")
  width <- n_stats$samples * (max_width / 100) / (1 - max_width / 100)
  
  return(list(labels = labels,
              header = header,
              width = width))
}

scale_values_plot_space <- function(x, min_ps, max_ps = NULL, min_val = 0, max_val = NULL, extent = 0.9) {
  if(is.null(max_ps)) { max_ps <- min_ps + 1 }
  if(min_val > 0) { x <- x - min_val; x[x < 0] <- 0 }
  if(is.null(max_val)) { max_val = max(x, na.rm = TRUE) }
  x / max_val * extent * (max_ps - min_ps) + min_ps
}

ggplot_header_labels <- function(p, header_labels, header_polygons = NULL, font_size) {
  p <- p + 
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
              size = pt2mm(font_size))
  
  
  if(!is.null(header_polygons)) {
    p <- p +
      geom_polygon(data = header_polygons,
                   aes(x = poly.x, 
                       y = poly.y, 
                       fill = color, 
                       group = id) )
  }
  
  p
}


ggplot_scale_bars <- function(p, n_genes, n_samples, extent = 0.9) {
  # Calculate Plot Scale bars
  scale_bars <- data.frame(ymin = 1:n_genes,
                           ymid = 1:n_genes + extent / 2,
                           ymax = 1:n_genes + extent,
                           xmin = -n_samples * 0.01,
                           xmax = 0)
  
  p +
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
                 size = 0.2)
}

ggplot_max_vals <- function(p, n_stats, max_val_dfs, font_size) {
  p +
    geom_rect(aes(xmin = n_stats$samples + 1, 
                  xmax = n_stats$samples + max_val_dfs$width, 
                  ymin = 1, 
                  ymax = 1), 
              fill = "#FFFFFF") +
    geom_text(data = max_val_dfs$header,
              aes(x = x, 
                  y = y, 
                  label = label),
              angle = 90, 
              hjust = 0, 
              vjust = 0.5, 
              size = pt2mm(font_size) ) +
    geom_text(data = max_val_dfs$labels,
              aes(x = x, 
                  y = y, 
                  label = label),
              hjust = 0, 
              vjust = 0.5, 
              size = pt2mm(font_size) , 
              parse = TRUE)
}