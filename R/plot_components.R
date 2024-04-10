#' Convert integers to scientific notation labels
#' 
#' @param in_num a numeric vector
#' @param sig_figs a number indicating how many significant figures should be displayed.
#' @param type Either "ggplot2", "DT", or "text", which will set how the values are returned. 
#' "ggplot2" gives formatting good for ggplot2 labels, e.g. 4.7\%*\%10^-4. "DT" returns 
#' formatting good for use with the DT package, e.g. 4.7\\u271510<sup>-4</sup>. "text" returns 
#' basic text formatting, e.g. 4.7E-4
#' 
#' @return a character vector with numeric values reformatted for display
#' 
#' @examples
#' my_numbers <- c(100,15.359,32687,.000468)
#' 
#' sci_label(my_numbers)
#' 
#' sci_label(my_numbers, sig_figs = 3)

sci_label <- function(in_num, 
                      sig_figs = 2, 
                      type = "ggplot2") {
  
  labels <- character()
  
  for(i in 1:length(in_num)) {
    x <- in_num[i]
    
    if(x < 0) {

      neg <- "-"
      x <- abs(x)
    } else {
      neg <- ""
    }
    
    # Format the string to adjust for number of sig_figs
    if (x == 0) {
      # If the value is 0, build 0.(0)N based on the number of sig figs requested.
      first <- paste0("0", ".", paste0(rep("0", sig_figs - 1), collapse = "") )
    } else if (log10(x) %% 1 == 0) {
      first <- substr(x, 1, 1)
      if (sig_figs > 1) {
        first <- paste0(first, ".", paste0(rep("0", sig_figs - 1), collapse = ""))
      }
    } else {
      first <- round(x / (10 ^ floor(log10(x))), sig_figs - 1)
      if (nchar(first) < sig_figs + 1) {
        if (first %% 1 == 0) {
          first <- paste0(first, ".", paste0(rep("0", sig_figs - 1), collapse = ""))
        } else {
          # +1 because of decimal place
          first <- paste0(first, paste0(rep("0",sig_figs + 1 - nchar(first)), collapse = ""))
        }
        
      }
    }
    # Add suffixes based on type parameter
    if (x == 0) {
      if (type == "text") {
        label <- paste0(first,"E0")
      } else if (type == "ggplot2") {
        label <- paste0(first, "%*%10^0" )
      } else if (type == "DT") {
        label <- paste0(first, "\u271510<sup>0</sup>" )
      }
    } else {
      if (type == "text") {
        label <- paste0(first,"E",floor(log10(x)))
      } else if (type == "ggplot2") {
        label <- paste0(first, "%*%10^", floor(log10(x)))
      } else if (type == "DT") {
        label <- paste0(first, "\u271510<sup>", floor(log10(x)),"</sup>")
      }
    }
    label <- paste0(neg, label)
    labels <- c(labels, label)
  }
  return(labels)
}

#' Remove the X-axis (and most other margins)
#' 
#' Makes plots more suitable for use with Illustrator by removing most margins
#' and the X-axis (which is usually replaced by something else in my plots).
#' 
#' To further remove the space below the x-axis, use labs(x = NULL)
#' 
#' Based on theme_classic() from ggplot2.
#' 
#' @examples
#' ggplot(mtcars) +
#'  geom_point(aes(x = mpg, y = wt)) +
#'  theme_no_x() +
#'  labs(x = NULL)
theme_no_x <- function(base_size = 12, base_family = "") {
  ggplot2::`%+replace%`(
    ggplot2::theme_classic(base_size = base_size, base_family = base_family),
    ggplot2::theme(plot.margin = unit(c(rep(0,4)),"line"),
                   axis.text = ggplot2::element_text(size = rel(1)),
                   axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_blank(),
                   axis.ticks.margin = ggplot2::unit(0,"cm"))
  )
}

#' Build polygons from plot data for fancy headers built into the plot area
#' 
#' @param data A data.frame containing joined annotations and expression data
#' @param grouping The base to use for grouping samples.
#' @param ymin The minimum y value for the bottom of the headers. In gene plot contexts, this will usually be the number of genes + 1
#' @param label_height Percentage of the plot area that the headers should take up. Default = 25. 
#' @param poly_type Either "angle" or "square". Angle will scale the header bars to be equal width, with a layer of 
#' angled polygons connecting them to the data. Square will construct headers with rectangular bars that match the width 
#' of the samples in each group.
#' 
#' 
build_header_polygons <- function(data,
                                  anno,
                                  grouping,
                                  group_order = NULL,
                                  ymin, 
                                  label_height = 25, 
                                  fraction_of_label = 0.1,
                                  poly_type = "angle") {
  # Two label types: 
  # "angle" will draw a polygon with the base lined up with samples, and the top
  # divided evenly for each cluster.
  # "square" will draw a rectangle with the top points matching bottom points
  
  group_id <- paste0(grouping, "_id")
  parsed_group_id <- rlang::parse_expr(group_id)
  group_color <- paste0(grouping, "_color")
  n_groups <- length(unique(anno[[group_id]]))
  n_samples <- nrow(anno)
  
  ## Note on plot dimensions
  # The range of the plot area (not including labels) will be
  # y-axis: 1:ymin (ngenes + 1)
  # x-axis: 0:(n_samples)
  
  # Calculate the height of the label in plot dimensions:
  labheight <- (ymin - 1) * (label_height / 100) / (1 - label_height / 100)
  
  # Build cell type label polygons
  # polygon points are built in this order: 1 = bottom-right, 2 = bottom-left, 3 = top-left, 4 = top-right
  # For angled labels, the bottom two x positions are calculated based on the number of samples
  # in each cluster. The top positions are evenly spaced based on the number of clusters.
  # Add an x position to each group
  
  data <- dplyr::left_join(anno, data, by = "sample_name")

  if(!"xpos" %in% names(data)) {
    if(!is.null(group_order)) {
      group_order_df <- build_vec_pos(group_order,
                                      vec_name = group_id,
                                      sort = "none",
                                      axis_name = "group_order")
      
      data <- data %>%
        dplyr::left_join(group_order_df, by = group_id) %>%
        dplyr::arrange(group_order) %>%
        dplyr::mutate(xpos = 1:n())
      
    } else {
      # Otherwise, arrange using the group_id for the group_by parameter, and use that order.
      
      data <- data %>%
        dplyr::arrange(!!parsed_group_id) %>%
        dplyr::mutate(xpos = 1:n())
    }
  }
  
  poly.data <- data %>% 
    dplyr::group_by(!!parsed_group_id) %>%
    dplyr::summarise(color = .data[[group_color]][1],
                     x1 = max(xpos),
                     x2 = min(xpos) - 1) %>%
    dplyr::arrange(x1) %>%
    dplyr::mutate(x3 = (n_samples) * (1:n_groups - 1) / n_groups,
                  x4 = (n_samples) * (1:n_groups) / n_groups,
                  # ymin is the top of the plot body
                  y1 = ymin,
                  y2 = ymin,
                  # The angled portion of the label will be 10% of the total label height 
                  y3 = ymin + labheight * fraction_of_label,
                  y4 = ymin + labheight * fraction_of_label)
  
  # For a simpler square label, set the top and bottom x-positions to be the same
  if (poly_type == "square") {
    poly.data <- poly.data %>%
      dplyr::mutate(x3 = x2,
                    x4 = x1)
  }
  
  # Restructure the polygons for ggplot2's geom_poly().
  # The data should have a single x and y column in order, with id and color for each polygon
  poly <- data.frame(id = rep(poly.data[[group_id]], each = 4),
                     color = rep(poly.data$color, each = 4))
  poly.x <- numeric()
  poly.y <- numeric()
  for (i in 1:nrow(poly.data)) {
    poly.x <- c(poly.x,
                poly.data$x1[i],
                poly.data$x2[i],
                poly.data$x3[i],
                poly.data$x4[i])
    poly.y <- c(poly.y,
                poly.data$y1[i],
                poly.data$y2[i],
                poly.data$y3[i],
                poly.data$y4[i])
  }
  poly <- cbind(poly,
                poly.x = poly.x,
                poly.y = poly.y)
  
  poly
}

#' Build colorful, rectangular labels for plot headers in plot space
#' 
#' @param data A data.frame containing joined annotations and expression data
#' @param grouping The base to use for grouping samples.
#' @param group_order Alternate order to use for grouping if necessary. Default is NULL, which will use grouping_id.
#' @param ymin The minimum y value for the bottom of the headers. In gene plot contexts, this will usually be the number of genes + 1
#' @param label_height Percentage of the plot area that the headers should take up. Default = 25. 
#' @param label_type Either "simple", "angle", or "square". Simple is for use with grouped plots. 
#' Angle will scale the header bars to be equal width, with a layer of 
#' angled polygons connecting them to the data. Square will construct headers with rectangular bars that match the width 
#' of the samples in each group.
#' 
build_header_labels <- function(data, 
                                grouping, 
                                group_order = NULL,
                                ymin,  
                                label_height = 25, 
                                label_type = "simple") {
  
  # Three label types: 
  # simple, which is for use with cluster-based plots
  # angle, for cell-based plots with "angle"-type polygonal labels
  # square, for cell-based plots with "square"-type labels
  
  ## Note on plot dimensions
  # The range of the plot area (not including labels) will be
  # y-axis: 1:(ngenes + 1)
  # x-axis: 0:(nsamples) (for cell-based plots)
  # x-axis: 1:(nclust + 1) (for cluster-based plots)
  
  labheight <- (ymin - 1)*(label_height/100)/(1 - label_height/100)
  
  group_id <- paste0(grouping, "_id")
  parsed_group_id <- rlang::parse_expr(group_id)
  group_label <- paste0(grouping, "_label")
  parsed_group_label <- rlang::parse_expr(group_label)
  group_color <- paste0(grouping, "_color")
  parsed_group_color <- rlang::parse_expr(group_color)
  
  n_samples <- nrow(data)
  n_clust <- length(unique(data[[group_id]]))
  
  # Add an x position to each group
  if(!"xpos" %in% names(data)) {
    if(!is.null(group_order)) {
      group_order_df <- build_vec_pos(group_order,
                                      vec_name = group_id,
                                      sort = "none",
                                      axis_name = "xpos")
      
      data <- data %>%
        dplyr::left_join(group_order_df, by = group_id)
      
    } else {
      # Otherwise, arrange using the group_id for the group_by parameter, and use that order.
      data <- data %>%
        dplyr::arrange(!!parsed_group_id) %>%
        dplyr::mutate(xpos = 1:n())
    }
  }
  
  data <- data %>%
    dplyr::group_by(!!parsed_group_id, !!parsed_group_label, !!parsed_group_color) %>%
    dplyr::summarise(minx = min(xpos),
                     maxx = max(xpos)) %>%
    dplyr::arrange(minx)
  
  if (label_type == "simple") {
    xlab.rect <- data.frame(xmin = data$minx - 0.5,
                            xmax = data$maxx + 0.5,
                            ymin = ymin,
                            ymax = ymin + labheight,
                            color = data[[group_color]],
                            label = data[[group_label]] )
  } else if (label_type == "angle") {
    xlab.rect <- data.frame(xmin = (n_samples) * (1:n_clust - 1) / n_clust,
                            xmax = (n_samples) * (1:n_clust) / n_clust,
                            # 10% of the label height is reserved for angled polygons
                            ymin = ymin + labheight * 0.1,
                            ymax = ymin + labheight,
                            color = data[[group_color]],
                            label = data[[group_label]] )
  } else if (label_type == "square") {
    xlab.rect <- data %>% 
      dplyr::group_by(!!parsed_group_id) %>%
      dplyr::summarise(xmin = minx - 1,
                       xmax = maxx,
                       ymin = ymin + labheight * 0.1,
                       ymax = ymin + labheight,
                       color = .[[group_color]][1],
                       label = .[[group_label]][1])
  }

  xlab.rect  
}

#' Covert hclust objects to segments for use in ggplots
#' 
#' @param hc a hclust object (hierarchical clustering, stats package). First convert df to distances (dist(df)) then hclust(df)
#' @param tree.dir a character object with the direction the tree points to, from root to leaves. options are "down" (default), "up","left", "right".
#' @param dir.lims a 2-member vector with the space in the direction of plotting that the dendrogram will occupy. default = c(0,1)
#' @return a data.frame with segment values for ggplot2's geom_seg. columns: "x","xend","y","yend".
hclust_to_seg <- function(hc, tree.dir = "down", dir.lims = c(0,1)) {
  
  hc.dendro <- as.dendrogram(hc)
  hc.segs <- as.data.frame(ggdendro::segment(ggdendro::dendro_data(hc.dendro)))
  
  ymin = min(dir.lims)
  ymax = max(dir.lims)
  
  yheight = ymax - ymin
  
  norm.segs <- hc.segs %>%
    dplyr::mutate(y = (y/max(y))*yheight + ymin) %>%
    dplyr::mutate(yend = (yend/max(yend))*yheight + ymin)
  
  if (tree.dir == "down") {
    
    plot.segs <- norm.segs
    
  } else if (tree.dir == "up") {
    
    ycenter = (ymin + ymax) / 2
    
    plot.segs <- norm.segs %>%
      dplyr::mutate(y = ycenter + (ycenter - y),
                    yend = ycenter + (ycenter - yend))
    
  } else if (tree.dir == "left") {
    
    plot.segs <- norm.segs
    names(plot.segs) <- c("y","x","yend","xend")
    
  } else if (tree.dir == "right") {
    xcenter = (ymin + ymax) / 2
    
    plot.segs <- norm.segs
    names(plot.segs) <- c("y","x","yend","xend")
    
    plot.segs <- plot.segs %>%
      dplyr::mutate(x = xcenter + (xcenter - x),
                    xend = xcenter + (xcenter - xend))
    
  }
  
  return(plot.segs)
  
}

#' Jitter x-y coordinates in a spiral pattern
#'
#' @param x coordinate x
#' @param y coordinate y
#' @param n length(x)
#' @param max_n 
#' @param radius 
#' @param aspect 
#' @param ratio golden
#' @examples 
#' spiral_jitter(x,y)
spiral_jitter <- function(x, y, n = NULL, max_n = NULL, radius = 1, aspect = 1, ratio = "golden") {
  
  if (is.null(n)) {
    n <- length(x)
  }
  
  # pre-calculated golden ratio
  if (ratio == "golden") {
    ratio <- (sqrt(5) + 1) / 2
  }
  
  # calculate the angle based on the ratio
  angle <- 360 / (ratio ^ 2)
  
  # scale the spacing between points
  if (is.null(max_n)) {
    # If no maximum is provided, then c is based
    # on the number of points (n)
    c <- radius / sqrt(n)
  } else {
    # If a maximum is provided, then c is based
    # on the maximum number of points
    c <- radius / sqrt(max_n)
  }
  
  # set up vectors for the jittered results
  x_j <- rep(x,n)
  y_j <- rep(y,n)
  
  # Jitter each point
  for (m in 1:n) {
    # calculate position using polar coordinates
    r <- c * sqrt(m)
    theta <- angle * m
    
    # convert polar coordinates to cartesian coordinates
    if (aspect > 1) {
      # if the aspect is larger than 1 (wider than tall),
      # scale the y positions to compensate
      x_j[m] <- x + r * cos(theta)
      y_j[m] <- y + r * sin(theta) * 1/aspect
    } else {
      # if the aspect is smaller than 1 (taller than wide),
      # scale the x positions to compensate
      x_j[m] <- x + r * cos(theta) * aspect
      y_j[m] <- y + r * sin(theta) 
      
    }
    
  }
  
  results <- data.frame(x = x_j, y = y_j)
  
  return(results)
  
}

#' Compute basic count statistics for use in generating plots
#' 
#' @param plot_data a data.frame with joined annotations and gene expression data
#' @param group_cols a group_cols list generated by group_columns()
#' @param genes a vector of genes used for the plot
#' 
#' @return a list of 3 count statistics: genes, groups, and samples
get_n_stats <- function(plot_data, group_cols, genes) {
  list(genes = length(genes),
       groups = length(unique(plot_data[[group_cols$id]])),
       samples = nrow(plot_data))
}

#' Generate a list of anno columns based on grouping
#' 
#' @param grouping the base of the grouping to use.
#' 
#' @return a list with 3 character objects: id, label, and color
#' 
group_columns <- function(grouping) {
  list(id = paste0(grouping, "_id"),
       label = paste0(grouping, "_label"),
       color = paste0(grouping, "_color"))
}

#' Filter gene data based on annotations, genes, and groups supplied with group_order.
#' 
#' @param data a data.frame with gene expression data and a sample_name
#' @param genes a vector of genes to plot
#' @param anno a data.frame of sample annotations
#' @param group_cols a group_cols list generated by group_columns()
#' @param group_order a set of group_ids to use for the plot. Default is NULL, which ignores filtering on this parameter.
#' @param sample_key_column the name of the key column for sample identification. Default is "sample_name"
#' 
#' @return a data.frame of data filtered based on genes, the sample_names in anno, and (if supplied) the groups in group_order.
filter_gene_data <- function(data, 
                             genes, 
                             anno, 
                             group_cols, 
                             group_order = NULL, 
                             sample_key_column = "sample_name") {
  
  if(!is.null(group_order)) {
    group_order <- group_order[group_order %in% anno[[group_cols$id]]]
    
    anno_filter <- anno[[group_cols$id]] %in% group_order
    anno <- anno[anno_filter, ]
  }
  
  row_filter <- match(anno[[sample_key_column]], data[[sample_key_column]])
  col_filter <- c(sample_key_column, genes)
  
  data[row_filter, col_filter]
  
}

#' Compute maximum expression values in data for a set of genes.
#' 
#' @param data a data.frame with gene expression data
#' @param genes the genes to plot
#' 
#' @return a numeric vector with maximum expression values named for each gene.
max_gene_vals <- function(data, genes) {
  max_vals <- purrr::map_dbl(genes, function(x) { max(data[[x]], na.rm = TRUE) })
  names(max_vals) <- genes
  
  max_vals
}

#' Scale gene data columns
#' 
#' @param data a data.frame with gene expression data
#' @param genes the genes to plot
#' @param scale_type the type of scaling to supply. Currently only supports "log10" (default), and "linear" (no rescaling).
#' 
#' @return a data.frame with columns matching genes rescaled based on scale_type.
scale_gene_data <- function(data, genes, scale_type = "log10") {
  if(scale_type == "linear") {
    data
  } else if(scale_type == "log10") {
    data[, genes] <- log10(data[, genes] + 1)
    data
  }
}

#' Add a per-sample xpos column for sample data based on grouping
#' 
#' @param data a data.frame with both annotations and gene expression vales
#' @param group_cols a group_cols list generated by group_columns()
#' @param group_order an optional vector of group ids for ordering (default is NULL).
#' 
#' @return a data.frame of data with a column xpos appended.
#' 
add_sample_xpos <- function(data, group_cols, group_order = NULL) {
  if(is.null(group_order)) {
    parsed_id <- rlang::parse_expr(group_cols$id)
    data %>%
      dplyr::arrange(!!parsed_id) %>%
      dplyr::mutate(xpos = 1:n())
  } else {
    group_order_df <- data.frame(group = group_order) %>%
      dplyr::mutate(.plot_order = 1:n())
    names(group_order_df)[1] <- group_cols$id
    
    group_filter <- data[[group_cols$id]] %in% group_order
    
    data[group_filter, ] %>%
      dplyr::left_join(group_order_df, 
                       by = group_cols$id) %>%
      dplyr::arrange(.plot_order) %>%
      dplyr::mutate(xpos = 1:n()) %>%
      dplyr::select(-.plot_order)
  }
}

#' Adda per-group xpos column for sample data based on grouping
#' 
#' @param data a data.frame with both annotations and gene expression vales
#' @param group_cols a group_cols list generated by group_columns()
#' @param group_order an optional vector of group ids for ordering (default is NULL).
#' 
#' @return a data.frame of data with a column xpos appended.
#' 
add_group_xpos <- function(data, group_cols, group_order = NULL) {
  if(is.null(group_order)) {
    parsed_id <- rlang::parse_expr(group_cols$id)
    group_order_df <- data %>%
      dplyr::select(one_of(group_cols$id)) %>%
      unique() %>%
      dplyr::arrange(!!parsed_id) %>%
      dplyr::mutate(xpos = 1:n())
    
    data <- data %>%
      dplyr::left_join(group_order_df, by = group_cols$id)
  } else {
    group_order_df <- data.frame(group = group_order) %>%
      dplyr::mutate(xpos = 1:n())
    names(group_order_df)[1] <- group_cols$id
    
    data <- data %>%
      dplyr::left_join(group_order_df, by = group_cols$id)
  }
  
  data
}

#' Build data.frames for maximum value positions in plot space
#' 
#' @param n_stats count statistics made by get_n_stats()
#' @param width_stat The count used to determine x-position. In sample plots, use "samples". in group plots, use "groups".
#' @param max_vals maximum values per gene generated by max_gene_vals()
#' @param max_width the percentage of the plot that the max values should occupy.
#' 
#' @return a list with 3 data.frames: labels with label positions and values; header, with header position and label; width with the width in plot space. 
build_max_dfs <- function(n_stats, width_stat = "samples", max_vals, max_width) {
  xpos <- ifelse(width_stat == "samples",
                 n_stats[[width_stat]] * 1.01,
                 n_stats[[width_stat]] + 0.75)
  labels <- data.frame(x = xpos,
                       y = 1:n_stats$genes + 0.5,
                       label = sci_label(max_vals))
  header <- data.frame(x = xpos,
                       y = n_stats$genes + 1,
                       label = "Max value")
  width <- n_stats[[width_stat]] * (max_width / 100) / (1 - max_width / 100)
  
  return(list(labels = labels,
              header = header,
              width = width))
}

#' Scale values to plot space
#' 
#' @param x A vector of values to scale
#' @param min_ps The minimum value in plot space to use
#' @param max_ps The maximum value in plot space to use. Default is NULL, which will use min_ps + 1.
#' @param min_val The minimum value to plot. Default is 0. If NULL, will use min(x)
#' @param max_val The maximum value to plot. Default is NULL, which will use max(x)
#' @param extent The fraction of the space between min_ps and max_ps to use. Default is 0.9.
#' 
#' @return a vector of scaled values.
scale_values_plot_space <- function(x, 
                                    min_ps, 
                                    max_ps = NULL, 
                                    min_val = 0, 
                                    max_val = NULL, 
                                    extent = 0.9) {
  if(is.null(max_ps)) { max_ps <- min_ps + 1 }
  if(is.null(min_val)) { min_val <- min(x, na.rm = TRUE) }
  if(min_val > 0) { x <- x - min_val; x[x < 0] <- 0 }
  if(is.null(max_val)) { max_val = max(x, na.rm = TRUE) }
  x / max_val * extent * (max_ps - min_ps) + min_ps
}

#' Add header labels objects to a ggplot
#' 
#' @param p The plot to add headers to
#' @param header_labels A data.frame of header label coordinates generated by build_header_labels()
#' @param header_polygons An optional data.frame of header polygons generated by build_header_polygons(). Default is NULL.
#' @param font_size The font size in pt to use for header label plotting.
#' 
#' @return a ggplot2 object
#' 
ggplot_header_labels <- function(p, header_labels, header_polygons = NULL, font_size) {
  p <- p + 
    ggplot2::geom_rect(data = header_labels,
                       ggplot2::aes(xmin = xmin , 
                       xmax = xmax, 
                       ymin = ymin, 
                       ymax = ymax, 
                       fill = color) ) +
    ggplot2::geom_text(data = header_labels, 
                       ggplot2::aes(x = (xmin + xmax) / 2, 
                       y = ymin + 0.05, 
                       label = label),
              angle = 90, 
              vjust = 0.35, 
              hjust = 0, 
              size = pt2mm(font_size))
  
  
  if(!is.null(header_polygons)) {
    p <- p +
      ggplot2::geom_polygon(data = header_polygons,
                            ggplot2::aes(x = poly.x, 
                            y = poly.y, 
                            fill = color, 
                            group = id) )
  }
  
  p
}

#' Add scale bar objects to a ggplot
#' 
#' @param p The plot to add scale bars to
#' @param n_genes The number of genes in the plot
#' @param n_samples The number of samples in the plot
#' @param extent The amount of space in the plot row that each row of data uses.
#' 
#' @return a ggplot2 object
#' 
ggplot_scale_bars <- function(p, n_genes, n_samples, extent = 0.9) {
  # Calculate Plot Scale bars
  scale_bars <- data.frame(ymin = 1:n_genes,
                           ymid = 1:n_genes + extent / 2,
                           ymax = 1:n_genes + extent,
                           xmin = -n_samples * 0.01,
                           xmax = 0)
  
  p +
    ggplot2::geom_hline(data = scale_bars,
                        ggplot2::aes(yintercept = ymin), 
                                     linewidth = 0.2) +
    ggplot2::geom_segment(data = scale_bars, 
                          ggplot2::aes(x = xmin,
                                      xend = xmax,
                                      y = ymid, 
                                      yend = ymid),
                          linewidth = 0.2) +
    ggplot2::geom_segment(data = scale_bars,
                          ggplot2::aes(x = xmin, 
                                       xend = xmax, 
                                       y = ymax, 
                                       yend = ymax),
                          linewidth = 0.2) +
    ggplot2::geom_segment(data = scale_bars,
                          ggplot2::aes(x = xmax, 
                                       xend = xmax, 
                                       y = ymin, 
                                       yend = ymax),
                          linewidth = 0.2)
}

#' Add max value labels to a ggplot
#' 
#' @param p The plot to add max values to
#' @param n_stats Tcount statistics made by get_n_stats()
#' @param width_stat The count used to determine x-position. In sample plots, use "samples". in group plots, use "groups".
#' @param max_value_dfs The max value list generated by build_max_dfs()
#' @param font_size The font size in pt to use for header label plotting.
#' 
#' @return a ggplot2 object
#' 
ggplot_max_vals <- function(p, n_stats, width_stat = "samples", max_val_dfs, font_size) {
  
  right_pad <- data.frame(xmin = n_stats[[width_stat]] + 1,
                          xmax = n_stats[[width_stat]] + max_val_dfs$width)
  
  p +
    ggplot2::geom_rect(data = right_pad,
                       ggplot2::aes(xmin = xmin,
                                    xmax = xmax,
                                    ymin = 1,
                                    ymax = 1),
                       fill = "#FFFFFF") +
    ggplot2::geom_text(data = max_val_dfs$header,
                       ggplot2::aes(x = x, 
                                    y = y, 
                                    label = label),
                       angle = 90, 
                       hjust = 0, 
                       vjust = 0.5, 
                       size = pt2mm(font_size) ) +
    ggplot2::geom_text(data = max_val_dfs$labels,
                       ggplot2::aes(x = x, 
                                    y = y, 
                                    label = label),
                       hjust = 0, 
                       vjust = 0.5, 
                       size = pt2mm(font_size) , 
                       parse = TRUE)
}