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
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) %+replace%
    ggplot2::theme(plot.margin = unit(c(rep(0,4)),"line"),
                   axis.text = element_text(size = rel(1)),
                   axis.text.x = element_blank(),
                   axis.ticks.x = element_blank(),
                   axis.title.x = element_blank(),
                   axis.ticks.margin = unit(0,"cm"),
                   axis.ticks.x = element_blank())
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
  
  library(dplyr)
  
  group_id <- paste0(grouping, "_id")
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
        left_join(group_order_df, by = group_id) %>%
        arrange(group_order) %>%
        mutate(xpos = 1:n())
      
    } else {
      # Otherwise, arrange using the group_id for the group_by parameter, and use that order.
      data <- data %>%
        dplyr::arrange({{ group_id }}) %>%
        dplyr::mutate(xpos = 1:n())
    }
  }
  
  poly.data <- data %>% 
    dplyr::group_by(!!rlang::parse_expr(group_id)) %>%
    dplyr::summarise(color = .data[[group_color]][1],
                     x1 = max(xpos),
                     x2 = min(xpos) - 1) %>%
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
  group_label <- paste0(grouping, "_label")
  group_color <- paste0(grouping, "_color")
  
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
        dplyr::arrange(!!rlang::parse_expr(group_id)) %>%
        dplyr::mutate(xpos = 1:n())
    }
  }
  
  grouping <- lapply(c(group_id, group_label, group_color),
                     rlang::parse_expr)
  
  data <- data %>%
    dplyr::group_by(!!!grouping) %>%
    dplyr::summarise(minx = min(xpos),
                     maxx = max(xpos))
  
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
      dplyr::group_by(!!parse_expr(group_id)) %>%
      summarise(xmin = minx - 1,
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
  hc.segs <- dendextend::as.ggdend(hc.dendro)$segments
  
  ymin = min(dir.lims)
  ymax = max(dir.lims)
  
  yheight = ymax - ymin
  
  norm.segs <- hc.segs %>%
    mutate(y = (y/max(y))*yheight + ymin) %>%
    mutate(yend = (yend/max(yend))*yheight + ymin)
  
  if (tree.dir == "down") {
    
    plot.segs <- norm.segs
    
  } else if (tree.dir == "up") {
    
    ycenter = (ymin + ymax) / 2
    
    plot.segs <- norm.segs %>%
      mutate(y = ycenter + (ycenter - y),
             yend = ycenter + (ycenter - yend))
    
  } else if (tree.dir == "left") {
    
    plot.segs <- norm.segs
    names(plot.segs) <- c("y","x","yend","xend")
    
  } else if (tree.dir == "right") {
    xcenter = (ymin + ymax) / 2
    
    plot.segs <- norm.segs
    names(plot.segs) <- c("y","x","yend","xend")
    
    plot.segs <- plot.segs %>%
      mutate(x = xcenter + (xcenter - x),
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
spiral_jitter <- function(x, 
                          y, 
                          n = NULL, 
                          max_n = NULL, 
                          radius = 1, 
                          aspect = 1, 
                          ratio = "golden") {
  
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
