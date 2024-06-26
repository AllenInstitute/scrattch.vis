# Error function
# used to generate sigmoidal curve 

#' Simple error function
#' 
#' @param x input for error function
#' @examples 
#' erf(0.5)
#' erf(0.1) 
erf <- function(x) {
  2 * pnorm(x * sqrt(2)) - 1
}

# Generates x-y coordinates for a sigmoidal line from x,y to xend,yend with
# the given number of steps.
# additional arguments define what curve to use
# and how far in the x-direction to use it

#' Build a sigmoidal line using sigfun
#' 
#' @param x starting x position. Default = 0
#' @param xend ending x position. Default = 1
#' @param y starting y position. Default = 0
#' @param yend ending y position. Default = 1
#' @param steps Number of points in the line. Default = 50
#' @param sigfun The sigmoidal function to use. Currently only supports "erf"
#' @param sigx Will generate values for -sigx to +sigx using sigfun, then rescale to x and y params.
#' 
#' @return a data.frame with x and y coordinates for the sigmoidal line
sigline <- function(x = 0, xend = 1, 
                    y = 0, yend = 1, 
                    steps = 50,
                    sigfun = "erf",
                    sigx = 1.5) {
  
  xsteps <- seq(-sigx, sigx, length.out = steps)
  
  if (sigfun == "erf") {
    ysteps <- erf(xsteps)
  }
  
  xsteps <- (xsteps + sigx)/(2*sigx)
  ysteps <- (ysteps + max(ysteps))/(2*max(ysteps))
  
  xscaled <- xsteps*(xend - x) + x
  yscaled <- ysteps*(yend - y) + y
  
  data.frame(x = c(x, xscaled, xend),
             y = c(y, yscaled, yend))
  
}

#' Expand a sigmoidal line to a ribbon by adding ymin
#' 
#' @param signline a data.frame with x and y coordinates for a sigmoidal line
#' @param height The desired height of the ribbon
#' @param from Whether the line should be treated as the "top", "bot"(tom), or "mid"(dle) of the ribbon
#' 
#' @return a data.frame with coordinates x, y, and ymin for the ribbon.
sigribbon <- function(sigline, height, from = "top") {
  library(dplyr)
  
  if (!is.numeric(height)) 
    stop("Your height input should be numeric.")
  
  if (height == Inf | is.na(height))
    stop("Your height input is NA or infinite.")
    
  if (from == "top") {
    ribbon <- sigline %>%
      mutate(ymin = y - height)
  } else if (from == "bot") {
    ribbon <- sigline %>%
      rename(ymin = y) %>%
      mutate(y = ymin + height)
  } else if (from == "mid") {
    ribbon <- sigline %>%
      mutate(y = y + height/2,
             ymin = y - height/2)
  }
  
  ribbon
  
}

#' Make group nodes for categorical data
#' 
#' @param anno sample annotations
#' @param group_by which columns of anno to use for sample grouping (should be at least 2). Group must have "_id", "_label", "_color" columns.
#' @examples 
#' make_group_nodes(tasic_2016_anno, c("primary_type", "secondary_type")
make_group_nodes <- function(anno,
                             group_by,
                             value_col,
                             xpos = NULL) {
  
  
  if (length(group_by) == 1) 
    stop("Group_by should have at least 2 elements.")
  
  library(dplyr)
  
  nodes <- data.frame(id = numeric(),
                      name = character(),
                      color = character(),
                      n = numeric(),
                      group = character(),
                      xpos = numeric())
  
  for (i in 1:length(group_by)) {
    base <- group_by[i]
    anno_id <- paste0(base,"_id")
    parsed_anno_id <- rlang::parse_expr(anno_id)
    anno_label <- paste0(base,"_label")
    anno_color <- paste0(base,"_color")
    
    grouping <- c(anno_id,anno_label,anno_color)
    parsed_grouping <- lapply(grouping, rlang::parse_expr)
    
    group_nodes <- anno %>%
      select(one_of(grouping)) %>%
      group_by(!!!parsed_grouping) %>%
      summarise(n = n()) %>%
      arrange(!!parsed_anno_id) %>%
      mutate(group = base) %>%
      ungroup()
    
    if (is.null(xpos)) {
      group_nodes <- mutate(group_nodes, xpos = i)
    } else {
      group_nodes <- mutate(group_nodes, xpos = xpos[i])
    }
    
    names(group_nodes) <- c("id","name","color","n","group","xpos")
    
    nodes <- rbind(nodes, group_nodes)
  }
  
  nodes
}


#' Make plot nodes
#'
#' @param group_nodes output from make_group_nodes()
#' @param pad the fraction of y-axis space to use as padding between nodes. Default = 0.1.
#' @param width the width of the nodes in plot space. Default = 0.1.
#'
#' @return a data.frame with coordinates compatible with geom_rect() in ggplot2
#' @export
#'
make_plot_nodes <- function(group_nodes,
                            # % of height to distribute for padding between nodes
                            pad = 0.1,
                            # plot space for total width
                            width = 0.1) {
  
  library(dplyr)
  
  total_n <- sum(group_nodes$n)/length(unique(group_nodes$group))
  total_pad <- total_n * pad
  
  group_rects <- group_nodes %>%
    group_by(group) %>%
    mutate(xmin = xpos - width/2,
           xmax = xpos + width/2,
           n_groups = n(),
           group_pad = ifelse(n_groups > 1, 
                              total_pad/(n() - 1), 
                              total_pad),
           n_cum = cumsum(n),
           ymin = ifelse(n_groups > 1,
                         lag(n_cum, default = 0) + (1:n() - 1)*group_pad,
                         group_pad / 2),
           ymax = ifelse(n_groups > 1,
                         n_cum + (0:(n()-1))*group_pad,
                         group_pad /2 + n))
  
  group_rects
  
}


#' Make group links
#'
#' @param anno annotation data
#' @param group_by grouping
#' @param plot_nodes output of make_plot_nodes()
#'
#' @return a data.frame with starts, ends, and group details for generating links
#' @export
#'
make_group_links <- function(anno,
                        group_by,
                        plot_nodes) {
  
  library(dplyr)
  
  if (length(group_by) == 1)
    stop("Group_by should have at least 2 elements.")
  
  pairs <- list()
  
  for (i in 2:length(group_by)) {
    pair <- group_by[(i - 1):i]
    pairs <- c(pairs, list(pair))
  }
  
  for (pair in pairs) {
    base <- pair
    anno_id <- paste0(base,"_id")
    parsed_anno_id <- lapply(anno_id, rlang::parse_expr)
    anno_label <- paste0(base,"_label")
    anno_color <- paste0(base,"_color")
    
    group1_nodes <- plot_nodes %>%
      filter(group == base[1]) %>%
      select(group, id, xmax, ymin, n)
    names(group1_nodes) <- c("group1",anno_id[1],"x","group1_min", "group1_n")
    
    group2_nodes <- plot_nodes %>%
      filter(group == base[2]) %>%
      select(group, id, xmin, ymin, n)
    names(group2_nodes) <- c("group2",anno_id[2],"xend","group2_min", "group2_n")
    
    grouping <- c(anno_id, anno_label, anno_color)
    parsed_grouping <- lapply(grouping, rlang::parse_expr)
    
    group_links <- anno %>%
      select(one_of(grouping)) %>%
      group_by(!!!parsed_grouping) %>%
      summarise(n = n()) %>%
      arrange(!!!parsed_anno_id) %>%
      mutate(group1 = base[1],
             group2 = base[2]) %>%
      ungroup() %>%
      left_join(group1_nodes, by = c(anno_id[1], "group1")) %>%
      left_join(group2_nodes, by = c(anno_id[2], "group2")) %>%
      group_by(!!parsed_anno_id[[1]]) %>%
      arrange(!!parsed_anno_id[[2]]) %>%
      mutate(y = group1_min + cumsum(n)) %>%
      ungroup() %>%
      group_by(!!parsed_anno_id[[2]]) %>%
      arrange(!!parsed_anno_id[[1]]) %>%
      mutate(yend = group2_min + cumsum(n)) %>%
      ungroup()
      
    names(group_links) <- c("group1_id","group2_id",
                            "group1_label","group2_label",
                            "group1_color","group2_color",
                            "n","group1","group2",
                            "x","group1_min","group1_n","xend","group2_min","group2_n",
                            "y","yend")
    
    group_links <- group_links %>%
      rowwise() %>%
      mutate(link_id = paste0(group1_label,"_",group1_id,"_to_",
                              group2_label,"_",group2_id)) %>%
      ungroup() %>%
      mutate(group1_frac = n/group1_n,
             group2_frac = n/group2_n)
    
    if (exists("all_links")) {
      all_links <- rbind(all_links, group_links)
    } else {
      all_links <- group_links
    }
    
  }
  
  all_links
  
}



#' Generate links in plot space between nodes for a river plot
#'
#' @param group_links output of make_group_links()
#' @param fill the group to use for fill colors
#'
#' @return a data.frame compatible with geom_ribbon() in ggplot2.
#' @export
#'
make_plot_links <- function(group_links,
                            fill = NULL) {
  
  library(dplyr)
  
  for (i in 1:nrow(group_links)) {
    
    link_line <- sigline(x = group_links$x[i], xend = group_links$xend[i],
                         y = group_links$y[i], yend = group_links$yend[i])
    
    link_ribbon <- sigribbon(link_line, h = group_links$n[i])
    
    if (is.null(fill)) {
      link_ribbon <- mutate(link_ribbon, fill = "#808080")
    } else if (fill == group_links$group1[i]) {
      link_ribbon <- mutate(link_ribbon, fill = group_links$group1_color[i])
    } else if (fill == group_links$group2[i]) {
      link_ribbon <- mutate(link_ribbon, fill = group_links$group2_color[i])
    } else {
      link_ribbon <- mutate(link_ribbon, fill = "#808080")
    }
    
    link_ribbon <- mutate(link_ribbon, link_id = i)
    
    if (exists("all_ribbons")) {
      all_ribbons <- rbind(all_ribbons, link_ribbon)
    } else {
      all_ribbons <- link_ribbon
    }
  }
  
  all_ribbons
  
}

#' Build river plots for annotations
#' 
#' @param anno The sample annotations to use. Must have _id, _label, and _color columns for each grouping.
#' @param group_by The bases to use for the river plot, from left to right.
#' @param show_labels Logical, whether or not to show labels. Default is TRUE.
#' @param label_pos Label position - "left", "center", or "right". Can be specified for each entry in group_by. Default = "center".
#' @param min_link_size Numeric, the minimum fraction of cells in either group that must be included in a link for display. Default is 0 (show all links).
#' @param pad The fraction of vertical space to use as padding between groups. Default = 0.1.
#' @param fill_group One group to use as a source for ribbon colors. Default is NULL.
#' 
#' @return A ggplot2 plot object.
#' 
build_river_plot <- function(anno, 
                             grouping, 
                             show_labels = TRUE,
                             label_pos = "center",
                             min_link_size = 0,
                             pad = 0.1, 
                             fill_group = NULL) {
  
  library(dplyr)
  library(ggplot2)
  
  group_nodes <- make_group_nodes(anno, grouping)
  plot_nodes <- make_plot_nodes(group_nodes, pad = pad)
  
  if (show_labels) {
    if (length(label_pos) == 1) {
      label_pos <- rep(label_pos, length(grouping))
    }
    
    align_df <- data.frame(group = grouping,
                           label_pos = label_pos)
    
    node_labels <- plot_nodes %>%
      left_join(align_df, by = "group") %>%
      group_by(group, name) %>%
      summarise(label_xpos = case_when(label_pos[1] == "left" ~ xpos[1] - 0.075,
                                       label_pos[1] == "center" ~ as.numeric(xpos[1]),
                                       label_pos[1] == "right" ~ xpos[1] + 0.075),
                label_hjust = case_when(label_pos[1] == "left" ~ as.numeric(1),
                                        label_pos[1] == "center" ~ 0.5,
                                        label_pos[1] == "right" ~ as.numeric(0)),
                label_ypos = (ymin[1] + ymax[1]) / 2)
    
  }
  
  group_links <- make_group_links(anno, grouping, plot_nodes) %>%
    filter(group1_frac > min_link_size | group2_frac > min_link_size)
  plot_links <- make_plot_links(group_links, fill = fill_group)
  
  p <- ggplot() +
    geom_rect(data = plot_nodes,
              aes(xmin = xmin, xmax = xmax,
                  ymin = ymin, ymax = ymax,
                  fill = color),
              color = "#808080") +
    geom_ribbon(data = plot_links,
                aes(x = x, ymax = y,
                    ymin = ymin,
                    group = link_id,
                    fill = fill),
                color = "#808080",
                alpha = 0.4) +
    scale_fill_identity() +
    scale_y_reverse() +
    theme_void()
  
  if (show_labels) {
    p <- p +
      geom_text(data = node_labels,
                aes(x = label_xpos,
                    y = label_ypos,
                    hjust = label_hjust,
                    label = name),
                vjust = 0.3)
  }
  
  p
  
}

#' River plot bokeh
#'
#' @param anno The sample annotations to use. Must have _id, _label, and _color columns for each grouping.
#' @param group_by The bases to use for the river plot, from left to right.
#' @param pad The fraction of vertical space to use as padding between groups. Default = 0.1.
#' @param fill_group One group to use as a source for ribbon colors. Default is NULL.
#'
#' @return an interactive rbokeh plot object.
#' @export
#'
build_river_plot_bokeh <- function(anno, group_by, pad = 0.1, fill_group = NULL) {
  library(rbokeh)
  
  group_nodes <- make_group_nodes(anno, group_by)
  plot_nodes <- make_plot_nodes(group_nodes, pad = pad)
  
  #must be sorted by link_id for polygon hover to work
  group_links <- make_group_links(anno, group_by, plot_nodes) %>%
    arrange(link_id)
  
  poly_links <- data.frame(x = numeric, y = numeric, link_id = character())
  for (i in 1:nrow(group_links)) {
    plot_links <- make_plot_links(group_links[i,], fill = fill_group)
    
    poly_link <- data.frame(x = c(rev(plot_links$x),plot_links$x),
                            y = c(rev(plot_links$ymin),plot_links$y),
                            link_id = rep(plot_links$link_id, 2))
    
    poly_links <- rbind(poly_links, poly_link)
      
  }
  
  # hover behavior is strange - it doesn't follow grouping,
  # so these have to be added in sequence
  poly_links$group1_label <- group_links$group1_label
  poly_links$group2_label <- group_links$group2_label
  poly_links$n <- group_links$n
  
  b <- figure(height = 750, width = 1000) %>%
    ly_rect(data = plot_nodes,
            xleft = xmin, xright = xmax,
            ybottom = -ymin, ytop = -ymax,
            color = color,
            hover = list("N Cells" = n),
            fill_alpha = 1) %>%
    ly_polygons(data = poly_links,
                xs = x, ys = -y,
                group = link_id,
                # hover = list("Group 1" = group1_label,
                #              "Group 2" = group2_label,
                #              "N Cells" = n),
                color = "#808080")
  
  b
  
}
