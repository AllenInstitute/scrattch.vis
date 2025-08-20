#' dend_anno_barplot
#' 
#' 
#' @param anno Sample annotations. The first column should be sample_name, and each annotation should have \_id, \_label, and \_color columns. Requires cluster_id which needs to be sequential in order of the dendrogram.
#' @param dend Dendrogram object. 
#' @param section_wedges Default is NULL. Use annotation to separate  Can be used to generate lines between sections to divide leaves of dendrogram, e.g. separating subclass, class. 
#' @param bar_variables base name of variables to be represented as bargraphs below dendrogram. Annotation variables need to be represented as \_id, \_label, \_color in anno.
#' @param nr_hist plotting of cluster size below dend. Default is TRUE.
#' @param return_type What values to return - can be "plot", "data", or "both". Default is "plot".
#' @examples
#' # First, load the example data
#' load("data/dend_anno_example/dend.rda")
#' load("data/dend_anno_example/anno.df.rda")
#'
#' # Then, run the function with the loaded data
#' dend.plot <- dend_anno_barplot(
#'   anno = anno.df,
#'   dend = dend,
#'   section_wedges = "class",
#'   bar_variables = c("sex", "batch"),
#'   nr_hist = TRUE
#' )
dend_anno_barplot <- function(anno,
                              dend,
                              section_wedges=NULL,
                              bar_variables=NULL,
                              nr_hist=TRUE,
                              panel_width=0.2
) {
  
  rect.offset= c()
  panel_pad <- 0.02
  
  # cluster_id order doesn't match dendrogram order
  n_clusters <- max(anno$cluster_id)
  
  # convert to ggdend
  dend_gg <- dendextend::as.ggdend(dend)
  # extract segments for separate plotting later
  dend_seg <- dend_gg$segments
  
  
  if(!is.null(section_wedges)){
    var_wedge <- paste0(section_wedges,"_id")
    sections <- anno %>%
      dplyr::group_by_at(var_wedge) %>%
      dplyr::summarise(x = min(cluster_id - 0.5),
                xend = n_clusters + 0.5)
    
    sections = sections %>% arrange(x)
    tmp = sections$x
    section.length = c(tmp[-1], max(sections$xend)) - tmp
    select = which(section.length >= 3)
    select = sort(union(select, select+1))
    sections = sections[select,]
    
    wedge_lines <- data.frame(x = unique(c(sections$x, sections$xend)),
                              y = 0,
                              yend = -1.8) %>%
      dplyr::mutate(xend = x)
  }
  
  # set padding of initial graph below dend
  offset= panel_pad
  rect.offset= c(rect.offset, offset)
  panel_width = panel_width 
  
  
  ############
  # generating stacked bar graphs for selected metadata
  ############ 
  
  if(!is.null(bar_variables)) {
    var_rects <- list()
    
    grouping_id <- paste0(bar_variables[1], "_id")
    parsed_grouping_id <- rlang::parse_expr(grouping_id)
    grouping_label <- paste0(bar_variables[1], "_label")
    parsed_grouping_label <- rlang::parse_expr(grouping_label)
    grouping_color <- paste0(bar_variables[1], "_color")
    
    # platform rectangles
    rect <- anno %>%
      dplyr::select(cluster_id, cluster_label, cluster_color, grouping_id,  grouping_label, grouping_color) %>%        
      dplyr::group_by(!!parsed_cluster_id, !!parsed_grouping_id, !!parsed_grouping_label) %>%
      dplyr::mutate(ly_n = dplyr::n()) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(cluster_id) %>%
      dplyr::arrange(!!parsed_grouping_id) %>%
      mutate(cluster_n = dplyr::n(),
             ly_frac = ly_n/cluster_n) %>%
      dplyr::unique() %>%
      dplyr::arrange(!!parsed_grouping_id) %>%
      dplyr::mutate(ly_cum_frac = cumsum(ly_frac)) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(cluster_id, !!parsed_grouping_id) %>%
      dplyr::group_by(cluster_id) %>%
      dplyr::mutate(xmin = cluster_id - 0.5,
                    xmax = cluster_id + 0.5,
                    ymax = -offset - lag(ly_cum_frac, default = 0)* panel_width,
                    ymin = -offset - ly_cum_frac * panel_width)
    
    var_rects[[1]] <- as.data.frame(rect)  
    
    if(length(bar_variables) >1){
      for(i in 2:length(bar_variables)) {
        
        # add panel padding for each bargraph to be plotted
        offset <- offset + panel_width + panel_pad
        rect.offset <- c(rect.offset, offset)
        panel_width <- 0.2
        
        grouping_id <- paste0(bar_variables[i], "_id")
        parsed_grouping_id <- rlang::parse_expr(grouping_id)
        grouping_label <- paste0(bar_variables[i], "_label")
        parsed_grouping_label <- rlang::parse_expr(grouping_label)
        grouping_color <- paste0(bar_variables[i], "_color")
        
        # platform rectangles
        rect <- anno %>%
          dplyr::select(cluster_id, cluster_label, cluster_color, grouping_id,  grouping_label, grouping_color) %>%        
          dplyr::group_by(!!parsed_cluster_id, !!parsed_grouping_id, !!parsed_grouping_label) %>%
          dplyr::mutate(ly_n = dplyr::n()) %>%
          dplyr::ungroup() %>%
          dplyr::group_by(cluster_id) %>%
          dplyr::arrange(!!parsed_grouping_id) %>%
          dplyr::mutate(cluster_n = n(),
                        ly_frac = ly_n/cluster_n) %>%
          dplyr::unique() %>%
          dplyr::arrange(!!parsed_grouping_id) %>%
          dplyr::mutate(ly_cum_frac = cumsum(ly_frac)) %>%
          dplyr::ungroup() %>%
          dplyr::arrange(cluster_id, !!parsed_grouping_id) %>%
          dplyr::group_by(cluster_id) %>%
          dplyr::mutate(xmin = cluster_id - 0.5,
                        xmax = cluster_id + 0.5,
                        ymax = -offset - lag(ly_cum_frac, default = 0)* panel_width,
                        ymin = -offset - ly_cum_frac * panel_width)
        
        var_rects[[i]] <- as.data.frame(rect)  
      }
    }
    
    offset <- offset + panel_width + panel_pad
    rect.offset <- c(rect.offset, offset)
    panel_width <- 1.5 * panel_width
  }
  
  
  
  ############
  # generating cell nr. histogram data
  ############ 
  
  if(nr_hist==TRUE){
    n_rects <- anno  %>%
      dplyr::group_by(cluster_id, cluster_color, cluster_label) %>%
      dplyr::summarise(n = dplyr::n()) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(adj_n = log10(n)) %>%
      dplyr::mutate(f = adj_n/3) %>%
      dplyr::mutate(xmin = cluster_id - 0.5,
                    xmax = cluster_id + 0.5,
                    ymin = -offset - f  * panel_width,
                    ymax = -offset)
    
    
    ## still to fix guide labels 
    n_guides <- data.frame(y = seq(-offset - panel_width, -offset, by = 1/5 * panel_width),
                           x = 0.5,
                           xend = n_clusters + 1,
                           label = seq(5, 0, by = -1)) %>%
      dplyr::mutate(yend = y)
    
    offset = offset + panel_width + panel_pad*2
  }
  
  
  
  dend_leaves <- dend_gg$labels %>%
    dplyr::mutate(cluster_label = label,
                   y = - offset, cex=0.3)
  
  ############
  # plotting of the dendrogram 
  ############
  
  
  # dend segments
  flat_plot <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = dend_seg,
      ggplot2::aes(
        x = x,
        xend = xend,
        y = y,
        yend = yend,
        size = lwd,
        color = col
      ),
      lineend = "square"
    )
  
  # different bar plots with metadata
  if(!is.null(bar_variables)) {
    lapply(var_rects, function(df) {
      
      flat_plot <- flat_plot + 
        ggplot2::geom_rect(data = df,
                 ggplot2::aes(xmin = xmin,
                      xmax = xmax,
                      ymin = ymin,
                      ymax = ymax,
                      fill = df[,6]) )
      
    } )
  }
  
  if(nr_hist==TRUE){
    flat_plot <- flat_plot +
      ggplot2::geom_rect(
        data = n_rects,
        ggplot2::aes(
          xmin = xmin,
          xmax = xmax,
          ymin = ymin,
          ymax = ymax,
          fill = cluster_color
        )
      ) 
  }
  
  # N Cells labels
  flat_plot <- flat_plot +  
    ggplot2::geom_text(
      data = n_guides,
      ggplot2::aes(
        x = 0,
        y = y,
        label = label
        ),
      size = 2,
      hjust = 1
    ) +
    ggplot2::geom_segment(
      data = n_guides,
      ggplot2::aes(
        x=x,
        xend = xend,
        y=y,
        yend = yend),
      linetype="dashed") +
    # Leaf Labels
    ggplot2::geom_text(
      data = dend_leaves,
      ggplot2::aes(
        x = x,
        y = y,
        label = label,
        color = col),
      angle = 90,
      hjust = 1,
      vjust = 0.3,
      size = 2) 
  
  if(!is.null(section_wedges)){ 
    #Vertical class separators
    flat_plot <- flat_plot +
      ggplot2::geom_segment(
        data = wedge_lines,
        ggplot2::aes(
          x = x,
          xend = xend,
          y = y,
          yend = yend)
        ) 
    
  }  
  
  flat_plot <- flat_plot +
    ggplot2::scale_size(range = c(0.5, 1)) +
    ggplot2::scale_color_identity() +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_y_continuous(expand = c(0.25,0)) +
    ggplot2::scale_x_continuous(limits = c(-1,n_clusters + 1)) +
    ggplot2::theme_void()
  
  return(flat_plot)
}


#' dend_anno_dotplot
#' 
#' 
#' @param plot_anno Sample annotations. The first column should be sample_name, and each annotation should have \_id, \_label, and \_color columns. Requires cluster_id which needs to be sequential in order of the dendrogram.
#' @param dend Dendrogram object. Needs to be labeled with cluster_label
#' @param y_group  base name of variable to be represented on y-axis of dotplot below dendrogram. Annotation variable needs to be represented as \_id, \_label, \_color in anno.
#' @param c_group base name of variables to color the dotplot by. Annotation variables need to be represented as \_id, \_label, \_color in anno.
#' @param prefix prefix that groups cluster_id/_label/_color. Used in case of mapping to multiple taxonomies. Default is NULL.
#' @param denom When calculating fraction of either x- or y-axis variables. Input can be "cluster" or y_group. Default is NULL.
#' @examples
#' # First, load the example data
#' load("data/dend_anno_example/dend.rda")
#' load("data/dend_anno_example/anno.df.rda")
#'
#' # Then, run the function with the loaded data
#' dend_dotplot <- dend_anno_dotplot(
#'   plot_anno = anno.df,
#'   dend = dend,
#'   y_group = "genotype",
#'   c_group = "genotype",
#'   denom = NULL
#' )
#' 
dend_anno_dotplot<- function(anno, 
                             dend, 
                             y_group, 
                             c_group="none", 
                             prefix=NULL, 
                             denom=NULL, 
                             xlab.size=2,
                             ylab.size=3,
                             panel_width=0.4,
                             panel_pad=0.1)  {
  
  
  rect.offset= c()
  
  # cluster_id order doesn't match dendrogram order
  n_clusters <- length(labels(dend))
  
  
  #########################
  ## build dend elements ##
  #########################
  
  ## dend is always 1 high?
  
  # convert to ggdend
  dend_gg <- dendextend::as.ggdend(dend)
  # extract segments for separate plotting later
  dend_seg <- dend_gg$segments
  
  dend_nodes <- dend_gg[["nodes"]]
  dend_nodes <- dend_nodes[!is.na(dend_nodes$pch),]
  
  ## build anno dotplot ##
  
  # set padding of dotplot below dend
  offset= 2*panel_pad
  rect.offset= c(rect.offset, offset)
  panel_width = panel_width 
  
  cl_id <- paste0(prefix, "cluster_id")
  
  x_id <- paste0(prefix, "cluster_id")
  x_cl <- paste0(prefix, "cl")
  x_label <- paste0(prefix, "cluster_label")
  x_color <- paste0(prefix, "cluster_color")
  x_name <- x_label
  
  y_id <- paste0(y_group,"_id")
  y_label <- paste0(y_group,"_label")
  y_color <- paste0(y_group,"_color")
  y_name <- y_group
  
  plot_anno <- anno
  
  x_order <- data.frame(id=labels(dend), xpos=1:length(labels(dend)) )
  names(x_order)[1] <- x_label
  
  plot_anno <- plot_anno %>% 
    dplyr::left_join(x_order)
  
  y_order <- data.frame(id = unique(plot_anno[[y_id]])) %>%
    dplyr::arrange(id) %>%
    dplyr::mutate(new_id=1:length(id)) %>%
    dplyr::mutate(ly_n=dplyr::n()) %>%
    dplyr::mutate(ly_frac = new_id/ly_n) %>%
    dplyr::mutate(ypos = -offset - lag(ly_frac, default = 0) *panel_width)
  
  names(y_order)[1] <- y_id
  plot_anno <- plot_anno %>% 
    dplyr::left_join(y_order)
  
  if(c_group == "none") {
    
    plot_anno <- plot_anno %>% 
      dplyr::mutate(point_color = "skyblue")
    
  } else {
    point_color <- paste0(c_group,"_color")
    
    plot_anno <- plot_anno %>% 
      dplyr::mutate(point_color = plot_anno[[point_color]])
    
  }
  
  y_label_name <- c("y_label" = y_label)
  y_labels <- plot_anno %>%
    dplyr::select(dplyr::one_of("ypos",y_label)) %>%
    dplyr::unique() %>%
    dplyr::arrange(ypos) %>%
    dplyr::rename(dplyr::all_of(y_label_name))
  
  
  plot_anno <- plot_anno %>%
    dplyr::group_by(xpos,ypos,point_color) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::ungroup()
  
  ## add grid lines behind
  #grid_lines <-
  
  offset = min(y_labels$ypos)  - 2*panel_pad
  
  dend_leaves <- dend_gg$labels %>%
    dplyr::mutate(cluster_label = label,
                  y =  offset, cex=0.3)
  
  ## Plotting ##
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = dend_seg,
                          ggplot2::aes(x = x,
                                       xend = xend,
                                       y = y,
                                       yend = yend,
                                       size = 0.5,
                                       color = "black"),
                          lineend = "square") +
    ggplot2::geom_point(data=dend_nodes,
                        ggplot2::aes(x=x,
                                     y=y,
                                     color=col), 
                        shape=dend_nodes$pch,
                        cex=1)  +
    ggplot2::geom_segment(data =y_labels, 
                          ggplot2::aes(x= 0,
                                       xend= n_clusters,
                                       y=  ypos,
                                       yend=  ypos,
                                       color="grey95"))+
    ggplot2::geom_segment(data =dend_leaves, 
                          ggplot2::aes(x= x,
                                       xend= x,
                                       y=  dend_nodes$y-0.1,
                                       yend= dend_nodes$y + min(y_labels$ypos),
                                       color="grey95"))
  
  
  
  
  
  if(!is.null(denom)) {
    if(denom == x_group) {
      # As fraction of all X categories
      
      denom_vals <- plot_anno %>%
        dplyr::group_by(xpos) %>%
        dplyr::summarise(denom = sum(n))
      
      plot_anno <- plot_anno %>%
        dplyr::left_join(denom_vals) %>%
        dplyr::mutate(frac = n/denom)
      
      p <- p +
        ggplot2::geom_point(data = plot_anno,
                            ggplot2::aes(x = xpos, 
                                         y = ypos, 
                                         size = frac, 
                                         color = point_color))
      
    } else if(denom == y_group) {
      # As fraction of all Y categories
      denom_vals <- plot_anno %>%
        dplyr::group_by(ypos) %>%
        dplyr::summarise(denom = sum(n))
      
      plot_anno <- plot_anno %>%
        dplyr::left_join(denom_vals) %>%
        dplyr::mutate(frac = n/denom)
      
      p <- p +
        ggplot2::geom_point(data = plot_anno,
                            ggplot2::aes(x = xpos, 
                                         y = ypos, 
                                         size = frac, 
                                         color = point_color))
    } else {
      print(2)
      p <- p +
        ggplot2::geom_point(data = plot_anno,
                            ggplot2::aes(x = xpos, y = ypos, size = n, color = point_color))
      
    } }       else {
      print(1)
      p <- p +
        ggplot2::geom_point(data = plot_anno,
                            ggplot2::aes(x = xpos, y = ypos, size = n, color = point_color))
    }
  
  
  
  p <- p +  
    # Leaf Labels
    ggplot2::geom_text(data = dend_leaves,
                       ggplot2::aes(x = x,
                                    y = y,
                                    label = label,
                                    color = col),
                       angle = 90,
                       hjust = 1,
                       vjust = 0.3,
                       size = xlab.size)  +
    ggplot2::geom_text(data = y_labels,
                       ggplot2::aes(x = -5,
                                    y = ypos,
                                    label = y_label,
                                    color = "black"),
                       hjust = 1,
                       vjust = 0.5,
                       size = ylab.size)
  
  p <- p +
    ggplot2::scale_size_area() +
    ggplot2::scale_color_identity() +
    ggplot2::scale_fill_identity() +
    ggplot2::theme_void()
  
  return(p)
  
  
} 
