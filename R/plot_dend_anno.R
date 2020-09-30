#' dend_anno_barplot
#' 
#' 
#' @param anno Sample annotations. The first column should be sample_name, and each annotation should have \_id, \_label, and \_color columns. Requires cluster_id which needs to be sequential in order of the dendrogram.
#' @param dend Dendrogram object. 
#' @param section_wedges Default is NULL. Use annotation to separate  Can be used to generate lines between sections to divide leaves of dendrogram, e.g. separating subclass, class. 
#' @param bar_variables base name of variables to be represented as bargraphs below dendrogram. Annotation variables need to be represented as \_id, \_label, \_color in anno.
#' @param nr_hist plotting of cluster size below dend. Default is TRUE.
#' @param return_type What values to return - can be "plot", "data", or "both". Default is "plot".
#' @example_data:
#'  
#' load("data/dend_anno_example/dend.rda")
#' load("data/dend_anno_example/anno.df.rda")
#' 
#' 
#' @usage dend.plot <- dend_anno_barplot(anno=anno.df, dend=dend, section_wedges="class", bar_variables=c("sex", "batch"), nr_hist=TRUE)
#' 
#'  
#'    


dend_anno_barplot <- function(anno,
                           dend,
                           section_wedges=NULL,
                           bar_variables=NULL,
                           nr_hist=TRUE,
                           panel_width=0.2
) {
  
  
  # required libraries
  #library(dendextend)
  library(dplyr)
  library(ggplot2)
  
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
      group_by_at(var_wedge) %>%
      summarise(x = min(cluster_id - 0.5),
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
      mutate(xend = x)
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
    grouping_label <- paste0(bar_variables[1], "_label")
    grouping_color <- paste0(bar_variables[1], "_color")
    cluster_id <- "cluster_id"
    # platform rectangles
    rect <- anno %>%
      select(cluster_id, cluster_label, cluster_color, grouping_id,  grouping_label, grouping_color) %>%        
      group_by_(cluster_id, grouping_id,  grouping_label) %>%
      mutate(ly_n = n()) %>%
      ungroup() %>%
      group_by(cluster_id) %>%
      arrange_(grouping_id) %>%
      mutate(cluster_n = n(),
             ly_frac = ly_n/cluster_n) %>%
      unique() %>%
      arrange_(grouping_id) %>%
      mutate(ly_cum_frac = cumsum(ly_frac)) %>%
      ungroup() %>%
      arrange_(cluster_id, grouping_id) %>%
      group_by(cluster_id) %>%
      mutate(xmin = cluster_id - 0.5,
             xmax = cluster_id + 0.5,
             ymax = -offset - lag(ly_cum_frac, default = 0)* panel_width,
             ymin = -offset - ly_cum_frac * panel_width)
    
    var_rects[[1]] <- as.data.frame(rect)  
    
    if(length(bar_variables) >1){
      for(i in 2:length(bar_variables)) {
        
        # add panel padding for each bargraph to be plotted
        offset = offset + panel_width + panel_pad
        rect.offset= c(rect.offset, offset)
        panel_width =0.2
        
        grouping_id <- paste0(bar_variables[i], "_id")
        grouping_label <- paste0(bar_variables[i], "_label")
        grouping_color <- paste0(bar_variables[i], "_color")
        cluster_id <- "cluster_id"
        # platform rectangles
        rect <- anno %>%
          select(cluster_id, cluster_label, cluster_color, grouping_id,  grouping_label, grouping_color) %>%        
          group_by_(cluster_id, grouping_id,  grouping_label) %>%
          mutate(ly_n = n()) %>%
          ungroup() %>%
          group_by(cluster_id) %>%
          arrange_(grouping_id) %>%
          mutate(cluster_n = n(),
                 ly_frac = ly_n/cluster_n) %>%
          unique() %>%
          arrange_(grouping_id) %>%
          mutate(ly_cum_frac = cumsum(ly_frac)) %>%
          ungroup() %>%
          arrange_(cluster_id, grouping_id) %>%
          group_by(cluster_id) %>%
          mutate(xmin = cluster_id - 0.5,
                 xmax = cluster_id + 0.5,
                 ymax = -offset - lag(ly_cum_frac, default = 0)* panel_width,
                 ymin = -offset - ly_cum_frac * panel_width)
        
        var_rects[[i]] <- as.data.frame(rect)  
      }
    }
    
    offset = offset + panel_width + panel_pad
    rect.offset= c(rect.offset, offset)
    panel_width =1.5*panel_width
  }
  
  
  
  ############
  # generating cell nr. histogram data
  ############ 
  
  if(nr_hist==TRUE){
    n_rects <- anno  %>%
      group_by(cluster_id, cluster_color, cluster_label) %>%
      summarise(n = n()) %>%
      ungroup() %>%
      mutate(adj_n = log10(n)) %>%
      mutate(f = adj_n/3) %>%
      mutate(xmin = cluster_id - 0.5,
             xmax = cluster_id + 0.5,
             ymin = -offset - f  * panel_width,
             ymax = -offset)
    
    
    ## still to fix guide labels 
    n_guides <- data.frame(y = seq(-offset - panel_width, -offset, by = 1/5 * panel_width),
                           x = 0.5,
                           xend = n_clusters + 1,
                           label = seq(5, 0, by = -1)) %>%
      mutate(yend = y)
    
    offset = offset + panel_width + panel_pad*2
  }
  
  
  
  dend_leaves <- dend_gg$labels %>%
    mutate(cluster_label = label,
           y = - offset, cex=0.3)
  
  ############
  # plotting of the dendrogram 
  ############
  
  
  # dend segments
  flat_plot <- ggplot() +
    geom_segment(data = dend_seg,
                 aes(x = x,
                     xend = xend,
                     y = y,
                     yend = yend,
                     size = lwd,
                     color = col),
                 lineend = "square") 
  
  
  # different bar plots with metadata
  if(!is.null(bar_variables)) {
    lapply(var_rects, function(df) {
      
      flat_plot <<- flat_plot + geom_rect(data = df,
                                          aes(xmin = xmin,
                                              xmax = xmax,
                                              ymin = ymin,
                                              ymax = ymax,
                                              fill = df[,6]) )
      
    } )
  }
  
  if(nr_hist==TRUE){
    flat_plot <- flat_plot +
      geom_rect(data = n_rects,
                aes(xmin = xmin,
                    xmax = xmax,
                    ymin = ymin,
                    ymax = ymax,
                    fill = cluster_color)) 
  }
  
  # N Cells labels
  flat_plot <- flat_plot +  
    geom_text(data = n_guides,
              aes(x = 0,
                  y = y,
                  label = label),
              size = 2,
              hjust = 1) +
    geom_segment(data = n_guides,
                 aes(x=x,
                     xend = xend,
                     y=y,
                     yend = yend),
                 linetype="dashed") +
    # Leaf Labels
    geom_text(data = dend_leaves,
              aes(x = x,
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
      geom_segment(data = wedge_lines,
                   aes(x = x,
                       xend = xend,
                       y = y,
                       yend = yend)) 
    
  }  
  
  flat_plot <- flat_plot +
    scale_size(range = c(0.5, 1)) +
    scale_color_identity() +
    scale_fill_identity() +
    scale_y_continuous(expand = c(0.25,0)) +
    scale_x_continuous(limits = c(-1,n_clusters + 1)) +
    theme_void()
  
  
  
  
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
#' @example_data:
#'  
#' load("data/dend_anno_example/dend.rda")
#' load("data/dend_anno_example/anno.df.rda")
#' 
#' 
#' @usage dend_dotplot <- dend_anno_dotplot(plot_anno=anno.df, dend=dend, y_group="genotype", c_group="genotype", denom=NULL)
#' 
#'  
#'   




dend_anno_dotplot <- function(plot_anno, 
                              dend, 
                              y_group, 
                              c_group="none", 
                              prefix=NULL, 
                              denom=NULL )  {
  
  #####       
  # dotplot
  #####       
  
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  
  ########
  ## anno dotplot
  ########
  
  cl_id <- paste0(prefix, "cluster_id")
  
  
  x_id <- paste0(prefix, "cluster_id")
  #x_cl <- paste0(prefix, "cl")
  x_label <- paste0(prefix, "cluster_label")
  x_color <- paste0(prefix, "cluster_color")
  x_name <- x_label
  
  
  y_id <- paste0(y_group,"_id")
  y_label <- paste0(y_group,"_label")
  y_color <- paste0(y_group,"_color")
  y_name <- y_group
  
  # plot_anno <- anno
  
  #x_order <- data.frame(id=train.cl.df$cluster_id,xpos = train.cl.df$cluster_id)
  x_order <- data.frame(id=labels(dend), xpos=1:length(labels(dend)) )
  names(x_order)[1] <- x_label
  
  plot_anno <- plot_anno %>% 
    left_join(x_order)
  
  y_order <- data.frame(id = unique(plot_anno[[y_id]])) %>%
    arrange(id) %>%
    mutate(ypos = n():1)
  
  names(y_order)[1] <- y_id
  plot_anno <- plot_anno %>% 
    left_join(y_order)
  
  
  
  if(c_group == "none") {
    
    plot_anno <- plot_anno %>% 
      mutate(point_color = "skyblue")
    
  } else {
    point_color <- paste0(c_group,"_color")
    
    plot_anno <- plot_anno %>% 
      mutate(point_color = plot_anno[[point_color]])
    
  }
  
  
  x_labels <- plot_anno %>%
    select(one_of("xpos",x_label)) %>%
    unique() %>%
    arrange(xpos) %>%
    rename_("x_label" = x_label)
  
  
  y_labels <- plot_anno %>%
    select(one_of("ypos",y_label)) %>%
    unique() %>%
    arrange(ypos) %>%
    rename_("y_label" = y_label)
  
  
  plot_anno <- plot_anno %>%
    group_by(xpos,ypos,point_color) %>%
    summarise(n = n()) %>%
    ungroup()
  
  
  # convert to ggdend
  dend_gg <- dendextend::as.ggdend(dend)
  # extract segments for separate plotting later
  dend_seg <- dend_gg$segments
  
  dend_nodes <- dend_gg[["nodes"]]
  dend_nodes <- dend_nodes[!is.na(dend_nodes$pch),]
  
  dend_seg$y <- dend_seg$y +max(plot_anno$ypos) +1
  dend_seg$yend <- dend_seg$yend +max(plot_anno$ypos) +1
  dend_nodes$y <- dend_nodes$y +max(plot_anno$ypos) +1
  
  
  p <- ggplot() +
    scale_color_identity() +
    theme_bw(14) +
    theme(panel.border = element_blank()) +
    theme(line=element_line(size=0.2))
  
  
  if(!is.null(denom)) {
    if(denom == x_group) {
      # As fraction of all X categories
      
      denom_vals <- plot_anno %>%
        group_by(xpos) %>%
        summarise(denom = sum(n))
      
      plot_anno <- plot_anno %>%
        left_join(denom_vals) %>%
        mutate(frac = n/denom)
      
      p <- p +
        geom_point(data = plot_anno,
                   aes(x = xpos, 
                       y = ypos, 
                       size = frac, 
                       color = point_color))
      
    } else if(denom == y_group) {
      # As fraction of all Y categories
      denom_vals <- plot_anno %>%
        group_by(ypos) %>%
        summarise(denom = sum(n))
      
      plot_anno <- plot_anno %>%
        left_join(denom_vals) %>%
        mutate(frac = n/denom)
      
      p <- p +
        geom_point(data = plot_anno,
                   aes(x = xpos, 
                       y = ypos, 
                       size = frac, 
                       color = point_color))
    } else {
      print(2)
      p <- p +
        geom_point(data = plot_anno,
                   aes(x = xpos, y = ypos, size = n, color = point_color))
      
    } }       else {
      print(1)
      p <- p +
        geom_point(data = plot_anno,
                   aes(x = xpos, y = ypos, size = n, color = point_color))
    }
  
  
  p <- p +
    scale_size_area() +
    scale_x_continuous(x_name, breaks = x_labels$xpos, labels = x_labels$x_label) +
    scale_y_continuous(y_name, breaks = y_labels$ypos, labels = y_labels$y_label) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3))
  
  ########
  ## dendrogram
  ########
  
  p<- p +
    geom_segment(data = dend_seg,
                 aes(x = x,
                     xend = xend,
                     y = y,
                     yend = yend,
                     size = 0.5,
                     color = col),
                 lineend = "square") +
    geom_point(data=dend_nodes,
               aes(x=x,
                   y=y,
                   color=col), 
               shape=dend_nodes$pch,
               cex=1)
  
  

  return(p)
} 