#' Dendrogram plots for grouped samples with annotation barplots
#' 
#' @param anno Sample annotations. The first column should be sample_name, and each annotation should have \_id, \_label, and \_color columns
#' @param dend Cluster/group dendrogram. The labels of the dendrogram should match one of the \_label columns in anno.
#' @param dend_group The base of the group that matches dendrogram to anno (e.g. "cluster" if "cluster_label" matches anno).
#' @param leaf_labels Logical, whether or not to plot the leaf labels. Default is TRUE.
#' @param font_size Numeric, the size of the leaf labels to plot. Default is 7.
#' @param grouping The base(s) of groups to use for generating barplots (e.g. "cluster" for "cluster_label" in anno). Default is NULL.
#' @param count_bars Logical, whether or not to display a barplot showing the size of each dend_group.
#' @param count_cutoff Numeric, the maximum size to display in the count_bars. Clusters with > count_cutoff will have truncated bars.
#' @param count_log Logical, whether or not to log scale the count_bars. Helpful for showing small groups.
#' @param panel_pad Numeric, the fraction of the dendrogram/bar line to use as padding between panels.
#' 
#' @return a ggplot2 plot object
#' 
#' @export
#' 
group_dend_bar_plot <- function(anno,
                                dend,
                                dend_group,
                                leaf_labels = TRUE,
                                font_size = 7,
                                grouping = NULL,
                                count_bars = TRUE,
                                count_cutoff = NULL,
                                count_log = FALSE,
                                panel_pad = 0.05) {
  
  if(!paste0(dend_group,"_id") %in% names(anno)) {
    stop("dend_group error:",paste0(dend_group,"_id"), "not found in anno.")
  }
  
  dend_gg <- dendextend::as.ggdend(dend)
  
  dend_group_id <- rlang::sym(paste0(dend_group,"_id"))
  dend_group_label <- rlang::sym(paste0(dend_group,"_label"))
  dend_group_color <- rlang::sym(paste0(dend_group,"_color"))
  
  dend_seg <- dend_gg$segments
  dend_ymax <- max(dend_seg$y)
  dend_seg$y <- dend_seg$y / dend_ymax
  dend_seg$yend <- dend_seg$yend / dend_ymax
  
  n_dend_groups <- length(labels(dend))
  
  p <- ggplot2::ggplot() +
    ggplot2::scale_size(range = c(0.5, 1)) +
    ggplot2::scale_color_identity() +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_y_continuous(expand = c(0.25,0)) +
    ggplot2::scale_x_continuous(limits = c(-1, n_dend_groups + 1)) +
    ggplot2::theme_void()
  
  p <- p +
    ggplot2::geom_segment(data = dend_seg,
                          ggplot2::aes(x = x,
                                       xend = xend,
                                       y = y,
                                       yend = yend),
                          lineend = "square")
  
  anno <- anno %>%
    dplyr::filter(!!dend_group_label %in% labels(dend))
  
  if(!is.null(grouping)) {
    grouping_syms <- purrr::map(
      grouping,
      function(x) {
        group_id <- rlang::sym(paste0(x,"_id"))
        group_label <- rlang::sym(paste0(x,"_label"))
        group_color <- rlang::sym(paste0(x,"_color"))
        
        if(!paste0(x,"_id") %in% names(anno)) {
          stop("grouping error:",paste0(x,"_id"), "not found in anno.")
        }
        
        list(group_id = group_id,
             group_label = group_label,
             group_color = group_color)
      }
    )
    
    group_rects <- list()
    
    for(i in 1:length(grouping)) {
      group_syms <- grouping_syms[[i]]
      
      rects <- anno %>%
        select(!!dend_group_id, !!dend_group_label, !!dend_group_color, 
               !!group_syms$group_id, !!group_syms$group_label, !!group_syms$group_color) %>%
        group_by(!!dend_group_id, !!group_syms$group_id, !!group_syms$group_label) %>%
        mutate(dg_n = n()) %>%
        ungroup() %>%
        group_by(!!dend_group_id) %>%
        arrange(!!group_syms$group_id) %>%
        mutate(dend_group_n = n(),
               dg_frac = dg_n / dend_group_n) %>%
        unique() %>%
        arrange(!!group_syms$group_id) %>%
        mutate(dg_cum_frac = cumsum(dg_frac)) %>%
        ungroup() %>%
        arrange(!!dend_group_id, !!group_syms$group_id) %>%
        group_by(!!dend_group_id) %>%
        mutate(xmin = !!dend_group_id - 0.5,
               xmax = !!dend_group_id + 0.5,
               ymax = -1 * (i - 1) - lag(dg_cum_frac, default = 0) - panel_pad * i,
               ymin = -1 * (i - 1) - dg_cum_frac - panel_pad * i)
      
      group_rects[[i]] <- rects
      
      p <- p +
        ggplot2::geom_rect(data = rects,
                           ggplot2::aes(xmin = xmin,
                                        xmax = xmax,
                                        ymin = ymin,
                                        ymax = ymax,
                                        fill = !!group_syms$group_color))
    }
    
    
  }
  
  if(count_bars) {
    if(!is.null(grouping)) {
      y_adj <- length(grouping)
    } else {
      y_adj <- 0
    }
    
    if(count_log) {
      count_rects <- anno  %>%
        group_by(!!dend_group_id, !!dend_group_color, !!dend_group_label) %>%
        summarise(dg_n = log10(n() + 1)) %>%
        ungroup()
      
      if(!is.null(count_cutoff)) {
        count_cutoff <- log10(count_cutoff + 1)
      }
      
    } else {
      count_rects <- anno  %>%
        group_by(!!dend_group_id, !!dend_group_color, !!dend_group_label) %>%
        summarise(dg_n = n()) %>%
        ungroup() 
    }
    
    if(!is.null(count_cutoff)) {
      count_rects <- count_rects %>%
        mutate(adj_n = ifelse(dg_n <= count_cutoff, 
                              dg_n, 
                              count_cutoff)) %>%
        mutate(f = adj_n / count_cutoff)
      
    } else {
      count_rects <- count_rects %>%
        mutate(f = dg_n / max(dg_n)) 
    }
    
    count_rects <- count_rects %>%
      mutate(xmin = !!dend_group_id - 0.5,
             xmax = !!dend_group_id + 0.5,
             ymin = -1 * y_adj - f - panel_pad * (y_adj + 1),
             ymax = -1 * y_adj - panel_pad * (y_adj + 1))
    
    p <- p +
      ggplot2::geom_rect(data = count_rects,
                         ggplot2::aes(xmin = xmin,
                                      xmax = xmax,
                                      ymin = ymin,
                                      ymax = ymax,
                                      fill = !!dend_group_color))
    
  }
  
  if(leaf_labels) {
    y_adj <- 1
    
    if(!is.null(grouping)) {
      y_adj <- y_adj + length(grouping)
    }
    if(count_bars) {
      y_adj <- y_adj + 1
    }
    
    leaf_anno <- anno %>%
      select(!!dend_group_id, !!dend_group_label, !!dend_group_color) %>%
      unique() %>%
      mutate(y = -1 * (y_adj - 1) - (panel_pad * (y_adj)))
    
    p <- p +
      ggplot2::geom_text(data = leaf_anno,
                         ggplot2::aes(x = !!dend_group_id,
                                      y = y,
                                      label = !!dend_group_label,
                                      color = !!dend_group_color),
                         angle = 90,
                         hjust = 1,
                         vjust = 0.3,
                         size = pt2mm(font_size))
  }
  
  return(p)
}
