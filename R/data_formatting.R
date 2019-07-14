#' Convert a matrix to a data.frame for plotting
#'
#' @param mat a matrix or sparse matrix dgCMatrix object from the Matrix package. 
#' @param cols_are whether columns are "gene_names" or "sample_names"
#' @return dataframe 
#' @export
#'
#' @examples
#' mat_to_data_df(matrix, cols_are = "gene_names")
#' mat_to_data_df(matrix, cols_are = "sample_names)
mat_to_data_df <- function(mat,
                           cols_are = "gene_names") {

  if (grepl("sample", cols_are) == TRUE) {
    if (class(mat) == "matrix") {
      mat <- t(mat)
    } else {
      mat <- Matrix::t(mat)
    }
  }
  

  df <- cbind(sample_name = rownames(mat),
              as.data.frame(as.matrix(mat)),
              stringsAsFactors = FALSE)
  
  rownames(df) <- NULL
  
  return(df)
  
}


#' Melt a data_df to prepare for plot parameters
#' 
#' @param data_df a data.frame with sample_name or grouping as well as value columns
#' @param grouping grouping for the samples. Default = "sample_name", but if this is output from group_stats, provide the same grouping.
#' @param value_cols The value columns to include. If NULL (default), will automatically choose all non-grouping columns.
#' @return melted dataframe
#' @examples 
#' melt_data_df(df,grouping = "sample_name", value_cols = NULL)
#' @export
melt_data_df <- function(data_df, 
                         grouping = "sample_name", 
                         value_cols = NULL) {
  
  if (is.null(value_cols)) {
    value_cols <- names(data_df)[!names(data_df) %in% grouping]
  }
  
  melted <- reshape2::melt(data_df, 
                           id.vars = grouping, 
                           measure.vars = value_cols)
  
  names(melted)[names(melted) == "variable"] <- "gene_name"
  
  melted$sample_name <- as.character(melted$sample_name)
  melted$gene_name <- as.character(melted$gene_name)
  
  return(melted)
  
}

#' Compute stats for samples grouped by one or more annotations
#' 
#' @param data_df a data.frame with sample_name and expression values
#' @param value_cols The values to use for computation. If NULL (default), will select all columns in data_df except sample_name.
#' @param anno a data.frame with sample_name and sample annotations
#' @param grouping one or more column names of anno to use for sample grouping
#' @param stat Which statistic to compute for grouped samples. \cr
#' options are: 
#' \itemize{
#'   \item "median"
#'   \item "mean"
#'   \item "tmean" (25\% trimmed mean)
#'   \item "prop_gt0" (proportion of samples > 0)
#'   \item "prop_gt1" (proportion of samples > 1)
#'   \item "min"
#'   \item "max"
#'   }
#' @examples group_stats(melt_df, value_cols = NULL, anno, grouping = "mouse_line", stat = "mean")
group_stats <- function(data_df,
                        value_cols = NULL,
                        anno,
                        grouping,
                        stat = c("median","mean","tmean","prop_gt0","prop_gt1","min","max")) {

  if (is.character(grouping) == TRUE) {
    grouping_quo <- rlang::syms(grouping)
  } else {
    grouping_quo <- dplyr::quo(grouping)
  }
  
  if (is.null(value_cols) == TRUE) {
    value_cols <- names(data_df)[-1]
  }
  
  data_df <- data_df %>%
    dplyr::select(dplyr::one_of(c("sample_name", value_cols)))
  
  anno_data_df <- anno %>%
    dplyr::select(dplyr::one_of(c("sample_name", grouping))) %>%
    dplyr::left_join(data_df, by = "sample_name")
  
  results_df <- anno_data_df %>%
    dplyr::select(dplyr::one_of(grouping)) %>%
    unique()
  
  for (v in value_cols) {
    val_col <- rlang::sym(v)
    val_df <- anno_data_df %>%
      dplyr::group_by(!!!grouping_quo) %>%
      dplyr::summarise(val_stat = text_stat(!!val_col, stat))
    names(val_df)[names(val_df) == "val_stat"] <- v
    results_df <- dplyr::left_join(results_df, val_df, by = grouping)
  }
  
  rownames(results_df) <- NULL
  
  return(results_df)
}



#' Convert expression data to heatmap colors for plotting
#' 
#' @param df data.frame with sample_name and expression values
#' @param value_cols The value columns to convert. Default is (null), which will use all columns except sample_name
#' @param colorset A set of colors to use for the heatmap palette. NULL will use the default for values_to_colors()
#' @param scale color value scaling, passed to values_to_colors(). Default = "linear".
#' @param per_col Logical, whether or not to normalize the colorscale per-value, or across all values
#' @param min_val Minimum value for color scale. If NULL, will be automatically computed from values. Default = 0.
#' @param max_val Maximum value for color scale. If NULL, will be automatically computed from values. Default = 0.
#' 
data_df_to_colors <- function(df,
                              value_cols = NULL,
                              colorset = NULL,
                              scale = "linear",
                              per_col = FALSE,
                              min_val = 0,
                              max_val = NULL) {
  
  library(purrr)
  
  if (is.null(value_cols)) {
    value_cols <- names(df)[-1]
  }
  
  if (scale == "log2") {
    df[[value_cols]] <- log2(df[[value_cols]])
  } else if (scale == "log10") {
    df[[value_cols]] <- log10(df[[value_cols]])
  }
  
  if (is.null(max_val) & per_col == FALSE) {
    max_val <- max(unlist(df[, value_cols]), na.rm = TRUE)
  }
  
  df[,value_cols] <- map(value_cols, 
                         function(x) {
                           vals <- unlist(df[[x]])
                           
                           if (is.null(colorset)) {
                             values_to_colors(vals,
                                              min_val = min_val,
                                              max_val = max_val)
                           } else {
                             values_to_colors(vals,
                                              min_val = min_val,
                                              max_val = max_val,
                                              colorset = colorset)
                           }
                           
                           
                         })
  
  return(df)
  
}


#' Build plot positions for a character vector
#' 
#' @param vec a character vector (like gene names)
#' @param sort how to resort the vector, if necessary. Default is "none", which keeps input order. 
#' Options: "none","rev","random","alpha".
#' @param axis_name which axis to use. Default is "y"
#' 
#' @return a data.frame with columns for gene_name and axis values
build_vec_pos <- function(vec,
                          vec_name = "gene_name",
                          sort = "none",
                          axis_name = "y") {
  
  if(sort == "rev") {
    vec <- rev(vec)
  } else if (sort == "random") {
    vec <- sample(vec, size = length(vec), replace = FALSE)
  } else if (sort == "alpha") {
    vec <- vec[order(vec)]
  }
  
  pos_df <- data.frame(vec_name = vec,
                       pos = 1:length(vec),
                       stringsAsFactors = FALSE)
  
  names(pos_df) <- c(vec_name, axis_name)

  return(pos_df)
}
