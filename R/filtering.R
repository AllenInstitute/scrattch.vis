#' Filter a data.frame with parameters provided in a list
#' 
#' @param df the data.frame to filter
#' @param filter_list a list of lists, named for target columns to filter, and containing filter parameters: values, and match_type. 
#' Match types can be "exact" or "grep" for text; 
#'  "lt", "lteq", "gt", "gteq", or "eq" for numbers.
#' 
#' @return a filtered data.frame
#' 
#' @examples 
#' library(tasic2016data)
#' 
#' anno <- tasic_2016_anno
#' 
#' filter_list <- list(pass_qc_checks = list(values = "Y",
#'                                           match_type = "exact"),
#'                     primary_type = list(values = c("Pvalb","Vip"),
#'                                         match_type = "grep"))
#' 
#' filtered_anno <- anno %>%
#'   filter_using_list(filter_list)
#' 
filter_using_list <- function(df, 
                             filter_list,
                             verbose = TRUE) {
  
  filter_names <- names(filter_list)
  missing_cols <- setdiff(filter_names, names(df))
  
  filter_list <- filter_list[filter_names %in% names(df)]
  filter_names <- filter_names[filter_names %in% names(df)]
  
  if(length(filter_names) == 0) {
    message("No columns of df match the names of the filter_list. Nothing filtered.")
    return(df)
  }
  
  if(length(missing_cols) > 0) {
    missing_cols <- paste(missing_cols, collapse = ",")
    message(paste0("No match found for ", missing_cols,". Filtering for other matches."))
  }

  for(i in seq_along(filter_names)) {
    filter_name <- filter_names[i]
    values <- filter_list[[i]]$values
    match_type <- filter_list[[i]]$match_type
    
    if(match_type == "exact") {
      df <- df[df[[filter_name]] %in% values,]
      filter_message <- paste(nrow(df),"matches for",filter_name)
      if(verbose) message(filter_message)
    } else if(match_type == "grep") {
      match_string <- paste(values, collapse = "|")
      df <- df[grepl(match_string, df[[filter_name]]),]
      filter_message <- paste(nrow(df),"matches for",filter_name)
      if(verbose) message(filter_message)
    } else if(match_type == "lt") {
      df <- df[df[[filter_name]] < values,] 
      filter_message <- paste(nrow(df),"matches for",filter_name)
      if(verbose) message(filter_message)
    } else if(match_type == "lteq") {
      df <- df[df[[filter_name]] <= values,] 
      filter_message <- paste(nrow(df),"matches for",filter_name)
      if(verbose) message(filter_message)
    } else if(match_type == "gt") {
      df <- df[df[[filter_name]] > values,]
      filter_message <- paste(nrow(df),"matches for",filter_name)
      if(verbose) message(filter_message)
    } else if(match_type == "gteq") {
      df <- df[df[[filter_name]] >= values,] 
      filter_message <- paste(nrow(df),"matches for",filter_name)
      if(verbose) message(filter_message)
    } else if(match_type == "lt") {
      df <- df[df[[filter_name]] < values,] 
      filter_message <- paste(nrow(df),"matches for",filter_name)
      if(verbose) message(filter_message)
    }
  }
  
  return(df)
  
}