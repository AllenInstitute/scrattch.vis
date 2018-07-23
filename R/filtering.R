filter_from_list <- function(df, 
                             filter_list, 
                             match_type = "exact") {
  
  filter_names <- names(filter_list)
  missing_cols <- setdiff(filter_names, names(df))
  
  filter_list <- filter_list[filter_names %in% names(df)]
  filter_names <- filter_names[filter_names %in% names(df)]
  
  if(length(match_type) > 1) {
    if(length(match_type) >= length(filter_names)) {
      match_type <- match_type[filter_names %in% names(df)]
    } else {
      match_logic <- filter_names %in% names(df)
      match_logic <- match_logic[1:length(filter_names)]
      missing_type <- rep("exact", length(filter_names) - length(match_type))
      match_type <- c(match_type[match_logic], missing_type)
    }
  } else {
    match_type <- rep(match_type, length(filter_names))
  }
  
  if(length(filter_names) == 0) {
    message("No columns of df match the names of the filter_list. Nothing filtered.")
    return(df)
  } else {
    if(length(missing_cols) > 0) {
      missing_cols <- paste(missing_cols, collapse = ",")
      message(paste0("No match found for ", missing_cols,". Filtering for other matches."))
    }
    
    
    
    for(i in seq_along(filter_names)) {
      if(match_type[i] == "exact") {
        df <- df[df[[filter_names[i]]] %in% filter_list[[i]],]
        filter_message <- paste(nrow(df),"matches for",filter_names[i])
        message(filter_message)
      } else if(match_type[i] == "grep") {
        match_string <- paste(filter_list[[i]], collapse = "|")
        df <- df[grepl(match_string, df[[filter_names[i]]]),]
        filter_message <- paste(nrow(df),"matches for",filter_names[i])
        message(filter_message)
        
      }
    }
    
    return(df)
    
  }
  
}