#' Check input genes against a vector of gene names
#' 
#' @param genes An input character string (or character vector) of gene symbols
#' @param gene_reference A reference character vector of gene symbols
#' @param result Whether to return a vector of matched or unmatched genes, or a list of both. Options are "matched" (default), "unmatched", or "both".
#' 
check_genes <- function(genes, 
                        gene_reference,
                        result = "matched") {
  
  # Add some input QC
  if (class(genes) != "character")
    stop("Your genes input is not a character string or vector.")
  if (class(gene_reference) != "character")
    stop("Your reference input is not a character vector.")
  
  
  if (length(genes) == 1) {
    raw_genes <- unique(split_text(genes))
  } else {
    raw_genes <- unique(genes)
  }
  
  
  # Convert gene_reference to lowercase and remove "-" and " " for matching
  match_genes <- tolower(gsub("[- .]+","_",gene_reference))
  
  # Remove leading X from Riken genes
  match_genes[grepl("x[0-9]+.+rik",match_genes)] <- sub("^x","",match_genes[grepl("x[0-9]+.+rik",match_genes)])
  
  # for loop will retain the order of the genes.
  good_genes <- character()
  bad_genes <- character()
  
  for (x in raw_genes) {
    this_gene <- tolower(gsub("[- .]+","_", x))
    
    if (this_gene %in% match_genes) {
      good_genes <- c(good_genes, gene_reference[match_genes == this_gene][1])
    } else {
      bad_genes <- c(bad_genes, x)
    }
  }
  
  if (result == "matched") {
    unique(good_genes)
  } else if (result == "not_matched") {
    unique(bad_genes)
  } else if (result == "both") {
    
    list(matched = unique(good_genes),
         not_matched = unique(bad_genes))
    
  }
}



#' Evaluate a character string specifying integer values to a numeric vector
#' 
#' @param in_chr a character string
#' @return a numeric object
#' 
#' @examples 
#' test <- "1:5"
#' chr_to_num(test)
#' 
#' test <- "8,7,45,20"
#' chr_to_num(test)
chr_to_num <- function(in_chr) {
  
  if (length(in_chr) != 1) {
    stop("Your input is not a character string.")
  }
  
  round_c_chr <- paste0("round(c(", in_chr, "),0)")
  parsed_chr <- parse(text = round_c_chr)
  result <- eval(parsed_chr)
  return(result)
}


#' Mix two colors additively in RGB space
#' 
#' @param col1 A hex or R color
#' @param col2 Another hex or R color
#' @return The sum of col1 and col2 as a character hex color (e.g. "#FFFFFF")
#' 
#' @examples
#' color_sum("red","green")
#' 
#' color_sum("#1B9E77","#D95F02")
color_sum <- function(col1,col2) {
  
  rgbmat1 <- col2rgb(col1)/255
  rgbmat2 <- col2rgb(col2)/255
  
  mix <- rgbmat1 + rgbmat2
  
  rgb(mix[1],mix[2],mix[3])
  
}


#' Convert font sizes in pt to mm
#' 
#' @param pt A numeric font size in pt.
#' @return A numeric font size in mm.
#' 
#' @examples
#' pt2mm(12)
#' 
#' ggplot(mtcars) +
#'   geom_text(aes(x = mpg, y = wt, label = rownames(mtcars)),
#'             size = pt2mm(7))
pt2mm <- function(pt) {
  
  if (class(pt) != "numeric")
    stop("Your input should be numeric font size in pt.")
  
  mm <- pt / 2.834645669
  return(mm)
}


#' Convert the case of Riken genes no matter input case
#' 
#' @param in_chr a character vector of Riken gene symbols
#' @return a character vector with correck Riken capitalization
#' 
#' @examples
#' test <- c("6330527o06RiK","A930038C07RIK","a330070k13rik")
#' riken_case(test)
riken_case <- function(in_chr) {
  
  library(stringr)
  if (!stringr::str_detect(in_chr, "[rikRIKrIKRik]"))
    warning("Your input is not a Riken gene.")
  
  upper <- toupper(in_chr)
  result <- sub("RIK","Rik",upper)
  return(result)
}


#' Split a character string by commas, spaces, tabs, and line breaks into a character vector
#' 
#' @param in_string a chr string containing commas, spaces, and/or tabs
#' @return a character vector with each object separated by any combination of commas, spaces, and/or tabs
#' 
#' @examples
#' test <- "Hspa8, Scnn1a,Rbp4    Ptgs2; GeneA:GeneB"
#' split_text(test)
#' 
split_text <- function(in_string) {
  
  if (length(in_string) != 1)
    warning("Input might not be a character string.")
  
  out_chr <- strsplit(in_string,"[,.;: \t\r\n]+")[[1]]
  return(out_chr)
}


#' Convert the case of objects in a character vector to Title Case
#' 
#' @param in_chr a character vector
#' @return a character vector with Each Object In Title Case
#' 
#' @examples
#' test <- c("hspa8","scnn1a","fhqwghads")
#' title_case(test)
title_case <- function(in_chr) {

  lower <- tolower(in_chr)
  s <- strsplit(lower, " ")
  result <- paste(toupper(substring(s, 1,1)), substring(s, 2), sep="")
  
  if (length(result) < 2)
    warning("Make sure that your input is a character vector.")
  
  return(result)
    
}


#' Convert values to colors along a color ramp
#' 
#' @param x a numeric vector to be converted to colors
#' @param min_val a number that's used to set the low end of the color scale (default = 0)
#' @param max_val a number that's used to set the high end of the color scale. If NULL (default), 
#' use the highest value in x
#' @param colorset a set of colors to interpolate between using colorRampPalette 
#' (default = c("darkblue","dodgerblue","gray80","orangered","red"))
#' @param missing_color a color to use for missing (NA) values.
#' @return a character vector of hex color values generated by colorRampPalette. Color values will
#' remain in the same order as x.
values_to_colors <- function(x, 
                             min_val = NULL, 
                             max_val = NULL, 
                             colorset = c("darkblue","dodgerblue","gray80","orange","orangered"),
                             missing_color = "black") {
  
  heat_colors <- colorRampPalette(colorset)(1001)
  
  if (is.null(max_val)) {
    max_val <- max(x, na.rm = T)
  } else {
    x[x > max_val] <- max_val
  }
  if (is.null(min_val)) {
    min_val <- min(x, na.rm = T)
  } else {
    x[x < min_val] <- min_val
  }
  
  if (sum(x == min_val, na.rm = TRUE) == length(x)) {
    colors <- rep(heat_colors[1],length(x))
  } else {
    if (length(x) > 1) {
      if (var(x, na.rm = TRUE) == 0) {
        colors <- rep(heat_colors[500], length(x))
      } else {
        heat_positions <- unlist(round((x - min_val) / (max_val - min_val) * 1000 + 1, 0))
        
        colors <- heat_colors[heat_positions]
      }
    } else {
      colors <- heat_colors[500]
    }
  }
  
  if (!is.null(missing_color)) {
    colors[is.na(colors)] <- rgb(t(col2rgb(missing_color)/255))
  }
  
  colors
}

#' Compute Trimmed Mean of a numeric vector
#' 
#' This is a simple wrapper around mean(x, trim)
#' 
#' @param x A numeric vector
#' @param trim the amount to trim
#' 
tmean <- function(x, trim = 0.25, na.rm = TRUE) {
  mean(x, trim = trim, na.rm = na.rm)
}

#' Compute the proportion/fraction of values > 0 in a numeric vector
#' 
#' @param x A numeric vector
#' 
prop_gt0 <- function(x) {
  sum(x > 0)/length(x)
}

#' Compute the proportion/fraction of values > 1 in a numeric vector
#' 
#' @param x A numeric vector
#' 
prop_gt1 <- function(x) {
  sum(x > 1)/length(x)
}

#' Compute statistics for a numeric vector, selected with a simple text string
#' 
#' This is useful for computation based on parameters passed by other functions.
#' 
#' @param x A numeric vector
#' @param stat The statistic to compute. Options are:
#' \itemize{
#'   \item "median"
#'   \item "mean"
#'   \item "tmean" (25\% trimmed mean)
#'   \item "nzmean" (mean of non-zero values)
#'   \item "nzmedian" (median of non-zero values)
#'   \item "prop_gt0" (proportion of samples > 0)
#'   \item "prop_gt1" (proportion of samples > 1)
#'   \item "prop_gt_cutoff" (proportion of samples > cutoff)
#'   \item "min"
#'   \item "max"
#'   }
#' @param cutoff A cutoff for use in stats calculations. Most ignore this value. Default is 0.
text_stat <- function(x, stat, cutoff = NULL) {
  if(stat == "mean") {
    mean(x, na.rm = TRUE)
  } else if(stat == "tmean") {
    tmean(x, na.rm = TRUE)
  } else if(stat == "nzmean") {
    mean(x[x > 0], na.rm = TRUE)
  } else if(stat == "median") {
    median(x, na.rm = TRUE)
  } else if(stat == "nzmedian") {
    median(x[x > 0], na.rm = TRUE)
  } else if(stat == "prop_gt0") {
    sum(x > 0)/length(x)
  } else if(stat == "prop_gt1") {
    sum(x > 1)/length(x)
  } else if(stat == "prop_gt_cutoff") {
    sum(x > cutoff)/length(x)
  } else if(stat == "min") {
    min(x, na.rm = TRUE)
  } else if(stat == "max") {
    max(x, na.rm = TRUE)
  }
}