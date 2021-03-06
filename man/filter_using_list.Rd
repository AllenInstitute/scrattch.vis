% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filtering.R
\name{filter_using_list}
\alias{filter_using_list}
\title{Filter a data.frame with parameters provided in a list}
\usage{
filter_using_list(df, filter_list, verbose = TRUE)
}
\arguments{
\item{df}{the data.frame to filter}

\item{filter_list}{a list of lists, named for target columns to filter, and containing filter parameters: values, and match_type. 
Match types can be "exact" or "grep" for text; 
 "lt", "lteq", "gt", "gteq", or "eq" for numbers.}
}
\value{
a filtered data.frame
}
\description{
Filter a data.frame with parameters provided in a list
}
\examples{
library(tasic2016data)

anno <- tasic_2016_anno

filter_list <- list(pass_qc_checks = list(values = "Y",
                                          match_type = "exact"),
                    primary_type = list(values = c("Pvalb","Vip"),
                                        match_type = "grep"))

filtered_anno <- anno \%>\%
  filter_using_list(filter_list)

}
