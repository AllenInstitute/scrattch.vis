% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_components.R
\name{ggplot_max_vals}
\alias{ggplot_max_vals}
\title{Add max value labels to a ggplot}
\usage{
ggplot_max_vals(p, n_stats, width_stat = "samples", max_val_dfs, font_size)
}
\arguments{
\item{p}{The plot to add max values to}

\item{n_stats}{Tcount statistics made by get_n_stats()}

\item{width_stat}{The count used to determine x-position. In sample plots, use "samples". in group plots, use "groups".}

\item{font_size}{The font size in pt to use for header label plotting.}

\item{max_value_dfs}{The max value list generated by build_max_dfs()}
}
\value{
a ggplot2 object
}
\description{
Add max value labels to a ggplot
}
