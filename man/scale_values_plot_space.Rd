% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_components.R
\name{scale_values_plot_space}
\alias{scale_values_plot_space}
\title{Scale values to plot space}
\usage{
scale_values_plot_space(
  x,
  min_ps,
  max_ps = NULL,
  min_val = 0,
  max_val = NULL,
  extent = 0.9
)
}
\arguments{
\item{x}{A vector of values to scale}

\item{min_ps}{The minimum value in plot space to use}

\item{max_ps}{The maximum value in plot space to use. Default is NULL, which will use min_ps + 1.}

\item{min_val}{The minimum value to plot. Default is 0. If NULL, will use min(x)}

\item{max_val}{The maximum value to plot. Default is NULL, which will use max(x)}

\item{extent}{The fraction of the space between min_ps and max_ps to use. Default is 0.9.}
}
\value{
a vector of scaled values.
}
\description{
Scale values to plot space
}
