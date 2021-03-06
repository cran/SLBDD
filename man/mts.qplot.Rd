% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mts.qplot.R
\name{mts.qplot}
\alias{mts.qplot}
\title{Plot Timewise Quantiles in One Frame}
\usage{
mts.qplot(
  x,
  title = "mts quantile plot",
  prob = c(0.25, 0.5, 0.75),
  scaling = TRUE,
  xtime = NULL,
  plot = TRUE
)
}
\arguments{
\item{x}{T by k data matrix: T data points in rows with each row being data at a given time point,
and k time series in columns.}

\item{title}{Character with the title of the plot. Default title is "mts quantile plot".}

\item{prob}{Probability, the quantile series of which is to be computed.
Default values are 0.25, 0.5, 0.75.}

\item{scaling}{If scaling = TRUE (default), then each series is standardized based on its own range.
If scaling = FALSE, then the original series is used.}

\item{xtime}{A vector with the values for the x labels. Default values are 1, 2, 3, ...}

\item{plot}{Receives TRUE or FALSE values. If the value is TRUE, a quantile plot is presented.
Defaults is TRUE.}
}
\value{
A list containing:
\itemize{
\item standardized - A matrix containing standardized time series.
\item qseries - Matrix of timewise quantiles series of order prob.
}
}
\description{
Plot timewise quantiles in one frame.
}
\examples{
data(TaiwanAirBox032017)
output <-  mts.qplot(TaiwanAirBox032017[,1:5])

}
