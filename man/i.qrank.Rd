% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/i.qrank.R
\name{i.qrank}
\alias{i.qrank}
\title{Rank Individual Time Series According to a Given Timewise Quantile Series}
\usage{
i.qrank(x, prob = 0.5)
}
\arguments{
\item{x}{T by k data matrix: T data points in rows with each row being data at a given time point,
and k time series in columns.}

\item{prob}{Probability, the quantile series of which is to be computed.
Default value is 0.5.}
}
\value{
A list containing:
\itemize{
\item standardized - A matrix containing standardized time series.
\item qts - The timewise quantile of order prob.
\item ranks - Rank of the individual time series according to a the given timewise quantile series.
\item crit - Sum of absolute deviations of each individual series. Distance of each series to the quantile.
}
}
\description{
Use sum of absolute deviations to select the individual time series that is closest
to a given timewise quantile series.
}
\examples{
data(TaiwanAirBox032017)
output <-  i.qrank(TaiwanAirBox032017[,1:3])

}
