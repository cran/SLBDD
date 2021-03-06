% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/outlier.plot.R
\name{outlier.plot}
\alias{outlier.plot}
\title{Find Outliers Using an Upper and a Lower Timewise Quantile Series}
\usage{
outlier.plot(x, prob = 0.05, percent = 0.05, xtime = NULL)
}
\arguments{
\item{x}{T by k data matrix: T data points in rows with each row being data at a given time point,
and k time series in columns.}

\item{prob}{Tail probability. That is, the two quantile series is (prob, 1-prob).
prob is restricted to be in (0,0.15). Default value is 0.05.}

\item{percent}{The number of possible outliers in each side is T*k*prob*percent.}

\item{xtime}{A vector with the values for the x labels. Default values are 1, 2, 3, ...}
}
\value{
A list containing:
\itemize{
\item standardized - A matrix containing standardized time series.
\item qts - The timewise quantile of order prob.
\item minseries - The timewise minimum of the standardized time series.
\item maxseries - The timewise maximum of the standardized time series.
}
}
\description{
Use an upper and a lower timewise quantile series to highlight the possible outliers
in a collection of  time series.
}
\examples{
data(TaiwanAirBox032017)
output <-  outlier.plot(TaiwanAirBox032017[,1:3])

}
