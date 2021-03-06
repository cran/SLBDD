% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summaryoutliers.R
\name{SummaryOutliers}
\alias{SummaryOutliers}
\title{Summary Outliers}
\usage{
SummaryOutliers(
  x,
  type = c("LS", "AO", "TC"),
  tsmethod = "arima",
  args.tsmethod = list(order = c(5, 0, 0))
)
}
\arguments{
\item{x}{T by k data matrix: T data points in rows with each row being data at a given time point,
and k time series in columns.}

\item{type}{A character vector indicating the type of outlier to be considered by the
detection procedure. See 'types' in tso function.}

\item{tsmethod}{The framework for time series modeling. Default is "arima". See 'tsmethod' in tso function.}

\item{args.tsmethod}{An optional list containing arguments to be passed to the function invoking the
method selected in tsmethod. See 'args.tsmethod' in tso function. Default value is c(5,0,0).}
}
\value{
A list containing:
\itemize{
\item Otable - Summary of various types of outliers detected.
\item x.cleaned - Outlier-adjusted data.
\item xadja - T-dimensional vector containing the number of time series that have outlier
at a given time point.
}
}
\description{
Use the command "tso" of the R package "tsoutliers" to identify outliers
for each individual time series.
}
\examples{
data(TaiwanAirBox032017)
output <- SummaryOutliers(TaiwanAirBox032017[1:50,1:3])
}
