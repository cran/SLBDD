% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DLdata.R
\name{DLdata}
\alias{DLdata}
\title{Create an input data matrix for a Deep learning program that uses time series data.}
\usage{
DLdata(x, forerate = 0.2, locY = 1, lag = 1)
}
\arguments{
\item{x}{T by k data matrix: T data points in rows and k time series in columns.}

\item{forerate}{Fraction of sample size to form the forecasting (or testing) sample.}

\item{locY}{Locator for the dependent variable.}

\item{lag}{Number of lags to be used to form predictors.}
}
\value{
A list containing:
\itemize{
\item Xtrain - Standardized predictors matrix.
\item Ytrain - Dependent variable in training sample.
\item Xtest - Predictor in testing sample, standardized according to X_train.
\item Ytest - Dependent variable in the testing sample.
\item nfore - Number of forecasts.
}
}
\description{
R command to setup the training and forecasting data for deep learning.
}
\examples{
x <- matrix(rnorm(7000), nrow=700, ncol=100)
m1 <- DLdata(x, forerate=c(200/nrow(x)), lag=6, locY=6)
}
