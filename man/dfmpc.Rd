% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DFMPC.R
\name{dfmpc}
\alias{dfmpc}
\title{Dynamic Factor Model by Principal Components}
\usage{
dfmpc(x, stand = 0, mth = 4, r, lagk = 0)
}
\arguments{
\item{x}{T by k data matrix: T data points in rows with each row being data at a given time point,
and k time series in columns.}

\item{stand}{Data standardization. The default is stand = 0 and x is not transformed,
if stand = 1 each column of x has zero mean an if stand=2 also unit variance.}

\item{mth}{Method to estimate the number of factors and the common component (factors and loadings):
\itemize{
\item mth = 0 - the number of factors must be given by the user and
the model is estimated by Principal Components.
\item mth = 1 - the number of factors must be given by the user and
the model is estimated using Lam et al. (2011) methodology.
\item mth = 2 - the number of factors is estimated using Bai and Ng (2002) ICP1 criterion and
the model is estimated by Principal Components.
\item mth = 3 - the number of factors is estimated using Bai and Ng (2002) ICP1 criterion and
the model is estimated using Lam et al. (2011) methodology.
\item mth = 4 - the number of factors is estimated by applying once the Lam and Yao (2012) criterion and
the model is estimated using Lam et al. (2011) methodology (default method).
\item mth = 5 - the number of factors is estimated using Ahn and Horenstein (2013) test and
the model is estimated by Principal Components.
\item mth = 6 - the number of factors is estimated using  Caro and Peña (2020) test and
the model is estimated using Lam et al. (2011) methodology with the combined correlation matrix.
}}

\item{r}{Number of factors, default value is estimated by Lam and Yao (2012) criterion.}

\item{lagk}{Maximum number of lags considered in the combined matrix.
The default is lagk = 3.}
}
\value{
A list with the following items:
\itemize{
\item r - Estimated number of common factors, if mth=0, r is given by the user.
\item F -  Estimated common factor matrix (T x r).
\item L - Estimated loading  matrix (k x r).
\item E - Estimated noise matrix  (T x k).
\item VarF - Proportion of variability explained by the factor and the accumulated sum.
\item MarmaF - Matrix giving the number of AR, MA, seasonal AR and seasonal MA coefficients
for the Factors, plus the seasonal period and the number of non-seasonal and seasonal differences.
\item MarmaE - Matrix giving the number of AR, MA, seasonal AR and seasonal MA coefficients for the noises,
plus the seasonal period and the number of non-seasonal and seasonal differences.
}
}
\description{
The function estimates the Dynamic Factor Model by Principal Components and by the estimator of Lam et al. (2011).
}
\examples{
data(TaiwanAirBox032017)
dfm1 <- dfmpc(as.matrix(TaiwanAirBox032017[1:100,1:30]), mth=4)

}
\references{
Ahn, S. C. and Horenstein, A. R. (2013). Eigenvalue ratio test for the number
of factors. \emph{Econometrica}, 81(3):1203–1227.

Bai, J. and Ng, S. (2002). Determining the number of factors in approximate factor
models. \emph{Econometrica}, 70(1):191–221.

Caro, A. and Peña, D. (2020). A test for the number of factors in dynamic factor models.
UC3M Working papers. Statistics and Econometrics.

Lam, C. and Yao, Q. (2012). Factor modeling for high-dimensional time series:
inference for the number of factors. \emph{The Annals of Statistics}, 40(2):694–726.

Lam, C., Yao, Q., and Bathia, N. (2011). Estimation of latent factors for
high-dimensional time series. \emph{Biometrika}, 98(4):901–918.
}
