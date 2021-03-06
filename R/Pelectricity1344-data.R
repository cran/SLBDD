#' Electricity Prices in New England and USA
#'
#' Weakly series of electricity price each hour of each day during 678 weeks
#' in the 8 regions in New England. We have 1344 series corresponding to each
#' of the seven days, one of the 24 hours and for one of the regions \eqn{(7x24x8=1344)}.
#' The first series correspond to the price in the first region at 1 am CT of Thursday
#' 01/01/2004, the second to 2 am, same day and so on. Thus the first 192 series
#' (24 hours x 8 regions) are the price of all the hours of Thursday in the eight regions,
#' the next 192 are for Friday and so on. These series were used in the articles Alonso and Peña (2019),
#' and Peña, Tsay and Zamar (2019).
#' The series have been corrected of missing values at days of changing time in saving energy days.
#'
#' @docType data
#'
#' @usage data(PElectricity1344)
#'
#' @format An object of class \code{"data.frame"}.
#'
#' @keywords datasets
#'
#' @references Alonso, A. M. and Peña, D. (2019). Clustering time series by linear dependency.
#' \emph{Statistics and Computing}, 29(4):655–676.
#'
#' Peña, D. Tsay, R. and Zamar, R. (2019). Empirical Dynamic Quantiles for Visualization of High-
#' Dimensional Time Series, \emph{Technometrics}, 61:4, 429-444.
#'
"PElectricity1344"
