% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/orddid-checks.R
\name{ord_probit_check_cat}
\alias{ord_probit_check_cat}
\title{Input check for \code{fit_ord_probit}.}
\usage{
ord_probit_check_cat(Y)
}
\arguments{
\item{Y}{a vector of ategorical outcome.}
}
\description{
This function check if there is missing categories. For example, if Ymin = 0 and Ymax=J-1.
 the function checks there are J unique values in Y.
}
\keyword{internal}
