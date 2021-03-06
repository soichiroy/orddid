% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/orddid-diagnostic.R
\name{equivalence_test}
\alias{equivalence_test}
\title{Conduct an equivalence test of the distributional parallel trends assumption.}
\usage{
equivalence_test(object, alpha = 0.05, threshold = NULL)
}
\arguments{
\item{object}{A fitted object from \code{\link{ord_did}}.}

\item{alpha}{The level of a test. This value should take between 0 and 1. Default is 0.05.}

\item{threshold}{An equivalance threshold. 
If left as \code{NULL}, the data-driven threshold, estimated in \code{\link{calc_threshold}}, is used for the test.}
}
\value{
\code{equivalence_test()} returns a list of class `orddid.test', which contains the following items:
   \item{tv}{A vector of point-wise deviation between q1(v) and q0(v).}
   \item{tv_var}{A vector of variances for each t(v).}
   \item{tmax}{A maximum deviation of q1(v) and q0(v).}
   \item{v_range}{A range of v used to evaluate q1 and q0.}
   \item{Uv}{Point wise \code{1 - alpha} level upper confidence interval.}
   \item{Lv}{Point wise \code{1 - alpha} level lower confidence interval.}
   \item{Umax}{Maximum upper bound.}
   \item{Lmin}{Minimum lower bound.}
   \item{Upvalue}{Point-wise pvalues associated with the upper bounds.}
   \item{Lpvalue}{Point-wise pvalues associated with the lower bounds.}
   \item{pvalue}{P-value of the test.}
   \item{zscore}{Z-score of the test.}
   \item{reject}{Decision of the equivalance test. If \code{TRUE}, the test rejects the null of non-equivalance.}
}
\description{
\code{equivalence_test()} implements an equivalance test to assess the distributional parallel trends assumption
 using the data from the pre-treatment periods.
}
