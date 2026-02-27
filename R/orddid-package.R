#' Difference-in-Differences for Ordinal Outcomes
#'
#' This package provides tools for estimating treatment effects
#'  and diagnosing assumptions with ordinal outcomes in the difference-in-differences design.
#'
#' The package has the following main functions:
#' \itemize{
#'   \item \code{\link{ord_did}}: Implements the ordinal DID estimator.
#'   \item \code{\link{summary.orddid}}: Output the treatment effects and their uncertainties.
#'   \item \code{\link{equivalence_test}}: Test the assumption in the pre-treatment periods.
#'   \item \code{\link{plot.orddid.test}}: Visualize the result of the equivalance test.
#' }
#'
#' In addition to the above functions, \code{orddid} packages contains two example datasets:
#'  \code{\link{gun_twowave}} and \code{\link{gun_threewave}}.
"_PACKAGE"

# Suppress R CMD check notes for foreach iterator variables (NSE)
utils::globalVariables(c("i", "b"))
