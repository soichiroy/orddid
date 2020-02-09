
#' Ordinal Difference-in-Differences for Panel Data
#'
#' A function to implement the difference-in-differences for the ordinal outcome.
#'
#' @param Ynew A numeric vector of ordinal outcome for the post-treatment period.
#' @param Yold A numeric vector of ordinal outcome for the pre-treatment period.
#' @param treat A numeric vector of treatment indicator. 
#'  The treatment group should take 1 and the control group should take 0.
#' @param id_cluster A vector of cluster id.
#'   If left as \code{NULL}, bootstrap is implemented at the individual level.
#' @param cut A vector of cutoffs. Two numeric values should be specified. Default is \code{cut = c(0, 1)}.
#' @param n_boot The number of boostrapt iterations for estimating the variance. Default is \code{n_boot = 500}.
#' @param pre A boolean argument used to indicate if the data comes entirely from pre-treatment periods.
#'  This should be \code{TRUE} when the output is supplied to \code{\link{equivalence_test}}.
#' @param verbose If \code{TRUE}, print the progress of bootstrap iterations.
#' @return \code{ord_did} returns a list of class `orddid' containing the following components:
#' \item{fit}{A list with the output of the ordinal DID estimators, 
#'            which contains parameter estimates and predicted probabilities for each category.}
#' \item{boot}{A list with the output of bootstraps,
#'             which contains parameter estimates and predicted probabilities for each category.}
#' \item{boot_params}{A list with all objects generated during the bootstrap step.}
#' @export
ord_did <- function(Ynew, Yold, treat, id_cluster = NULL, cut = c(0, 1),
                    n_boot = 500, pre = FALSE, verbose = FALSE) {
  ## fit on the obs data
  fit <- ord_did_run(Ynew, Yold, treat, cut, pre)

  ## bootstrap (assuming a panel)
  boot <- ord_did_boot(Ynew, Yold, treat, cut, id_cluster, n_boot, verbose)

  ## return objects
  return_list <- list(
    'fit' = fit,
    'boot' = boot$boot_save,
    'boot_params' = boot$boot_params)

  ## data summary
  Yc <- c(Ynew, Yold)
  J  <- length(unique(Yc))
  n  <- length(treat)
  n1 <- sum(treat)

  ## add attributes to the returning object
  attr(return_list, "input")    <- list(n_boot = n_boot, pre = pre, cut = cut)
  attr(return_list, "n_choice") <- J
  attr(return_list, "n")        <- n
  attr(return_list, "n1")       <- n1

  class(return_list) <- c("orddid", "orddid.fit")
  return(return_list)
}
