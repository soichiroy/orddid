
#' Ordinal Difference-in-Differences for Panel Data
#'
#' \code{ord_did()} implements the difference-in-differences for the ordinal outcome.
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
#' @return \code{ord_did()} returns a list of class `orddid' containing the following components:
#' \item{fit}{A list with the output of the ordinal DID estimators,
#'            which contains parameter estimates and predicted probabilities for each category.}
#' \item{boot}{A list with the output of bootstraps,
#'             which contains parameter estimates and predicted probabilities for each category.}
#' \item{boot_params}{A list with all objects generated during the bootstrap step.}
#' @examples
#'\donttest{
#' ## load packages
#' library(orddid)
#' library(dplyr)
#'
#' ## load example data
#' data("gun_twowave")
#'
#' ## run
#' ## fit the ordinal DID
#' set.seed(1234)
#' fit <- ord_did(
#'   Ynew = gun_twowave %>% filter(year == 2012) %>% pull(guns),
#'   Yold = gun_twowave %>% filter(year == 2010) %>% pull(guns),
#'   treat = gun_twowave %>% filter(year == 2012) %>% pull(treat_100mi),
#'   id_cluster = gun_twowave %>% filter(year == 2010) %>% pull(reszip),
#'   n_boot = 10,
#'   pre = FALSE,
#'   verbose = FALSE
#' )
#'
#' ## view summary of the output
#' ## non-cumulative effects
#' summary(fit, cumulative = FALSE)
#'
#' ## cumulative effects
#' summary(fit)
#' }
#' @export
ord_did <- function(Ynew, Yold, treat, id_cluster = NULL, cut = c(0, 1),
                    n_boot = 500, pre = FALSE, verbose = FALSE) {

  ## quick input check
  ord_did_check_input(Ynew, Yold, treat, id_cluster, n_boot)

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
