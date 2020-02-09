
#' Ordinal Difference-in-Differences for Panel Data
#'
#' A function to implement the difference-in-differences for the ordinal outcome.
#'
#' @param Ynew a vector of outcome for the post-treatment period.
#' @param Yold a vector of outcome for the pre-treatment period.
#' @param treat a vector of treatment indicator.
#' @param id_cluster a vector of cluster id.
#'   If left \code{NULL}, bootstrap is blocked at the individaul level
#' @param cut a vector of cutoffs.
#' @param n_boot the number of boostrapt iterations for computing uncertainty. Default is 500.
#' @param pre a boolean to indicate if data comes from entirely pre-treatment periods.
#' @param verbose If \code{TRUE}, print the current iteration in the bootstrap.
#' @return A list
#' @export
ord_did <- function(Ynew, Yold, treat, id_cluster = NULL, cut = c(0, 1),
                    n_boot, pre = FALSE, verbose = FALSE) {
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
