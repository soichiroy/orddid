
#' Ordinal Difference-in-Differences for Panel Data
#'
#' A function to implement Ordinal Difference-in-Differences method.
#'
#' @param Ynew a vector of outcome for the post-treatment period.
#' @param Yold a vector of outcome for the pre-treatment period.
#' @param treat a vector of treatment indicator.
#' @param cut a vector of cutoffs.
#' @param n_boot the number of boostrapt iterations for computing uncertainty. Default is 500.
#' @param pre a boolean to indicate if data comes from entirely pre-treatment periods.
#' @param verbose If \code{TRUE}, print the current iteration in the bootstrap.
#' @return A list
#' @export
ord_did <- function(Ynew, Yold, treat, cut, n_boot, pre = FALSE, verbose = FALSE) {
  ## fit on the obs data
  fit <- ord_did_run(Ynew, Yold, treat, cut, pre)

  ## bootstrap (assuming a panel)
  boot <- ord_did_boot(Ynew, Yold, treat, cut, n_boot, verbose)
  
  ## return objects
  return_list <- list(
    'fit' = fit, 
    'boot' = boot$boot_save, 
    'boot_params' = boot$boot_params)
  
  ## data summary 
  Yc <- c(Ynew, Yold)
  J  <- length(unique(Yc))
  attr(return_list, "input") <- list(
    n_boot = n_boot, pre = pre, cut = cut)
  attr(return_list, "n_choice") <- J 
  
  class(return_list) <- c("orddid", "orddid.fit")
  return(return_list)
}
