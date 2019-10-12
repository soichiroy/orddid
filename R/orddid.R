
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
#' @return A list
#' @export
ord_did <- function(Ynew, Yold, treat, cut, n_boot, pre) {
  ## fit on the obs data
  fit <- ord_did_run(Ynew, Yold, treat, cut, pre)

  ## bootstrap (assuming a panel)
  boot <- ord_did_boot(Ynew, Yold, treat, cut, n_boot)
  
  ## return objects
  return(list('fit' = fit, 
    'boot' = boot$boot_save, 
    'boot_params' = boot$boot_params))
}
