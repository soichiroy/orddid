
#' Ordinal Difference-in-Differences
#'
#' A function to implement Ordinal Difference-in-Differences method.
#'
#' @param Ynew a vector of outcome for the post-treatment period.
#' @param Yold a vector of outcome for the pre-treatment period.
#' @param treat a vector of treatment indicator.
#' @param cut a vector of cutoffs.
#' @param n_boot the number of boostrapt iterations for computing uncertainty. Default is 500.
#' @return A list
#' @export
ord_did <- function(Ynew, Yold, treat, cut, n_boot) {
  ## fit on the obs data
  fit <- ord_did_run(Ynew, Yold, treat, cut)

  ## bootstrap (assuming a panel)
  boot_save <- list()
  boot_save[[1]] <- boot_save[[2]] <- matrix(NA, nrow = n_boot, ncol = length(cut)+1)
  boot_save[[3]] <- boot_save[[4]] <- rep(NA, n_boot)
  if (length(Ynew) == length(Yold)) {
    b <- 1
    while(b <= n_boot) {
      dat_tmp <- cbind(Ynew, Yold, treat)
      idx_use <- sample(1:nrow(dat_tmp), size = nrow(dat_tmp), replace = TRUE)
      tryCatch({
        fit_tmp <- ord_did_run(Ynew = dat_tmp[idx_use,1],
          Yold = dat_tmp[idx_use,2], treat = dat_tmp[idx_use, 3], cut = cut)

        boot_save[[1]][b,] <- fit_tmp$Y1  # Yobs
        boot_save[[2]][b,] <- fit_tmp$Y0  # Y(0)
        boot_save[[3]][b] <- fit_tmp$mu11
        boot_save[[4]][b] <- fit_tmp$ss11
        b <- b + 1
      }, error = function(cond) {
          return(NULL)
      })
    }
  } else {
    stop('length of two outcomes does not match')
  }

  return(list('fit' = fit, 'boot' = boot_save))
}
