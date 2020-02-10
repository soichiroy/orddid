#
# title: functions called by ord_did() in orddid.R
#
# author: Soichiro Yamauchi
#
#

#' Parameter Estimation for the Counterfactual Outcome
#'
#' @keywords internal
ord_did_run <- function(Ynew, Yold, treat, cut, pre = FALSE) {
  ## input check
  ord_did_run_check_input(Ynew, Yold, treat)
  dat <- create_group_dummy(Ynew, Yold, treat, pre = pre)

  ## fit probit regression
  fit_ct <- fit_ord_probit_gr(Y = dat$Y, id_group = dat$id_group, cut = cut)
  fit_tr <- fit_ord_probit(Y = Ynew[treat == 1], cut = cut)

  ## estimated params
  mu00 <- fit_ct$mu[1]; mu01 <- fit_ct$mu[2]; mu10 <- fit_ct$mu[3]
  sd00 <- fit_ct$sd[1]; sd01 <- fit_ct$sd[2]; sd10 <- fit_ct$sd[3]

  if (isTRUE(pre)) {
    ## fit_ct incldus "post-treatment x treatment" if pre = TRUE
    mu11 <- fit_ct$mu[4]
    ss11 <- fit_ct$sd[4]
  } else {
    ## identify mean and variance of Y11
    mu11 <- mu10 + (mu01 - mu00) * (sd10 / sd00)
    ss11 <- sd10 * sd01 / sd00
  }

  ##
  ## compute the probability Pr(Y(0) = j | D = 1)
  ##
  cuotff <- fit_ct$cutoff
  cut_new <- c(-Inf, cuotff, Inf)
  Ypred <- rep(NA, (length(cuotff)+1))
  for (j in 1:(length(cuotff)+1)) {
    Ypred[j] <- pnorm(cut_new[j+1], mean = mu11, sd = ss11)  -
                pnorm(cut_new[j], mean = mu11, sd = ss11)
  }

  ##
  ## compute observed Pr(Y(1) = j | D = 1)
  ##
  Yobs <- as.vector(prop.table(table(Ynew[treat==1])))

  ##
  ## estimatd parameters
  ##
  theta_est <- list(
    mu00 = mu00, sd00 = sd00,
    mu01 = mu01, sd01 = sd01,
    mu10 = mu10, sd10 = sd10,
    mu11 = mu11, sd11 = ss11
  )

  return(list("Y1" = Yobs, "Y0" = Ypred, 'mu11' = mu11, 'ss11' = ss11,
    fit_ct = fit_ct, fit_tr = fit_tr, theta = theta_est
  ))
}


#' Bootstrap function
#'
#' Conduct a block bootstrap at the unit level to compute variances and CIs
#' @keywords internal
ord_did_boot <- function(Ynew, Yold, treat, cut, id_cluster, n_boot, verbose) {

  ## create an object to save bootparams
  boot_save <- list()
  boot_save[[1]] <- boot_save[[2]] <- matrix(NA, nrow = n_boot, ncol = length(cut)+1)
  boot_save[[3]] <- boot_save[[4]] <- rep(NA, n_boot)
  boot_params_save <- list()

  ## check verbose
  if (isTRUE(verbose)) {
    iter_show <- round(n_boot * 0.1)
  }

  ## convert id_cluster to numeric
  if (!is.null(id_cluster)) {
    id_cluster <- as.numeric(as.factor(id_cluster))
  }


  ##
  ## setup for bootstrap
  ##
  # define an iterator
  b <- 1
  # prep data
  dat_tmp <- cbind(Ynew, Yold, treat)
  if (!is.null(id_cluster)) {
    id_unique <- unique(id_cluster)   # unique cluster id
    J         <- length(id_unique)    # number of clusters
    max_cluster_size <- max(table(id_cluster))
  }

  ## bootstrap -----------------------------------------------------------
  # reject the bootstrap replica when optimization fails
  # typically rejection happens when outcome is not well balanced
  while(b <= n_boot) {
    tryCatch({
      # sample bootstrap index
      if (is.null(id_cluster)) {
        dat_boot <- block_unit_sample(dat_tmp)
      } else {
        dat_boot <- block_sample(dat_tmp, id_cluster, id_unique, J, max_cluster_size)
      }


      # fit the model
      fit_tmp <- ord_did_run(
        Ynew  = dat_boot[, 1],
        Yold  = dat_boot[ ,2],
        treat = dat_boot[, 3],
        cut   = cut
      )

      # save estimates
      boot_save[[1]][b,]    <- fit_tmp$Y1  # Yobs
      boot_save[[2]][b,]    <- fit_tmp$Y0  # Y(0)
      boot_save[[3]][b]     <- fit_tmp$mu11
      boot_save[[4]][b]     <- fit_tmp$ss11
      boot_params_save[[b]] <- unlist(fit_tmp$theta)

      # update iterator
      b <- b + 1
    }, error = function(e) {
      NULL
    })

    ## verbose ----------------------------------------------------
    if (isTRUE(verbose)) {
      if ((b %% iter_show) == 0) {
          cat('\r', b, "out of", n_boot, "bootstrap iterations")
          flush.console()
      }
    }
    ## end of verbose ---------------------------------------------

  }
  ## end of bootstrap iterations -----------------------------------------

  # clear the console
  cat("\n")

  # concatenate all params
  boot_params <- do.call("rbind", boot_params_save)

  return(list("boot_params" = boot_params, "boot_save" = boot_save))
}



#' Block re-sampling for block bootstrap
#'
#' @param dat a data frame with \code{Ynew}, \code{Yold} and \code{treat}.
#' @param id_cluster a cluster id vector of length n.
#' @return a list of resampled data.
#' @keywords internal
block_sample <- function(dat, id_cluster, id_unique, J, max_cluster_size) {
  # sample cluster id & data
  id_cluster_boot <- sample(id_unique, size = J, replace = TRUE)
  dat_boot <- dat_block_boot(dat = dat, id_cluster = id_cluster,
    id_cluster_boot = id_cluster_boot, max_cluster_size = max_cluster_size)
  return(dat_boot)
}


#' Unit level block bootstrap
#' @keywords internal
block_unit_sample <- function(dat) {
  idx_use  <- sample(1:nrow(dat), size = nrow(dat), replace = TRUE)
  dat_boot <- dat[idx_use, ]
  return(dat_boot)
}
