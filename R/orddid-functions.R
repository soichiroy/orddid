#
# title: functions called by ord_did() in orddid.R
#
# author: Soichiro Yamauchi
#
#

#' Log-likelihood Function
#' @keywords internal
log_like_probit <- function(par, Y, cut) {
  ## prep parameters
  mu    <- par[1]
  ss    <- exp(par[2])
  j_min <- min(Y)
  j_max <- max(Y)
  item  <- 1
  ll    <- 0

  ## if the number of categories is more than 4,
  ## we estimate cutoff
  if (length(par[-c(1,2)]) > 0) {
    cut <- c(cut, cut[length(cut)] + cumsum(exp(par[-c(1, 2)])))
  }

  ## evalute the likelihood
  nj   <- sum(Y == j_min)
  ll   <- nj * log(pnorm((cut[item] - mu)/ss))
  item <- item + 1

  for (j in (j_min+1):(j_max-1)) {
    nj <- sum(Y == j)
    ll <- ll + nj * log(
      pnorm((cut[item] - mu) / ss) -
      pnorm((cut[item-1] - mu) / ss)
    )
    item <- item + 1
  }

  nj <- sum(Y == j_max)
  ll <- ll + nj * log(1 - pnorm((max(cut) - mu) / ss))

  return(-ll)
}

#' Probit Regression
#'
#' Fitting ordered probit regression
#' @param Y A vector of outcome
#' @keywords internal
fit_ord_probit <- function(Y, init = NULL, cut) {
  ## input check
  ord_probit_check_cat(Y)

  ## fit ordered probit
  if (is.null(init)) {
    par_init <- c(0.1, 1)
    ## initialize cutoff
    n_cat <- (max(Y) - min(Y)) + 1
    if (n_cat > 3) {
      n_cutoffs <- n_cat - 1
      ## first two cutoffs are fixed: default is c(0, 1)
      par_init <- c(par_init, rep(-0.5, n_cutoffs-2))
    }
  }
  fit <- optim(
    par = par_init, fn = log_like_probit,
    Y = Y, cut = cut, method = 'BFGS')
  ## return cutoff
  cutoff <- c(cut, cut[length(cut)] + cumsum(exp(fit$par[-c(1,2)])))
  return(list(mu = fit$par[1], sd = exp(fit$par[2]), ll = fit$value, cutoff = cutoff))
}



#' Parameter Estimation
#'
#' Fitting probit model for each period
#' @keywords internal
ord_did_run <- function(Ynew, Yold, treat, cut, pre = FALSE) {

  ## fit four probit
  ## Y[d = 0, t = 0]
  fit00 <- fit_ord_probit( Y = Yold[treat == 0], cut = cut )
  ## Y[d = 0, t = 1]
  fit01 <- fit_ord_probit( Y = Ynew[treat == 0], cut = cut )
  ## Y[d = 1, t = 0]
  fit10 <- fit_ord_probit( Y = Yold[treat == 1], cut = cut )
  # observed Y[t=1,d=1]
  fit_tr <- fit_ord_probit( Y = Ynew[treat == 1], cut = cut )

  ## estimate counter factual distribution
  if (isTRUE(pre)) {
    ## this is pre-treatment so fit_tr has mu11 and sd11
    mu11 <- fit_tr$mu
    ss11 <- fit_tr$sd
  } else {
    ## use the identification formula to recover counter-factual parameters
    mu11 <- fit10$mu + (fit01$mu - fit00$mu) *
            (fit10$sd / fit00$sd)
    ss11 <- fit10$sd * fit01$sd / fit00$sd
  }

  ## compute the probability
  cut_new <- c(-Inf, cut, Inf)
  Ypred <- rep(NA, (length(cut)+1))
  for (j in 1:(length(cut)+1)) {
    Ypred[j] <- pnorm(cut_new[j+1], mean = mu11, sd = ss11)  -
                pnorm(cut_new[j], mean = mu11, sd = ss11)
  }

  ## compute observed Pr(Y(1) = j | D)
  Yobs <- as.vector(prop.table(table(Ynew[treat==1])))

  ## estimatd parameters
  theta_est <- list(
    mu00 = fit00$mu, sd00 = fit00$sd,
    mu01 = fit01$mu, sd01 = fit01$sd,
    mu10 = fit10$mu, sd10 = fit10$sd,
    mu11 = mu11, sd11 = ss11
  )

  return(list("Y1" = Yobs, "Y0" = Ypred, 'mu11' = mu11, 'ss11' = ss11,
    fit00 = fit00, fit01 = fit01, fit10 = fit10, fit_tr = fit_tr,
    theta = theta_est
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

  ## assuming panel structure
  ## this check is minimum
  if (length(Ynew) == length(Yold)) {
    # define an iterator
    b <- 1

    # reject the bootstrap replica when optimization fails
    # typically rejection happens when outcome is not
    dat_tmp <- cbind(Ynew, Yold, treat)

    while(b <= n_boot) {
      tryCatch({
        # sample bootstrap index
        dat_boot <- block_sample(dat_tmp, id_cluster)

        # fit the model
        fit_tmp <- ord_did_run(
          Ynew  = dat_boot[, 1],
          Yold  = dat_boot[ ,2],
          treat = dat_boot[, 3],
          cut   = cut
        )

        # save
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

      ## verbose
      if (isTRUE(verbose)) {
        if ((b %% iter_show) == 0) {
            cat('\r', b, "out of", n_boot, "bootstrap iterations")
            flush.console()
        }
      }

    }

    # clear the console
    cat("\n")
    # concatenate all params
    boot_params <- do.call("rbind", boot_params_save)
  } else {
    stop('length of two outcomes does not match')
  }

  return(list("boot_params" = boot_params, "boot_save" = boot_save))
}



#' Block re-sampling for block bootstrap
#'
#' @param dat a data frame with \code{Ynew}, \code{Yold} and \code{treat}.
#' @param id_cluster a cluster id vector of length n.
#' @return a list of resampled data.
#' @keywords internal
block_sample <- function(dat, id_cluster) {
  if (is.null(id_cluster)) {
    idx_use  <- sample(1:nrow(dat), size = nrow(dat), replace = TRUE)
    dat_boot <- dat[idx_use, ]
  } else {
    id_unique <- unique(id_cluster)   # unique cluster id
    J         <- length(id_unique)    # number of clusters

    # sample cluster id & data
    id_cluster_boot <- sample(id_unique, size = J, replace = TRUE)
    dat_boot <- dat_block_boot(dat = dat, id_cluster = id_cluster,
      id_cluster_boot = id_cluster_boot, max_cluster_size = max(table(id_cluster))
    )
  }

  return(dat_boot)
}
