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

#' Log-likelihood Function
#'
#' @param par A vector of parameters.
#'  means: \code{par[1:3]}; scale: \code{par[4:6]}; cutoff: \code{par[7:]} (only when J > 3)
#' @param Y A vector of categorical outcome.
#' @param cut A set of two cutoffs fixed for identification. Default is \code{cut = c(0, 1)}.
#' @param id_group A vector of group indicator.
#' @keywords internal
log_like_probit_group <- function(par, Y, cut, id_group) {
  ## prep parameters
  n_group <- length(unique(id_group))
  mu    <- par[1:n_group]
  ss    <- exp(par[(n_group+1):(n_group*2)])
  j_min <- min(Y)
  j_max <- max(Y)

  ## if the number of categories is more than 4,
  ## we estimate cutoff
  if (length(par[-c(1:(n_group*2))]) > 0) {
    cut <- c(cut, cut[length(cut)] + cumsum(exp(par[-c(1:(n_group*2))])))
  }

  ## evalute the likelihood
  ll    <- 0
  for (g in 1:n_group) {
    ## get group params
    mu_g <- mu[g];
    ss_g <- ss[g]
    item  <- 1

    ## evaluate ll function 
    nj   <- sum(Y[id_group == g] == j_min)
    ll   <- ll + nj * log(pnorm((cut[item] - mu_g) / ss_g))
    item <- item + 1

    for (j in (j_min+1):(j_max-1)) {
      nj <- sum(Y[id_group == g] == j)
      ll <- ll + nj * log(
        pnorm((cut[item] - mu_g) / ss_g) -
        pnorm((cut[item-1] - mu_g) / ss_g)
      )
      item <- item + 1
    }

    nj <- sum(Y[id_group == g] == j_max)
    ll <- ll + nj * log(1 - pnorm((max(cut) - mu_g) / ss_g))
  }

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

  ## initial value setup
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

  ## fit ordered probit
  fit <- optim(
    par = par_init, fn = log_like_probit,
    Y = Y, cut = cut, method = 'BFGS'
  )

  ## return cutoff
  cutoff <- c(cut, cut[length(cut)] + cumsum(exp(fit$par[-c(1,2)])))
  return(list(mu = fit$par[1], sd = exp(fit$par[2]), ll = fit$value, cutoff = cutoff))
}


#' Probit regression with group dummies
#'
#' @param Y A vector of the ordinal outcome.
#' @param id_group A vector of group indicator.
#' @keywords internal
fit_ord_probit_gr <- function(Y, id_group, init = NULL, cut) {
  ## input check
  ord_probit_check_cat(Y)
  n_group <- length(unique(id_group))

  ## initial value setup
  if (is.null(init)) {
    par_init <- c(rep(0.1, n_group), rep(0.5, n_group))
    ## initialize cutoff
    n_cat <- (max(Y) - min(Y)) + 1
    if (n_cat > 3) {
      n_cutoffs <- n_cat - 1
      ## first two cutoffs are fixed: default is c(0, 1)
      par_init <- c(par_init, rep(-0.5, n_cutoffs-2))
    }
  }

  ## fit ordered probit
  fit <- optim(
    par = par_init, fn = log_like_probit_group,
    Y = Y, cut = cut, id_group = id_group,
    method = 'BFGS'
  )

  ## return object
  cutoff <- c(cut, cut[length(cut)] + cumsum(exp(fit$par[-c(1:(n_group*2))])))
  mu_vec <- fit$par[1:n_group]
  sd_vec <- fit$par[(n_group+1):(2*n_group)]

  return(list(mu = mu_vec, sd = exp(sd_vec), cutoff = cutoff, ll = fit$value))
}




#' Parameter Estimation for the Counterfactual Outcome
#'
#' @keywords internal
ord_did_run <- function(Ynew, Yold, treat, cut, pre = FALSE) {
  ## input check
  ord_did_check_input(Ynew, Yold, treat)
  dat <- create_group_dummy(Ynew, Yold, treat, pre = pre)

  ## fit probit regression

  fit_ct <- fit_ord_probit_gr(Y = dat$Y, id_group = dat$id_group, cut = cut)
  fit_tr <- fit_ord_probit(Y = Ynew[treat == 1], cut = cut)

  if (isTRUE(pre)) {
    ## fit_ct incldus "post-treatment x treatment" if pre = TRUE
    mu11 <- fit_ct$mu[4]
    ss11 <- fit_ct$sd[4]
  } else {
    ## estimated params
    mu00 <- fit_ct$mu[1]; mu01 <- fit_ct$mu[2]; mu10 <- fit_ct$mu[3]
    sd00 <- fit_ct$sd[1]; sd01 <- fit_ct$sd[2]; sd10 <- fit_ct$sd[3]

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

  # reject the bootstrap replica when optimization fails
  # typically rejection happens when outcome is not
  dat_tmp <- cbind(Ynew, Yold, treat)

  ## bootstrap -----------------------------------------------------------
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
