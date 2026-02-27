## ============================================================
## orddid-staggered-internals.R
## Three-step parametric estimator for staggered adoption
## ============================================================

#' Estimate all group-time effects for staggered design
#'
#' @param data long-format data frame
#' @param outcome column name for outcome
#' @param unit column name for unit id
#' @param time column name for time period
#' @param group column name for adoption group (Inf or NA = never-treated)
#' @param base_period "previous" or an integer
#' @return list with group_time results and fitted parameters
#' @keywords internal
.staggered_estimate_all <- function(data, outcome, unit, time, group,
                                    base_period = "previous") {

  Y     <- data[[outcome]]
  t_col <- data[[time]]
  g_col <- data[[group]]

  ## Identify never-treated and adoption groups
  g_col[is.na(g_col)] <- Inf
  periods        <- sort(unique(t_col))
  adopt_groups   <- sort(unique(g_col[is.finite(g_col)]))
  t0             <- min(periods)

  ## ---- Step 1: baseline from never-treated at t=0 ----
  idx_inf0   <- !is.finite(g_col) & (t_col == t0)
  Y_inf0     <- Y[idx_inf0]
  fit_inf0   <- .FitOrderedProbit(Y = Y_inf0, cutoff_val = NULL)
  kappa      <- fit_inf0$cutoff

  ## ---- Step 2: estimate (mu, sigma) for each (group, time) ----
  ## Store parameter estimates indexed by (group_label, period)
  params <- list()

  ## Never-treated, all periods
  for (tt in periods) {
    idx <- !is.finite(g_col) & (t_col == tt)
    Y_gt <- Y[idx]
    if (length(Y_gt) < 5) next
    if (tt == t0) {
      params[[paste0("Inf_", tt)]] <- list(
        mu = fit_inf0$mu, sd = fit_inf0$sd, group = Inf, time = tt
      )
    } else {
      fit_gt <- .FitOrderedProbit(Y = Y_gt, cutoff_val = kappa)
      params[[paste0("Inf_", tt)]] <- list(
        mu = fit_gt$mu, sd = fit_gt$sd, group = Inf, time = tt
      )
    }
  }

  ## Each adoption group, all periods
  for (gg in adopt_groups) {
    for (tt in periods) {
      idx <- (g_col == gg) & (t_col == tt)
      Y_gt <- Y[idx]
      if (length(Y_gt) < 5) next
      fit_gt <- .FitOrderedProbit(Y = Y_gt, cutoff_val = kappa)
      params[[paste0(gg, "_", tt)]] <- list(
        mu = fit_gt$mu, sd = fit_gt$sd, group = gg, time = tt
      )
    }
  }

  ## ---- Step 3: DID for post-treatment (g,t) pairs ----
  gt_results <- list()

  for (gg in adopt_groups) {
    post_periods <- periods[periods >= gg]

    for (tt in post_periods) {
      ## Choose base period s
      if (is.character(base_period) && base_period == "previous") {
        s <- max(periods[periods < gg])
      } else {
        s <- as.integer(base_period)
      }

      ## Retrieve parameters
      key_gs     <- paste0(gg, "_", s)
      key_inf_t  <- paste0("Inf_", tt)
      key_inf_s  <- paste0("Inf_", s)

      if (!all(c(key_gs, key_inf_t, key_inf_s) %in% names(params))) next

      p_gs    <- params[[key_gs]]
      p_inf_t <- params[[key_inf_t]]
      p_inf_s <- params[[key_inf_s]]

      ## DID formula for counterfactual:
      ## mu_gt(0) = mu_gs + sigma_gs * (mu_inf_t - mu_inf_s) / sigma_inf_s
      ## sigma_gt(0) = sigma_gs * sigma_inf_t / sigma_inf_s
      mu_cf    <- p_gs$mu + p_gs$sd * (p_inf_t$mu - p_inf_s$mu) / p_inf_s$sd
      sigma_cf <- p_gs$sd * p_inf_t$sd / p_inf_s$sd

      ## Counterfactual distribution
      prob_cf <- .EstimateOutcomeProportions(
        mu = mu_cf, sd = sigma_cf, cutoff = kappa
      )

      ## Observed distribution for (g, t)
      idx_obs <- (g_col == gg) & (t_col == tt)
      Y_obs   <- Y[idx_obs]
      J       <- length(kappa) + 1

      prob_obs <- numeric(J)
      cats     <- sort(unique(Y))
      for (j in seq_len(J)) {
        prob_obs[j] <- mean(Y_obs == cats[j])
      }

      ## Treatment effects
      zeta <- .ComputeTreatmentEffects(y1_prop = prob_obs, y0_prop = prob_cf)
      tau  <- .ComputeRelativeEffect(p1 = prob_obs, p0 = prob_cf)

      gt_results[[length(gt_results) + 1]] <- list(
        group     = gg,
        time      = tt,
        zeta      = zeta,
        tau       = tau,
        prob_obs  = prob_obs,
        prob_cf   = prob_cf,
        mu_cf     = mu_cf,
        sigma_cf  = sigma_cf
      )
    }
  }

  list(
    gt_results = gt_results,
    params     = params,
    kappa      = kappa,
    periods    = periods,
    adopt_groups = adopt_groups
  )
}


#' Aggregate group-time effects using weights
#'
#' @param gt_results list of group-time results from .staggered_estimate_all
#' @param weights "uniform" or named list
#' @return list with aggregate zeta and tau
#' @keywords internal
.staggered_aggregate <- function(gt_results, weights = "uniform") {
  n_gt <- length(gt_results)
  if (n_gt == 0) return(list(agg_zeta = NULL, agg_tau = NULL))

  J <- length(gt_results[[1]]$zeta)

  if (is.character(weights) && weights == "uniform") {
    w <- rep(1 / n_gt, n_gt)
  } else {
    w <- weights
    stopifnot(length(w) == n_gt, abs(sum(w) - 1) < 1e-8)
  }

  agg_zeta <- numeric(J)
  agg_tau  <- c(lower = 0, upper = 0)

  for (i in seq_len(n_gt)) {
    agg_zeta <- agg_zeta + w[i] * gt_results[[i]]$zeta
    agg_tau  <- agg_tau  + w[i] * gt_results[[i]]$tau
  }

  list(agg_zeta = agg_zeta, agg_tau = agg_tau)
}
