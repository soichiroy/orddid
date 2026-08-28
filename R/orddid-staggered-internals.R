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
#' @param tau_assumption relative-effect bound assumption passed to
#'   \code{.ComputeRelativeEffect()}
#' @return list with group_time results and fitted parameters
#' @keywords internal
.staggered_estimate_all <- function(data, outcome, unit, time, group,
                                    base_period = "previous",
                                    tau_assumption = "none",
                                    n_categories = NULL) {
  tau_assumption <- .orddid_match_tau_assumptions(tau_assumption)

  Y     <- data[[outcome]]
  t_col <- data[[time]]
  g_col <- data[[group]]

  if (is.null(n_categories)) {
    n_categories <- .orddid_pooled_category_count(list(outcome = Y))
  }

  ## Identify never-treated and adoption groups
  g_col[is.na(g_col)] <- Inf
  periods        <- sort(unique(t_col))
  adopt_groups   <- sort(unique(g_col[is.finite(g_col)]))
  t0             <- min(periods)
  category_levels <- sort(unique(Y[!is.na(Y)]))

  base_for_group <- function(gg) {
    if (is.character(base_period) && base_period == "previous") {
      candidates <- periods[periods < gg]
      if (length(candidates) == 0L) return(NA_real_)
      max(candidates)
    } else {
      as.numeric(base_period)
    }
  }
  group_bases <- vapply(adopt_groups, base_for_group, numeric(1))
  required_never_periods <- sort(unique(c(
    t0,
    group_bases[is.finite(group_bases)],
    unlist(lapply(adopt_groups, function(gg) periods[periods >= gg]))
  )))

  ## ---- Step 1: baseline from never-treated at t=0 ----
  idx_inf0   <- !is.finite(g_col) & (t_col == t0)
  Y_inf0     <- Y[idx_inf0]
  .orddid_require_full_category_support(
    cells = list(never_treated_baseline = Y_inf0),
    n_categories = n_categories,
    context = "The staggered never-treated baseline cell"
  )
  fit_inf0   <- .FitOrderedProbit(Y = Y_inf0, cutoff_val = NULL)
  kappa      <- fit_inf0$cutoff

  ## ---- Step 2: estimate (mu, sigma) for each (group, time) ----
  ## Store parameter estimates indexed by (group_label, period)
  params <- list()

  ## Never-treated, all periods
  for (tt in required_never_periods) {
    idx <- !is.finite(g_col) & (t_col == tt)
    Y_gt <- Y[idx]
    if (length(Y_gt) < 5) next
    .orddid_require_full_category_support(
      cells = list(never_treated = Y_gt),
      n_categories = n_categories,
      context = paste0("The staggered never-treated cell at time ", tt)
    )
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

  ## Each adoption group only needs its pre-treatment base-period fit. The
  ## treated post-period distribution is tabulated directly below against the
  ## fixed global category levels and may validly contain zero-count categories.
  for (i in seq_along(adopt_groups)) {
    gg <- adopt_groups[i]
    s <- group_bases[i]
    if (!is.finite(s)) next
    idx <- (g_col == gg) & (t_col == s)
    Y_gs <- Y[idx]
    if (length(Y_gs) < 5) next
    .orddid_require_full_category_support(
      cells = list(adoption_group_base = Y_gs),
      n_categories = n_categories,
      context = paste0("The staggered group ", gg, " base-period cell at time ", s)
    )
    fit_gs <- .FitOrderedProbit(Y = Y_gs, cutoff_val = kappa)
    params[[paste0(gg, "_", s)]] <- list(
      mu = fit_gs$mu, sd = fit_gs$sd, group = gg, time = s
    )
  }

  ## ---- Step 3: DID for post-treatment (g,t) pairs ----
  gt_results <- list()

  for (gg in adopt_groups) {
    post_periods <- periods[periods >= gg]

    for (tt in post_periods) {
      ## Choose base period s
      s <- base_for_group(gg)
      if (!is.finite(s)) next

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
      for (j in seq_len(J)) {
        prob_obs[j] <- mean(Y_obs == category_levels[j])
      }

      ## Treatment effects
      zeta <- .ComputeTreatmentEffects(y1_prop = prob_obs, y0_prop = prob_cf)
      tau  <- .ComputeRelativeEffect(
        p1 = prob_obs,
        p0 = prob_cf,
        tau_assumption = tau_assumption
      )

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
  tau_assumptions <- gt_results[[1]]$tau$tau_assumption

  if (is.character(weights) && weights == "uniform") {
    w <- rep(1 / n_gt, n_gt)
  } else {
    w <- weights
    stopifnot(length(w) == n_gt, abs(sum(w) - 1) < 1e-8)
  }

  agg_zeta <- numeric(J)
  for (i in seq_len(n_gt)) {
    agg_zeta <- agg_zeta + w[i] * gt_results[[i]]$zeta
  }

  agg_tau <- lapply(tau_assumptions, function(assump) {
    tau_rows <- lapply(gt_results, function(res) {
      res$tau[res$tau$tau_assumption == assump, , drop = FALSE]
    })

    tau_ok <- vapply(tau_rows, function(x) {
      nrow(x) == 1L && identical(x$status[1], "ok")
    }, logical(1))

    tibble::tibble(
      tau_assumption = assump,
      tau_lb = if (all(tau_ok)) {
        sum(w * vapply(tau_rows, function(x) x$effect_lb[1], numeric(1)))
      } else {
        NA_real_
      },
      tau_ub = if (all(tau_ok)) {
        sum(w * vapply(tau_rows, function(x) x$effect_ub[1], numeric(1)))
      } else {
        NA_real_
      },
      status = if (all(tau_ok)) "ok" else "infeasible"
    )
  })

  agg_tau <- dplyr::bind_rows(agg_tau)

  list(agg_zeta = agg_zeta, agg_tau = agg_tau)
}
