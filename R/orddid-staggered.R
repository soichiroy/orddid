## ============================================================
## orddid-staggered.R
## User-facing function for staggered adoption design
## ============================================================

#' Ordinal DID for Staggered Adoption Design
#'
#' Estimates group-time treatment effects for ordinal outcomes
#' with staggered treatment adoption using parametric (ordered probit)
#' identification.
#'
#' @param data Data frame in long format (one row per unit-period).
#' @param outcome Character; column name for the ordinal outcome.
#' @param unit Character; column name for the unit id.
#' @param time Character; column name for the time period.
#' @param group Character; column name for the adoption time.
#'   Use \code{Inf} or \code{NA} for never-treated units.
#' @param cluster Character; column name for bootstrap clustering.
#'   If \code{NULL}, clusters at the unit level.
#' @param n_boot Number of bootstrap iterations (default 500).
#' @param base_period \code{"previous"} (default) uses g-1 as base period;
#'   an integer uses that fixed period for all groups.
#' @param weights \code{"uniform"} (default) for equal weights on all (g,t) pairs.
#' @param use_parallel If \code{TRUE}, use parallel computing.
#'
#' @return An S3 object of class \code{"orddid_staggered"} with:
#' \item{group_time_effects}{Data frame of category-specific effects by (g,t).}
#' \item{group_time_tau}{Data frame of relative treatment effect bounds by (g,t).}
#' \item{aggregate_effects}{Weighted average effects with CIs.}
#' \item{inputs}{Function arguments for reproducibility.}
#'
#' @export
ord_did_staggered <- function(
  data,
  outcome,
  unit,
  time,
  group,
  cluster = NULL,
  n_boot = 500,
  base_period = "previous",
  weights = "uniform",
  use_parallel = FALSE
) {

  if (is.null(cluster)) cluster <- unit

  ## Point estimate
  point <- .staggered_estimate_all(
    data = data,
    outcome = outcome,
    unit = unit,
    time = time,
    group = group,
    base_period = base_period
  )

  gt_res <- point$gt_results
  n_gt   <- length(gt_res)

  if (n_gt == 0) {
    stop("No valid post-treatment (group, time) pairs found.")
  }

  J <- length(gt_res[[1]]$zeta)

  ## Build group_time_effects data frame
  gt_df <- do.call(rbind, lapply(gt_res, function(r) {
    data.frame(
      group    = r$group,
      time     = r$time,
      category = seq_len(J),
      effect   = r$zeta
    )
  }))

  ## Build group_time_tau data frame
  tau_df <- do.call(rbind, lapply(gt_res, function(r) {
    data.frame(
      group  = r$group,
      time   = r$time,
      tau_lb = r$tau[1],
      tau_ub = r$tau[2]
    )
  }))

  ## Aggregate
  agg <- .staggered_aggregate(gt_res, weights = weights)

  ## ---- Bootstrap ----
  cluster_col <- data[[cluster]]
  unique_ids  <- unique(cluster_col)
  n_units     <- length(unique_ids)
  n_fail      <- 0

  ## Storage for bootstrap draws
  boot_gt_zeta <- array(NA, dim = c(n_boot, n_gt, J))
  boot_gt_tau  <- array(NA, dim = c(n_boot, n_gt, 2))
  boot_agg_zeta <- matrix(NA, n_boot, J)
  boot_agg_tau  <- matrix(NA, n_boot, 2)

  for (b in seq_len(n_boot)) {
    ids_b <- sample(unique_ids, n_units, replace = TRUE)
    df_b  <- do.call(rbind, lapply(seq_along(ids_b), function(k) {
      data[cluster_col == ids_b[k], , drop = FALSE]
    }))

    res_b <- tryCatch(
      .staggered_estimate_all(
        data = df_b,
        outcome = outcome,
        unit = unit,
        time = time,
        group = group,
        base_period = base_period
      ),
      error = function(e) NULL
    )

    if (is.null(res_b) || length(res_b$gt_results) != n_gt) {
      n_fail <- n_fail + 1
      next
    }

    for (i in seq_len(n_gt)) {
      boot_gt_zeta[b, i, ] <- res_b$gt_results[[i]]$zeta
      boot_gt_tau[b, i, ]  <- res_b$gt_results[[i]]$tau
    }

    agg_b <- .staggered_aggregate(res_b$gt_results, weights = weights)
    boot_agg_zeta[b, ] <- agg_b$agg_zeta
    boot_agg_tau[b, ]  <- agg_b$agg_tau
  }

  ## Compute CIs for group-time effects
  valid <- !is.na(boot_gt_zeta[, 1, 1])
  alpha <- 0.05

  gt_df$se       <- NA
  gt_df$lower_ci <- NA
  gt_df$upper_ci <- NA

  if (sum(valid) > 10) {
    for (i in seq_len(n_gt)) {
      for (j in seq_len(J)) {
        row_idx <- (i - 1) * J + j
        vals <- boot_gt_zeta[valid, i, j]
        gt_df$se[row_idx]       <- stats::sd(vals)
        gt_df$lower_ci[row_idx] <- stats::quantile(vals, alpha / 2)
        gt_df$upper_ci[row_idx] <- stats::quantile(vals, 1 - alpha / 2)
      }
    }
  }

  ## CIs for group-time tau
  tau_df$ci_lower <- NA
  tau_df$ci_upper <- NA

  if (sum(valid) > 10) {
    for (i in seq_len(n_gt)) {
      tau_lb_vals <- boot_gt_tau[valid, i, 1]
      tau_ub_vals <- boot_gt_tau[valid, i, 2]
      se_lb <- stats::sd(tau_lb_vals)
      se_ub <- stats::sd(tau_ub_vals)
      tau_ci <- tryCatch(
        imbens_manski_ci(
          lb_hat = gt_res[[i]]$tau[1],
          ub_hat = gt_res[[i]]$tau[2],
          se_lb = se_lb,
          se_ub = se_ub,
          alpha = 0.05,
          truncate_diff = TRUE
        ),
        error = function(e) list(ci = c(NA, NA))
      )
      tau_df$ci_lower[i] <- tau_ci$ci[1]
      tau_df$ci_upper[i] <- tau_ci$ci[2]
    }
  }

  ## Aggregate CIs
  agg_zeta_ci <- matrix(NA, J, 2)
  agg_tau_ci  <- c(NA, NA)

  if (sum(valid) > 10) {
    for (j in seq_len(J)) {
      agg_zeta_ci[j, ] <- stats::quantile(boot_agg_zeta[valid, j],
                                           c(alpha / 2, 1 - alpha / 2))
    }
    agg_re_se <- apply(boot_agg_tau[valid, , drop = FALSE], 2, stats::sd)
    agg_tau_im <- tryCatch(
      imbens_manski_ci(
        lb_hat = agg$agg_tau[1],
        ub_hat = agg$agg_tau[2],
        se_lb = agg_re_se[1],
        se_ub = agg_re_se[2],
        alpha = 0.05,
        truncate_diff = TRUE
      ),
      error = function(e) list(ci = c(NA, NA))
    )
    agg_tau_ci <- agg_tau_im$ci
  }

  agg_effects <- list(
    agg_zeta    = agg$agg_zeta,
    agg_zeta_ci = agg_zeta_ci,
    agg_tau     = agg$agg_tau,
    agg_tau_ci  = agg_tau_ci
  )

  result <- list(
    group_time_effects = gt_df,
    group_time_tau     = tau_df,
    aggregate_effects  = agg_effects,
    inputs = list(
      outcome = outcome,
      unit = unit,
      time = time,
      group = group,
      cluster = cluster,
      n_boot = n_boot,
      base_period = base_period,
      weights = weights
    ),
    n_boot_fail = n_fail
  )

  class(result) <- "orddid_staggered"
  result
}
