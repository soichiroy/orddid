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
#' @param tau_assumption Character vector controlling how the relative-effect
#'   bounds are computed. \code{"none"} uses the unconstrained sharp bounds,
#'   while \code{"stochastic_dominance"} imposes positive regression dependence
#'   and computes the resulting bounds by linear programming. By default, both
#'   are returned.
#'
#' @return An S3 object of class \code{"orddid_staggered"} with:
#' \item{group_time_effects}{Data frame of category-specific effects by (g,t).}
#' \item{group_time_tau}{Data frame of relative treatment effect bounds by (g,t).}
#' \item{aggregate_effects}{Weighted average effects with CIs.}
#' \item{inputs}{Function arguments for reproducibility.}
#' \item{n_boot_valid}{Number of bootstrap draws used for inference.}
#' \item{n_boot_fail}{Number of rejected or failed bootstrap draws.}
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
  use_parallel = FALSE,
  tau_assumption = c("none", "stochastic_dominance")
) {
  tau_assumption <- .orddid_match_tau_assumptions(tau_assumption)

  if (is.null(cluster)) cluster <- unit

  n_categories <- .orddid_pooled_category_count(
    list(outcome = data[[outcome]])
  )

  ## Point estimate
  point <- .staggered_estimate_all(
    data = data,
    outcome = outcome,
    unit = unit,
    time = time,
    group = group,
    base_period = base_period,
    tau_assumption = tau_assumption,
    n_categories = n_categories
  )

  gt_res <- point$gt_results
  n_gt   <- length(gt_res)

  if (n_gt == 0) {
    stop("No valid post-treatment (group, time) pairs found.")
  }

  J <- length(gt_res[[1]]$zeta)
  gt_keys <- vapply(gt_res, function(r) {
    paste(r$group, r$time, sep = "::")
  }, character(1))
  tau_levels <- gt_res[[1]]$tau$tau_assumption
  n_tau <- length(tau_levels)

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
      group = r$group,
      time = r$time,
      tau_assumption = r$tau$tau_assumption,
      tau_lb = r$tau$effect_lb,
      tau_ub = r$tau$effect_ub,
      status = r$tau$status,
      stringsAsFactors = FALSE
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
  boot_gt_tau  <- array(NA, dim = c(n_boot, n_gt, n_tau, 2))
  boot_agg_zeta <- matrix(NA, n_boot, J)
  boot_agg_tau  <- array(NA, dim = c(n_boot, n_tau, 2))

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
        base_period = base_period,
        tau_assumption = tau_assumption,
        n_categories = n_categories
      ),
      error = function(e) NULL
    )

    boot_keys <- if (is.null(res_b)) {
      character()
    } else {
      vapply(res_b$gt_results, function(r) {
        paste(r$group, r$time, sep = "::")
      }, character(1))
    }
    if (is.null(res_b) || !identical(boot_keys, gt_keys)) {
      n_fail <- n_fail + 1
      next
    }

    for (i in seq_len(n_gt)) {
      boot_gt_zeta[b, i, ] <- res_b$gt_results[[i]]$zeta
      tau_tab <- res_b$gt_results[[i]]$tau
      tau_tab <- tau_tab[match(tau_levels, tau_tab$tau_assumption), , drop = FALSE]
      boot_gt_tau[b, i, , 1] <- tau_tab$effect_lb
      boot_gt_tau[b, i, , 2] <- tau_tab$effect_ub
    }

    agg_b <- .staggered_aggregate(res_b$gt_results, weights = weights)
    boot_agg_zeta[b, ] <- agg_b$agg_zeta
    agg_tau_tab <- agg_b$agg_tau
    agg_tau_tab <- agg_tau_tab[
      match(tau_levels, agg_tau_tab$tau_assumption),
      ,
      drop = FALSE
    ]
    boot_agg_tau[b, , 1] <- agg_tau_tab$tau_lb
    boot_agg_tau[b, , 2] <- agg_tau_tab$tau_ub
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
  tau_df$c_crit <- NA_real_
  tau_df$n_boot_valid <- 0L

  if (sum(valid) > 10) {
    for (i in seq_len(n_gt)) {
      tau_res <- .FinalizeRelativeEffects(
        relative_effects = gt_res[[i]]$tau,
        re_summary = .SummariseRelativeBootstrap(
          lapply(seq_len(n_boot), function(b) {
            tibble::tibble(
              tau_assumption = tau_levels,
              effect_lb = boot_gt_tau[b, i, , 1],
              effect_ub = boot_gt_tau[b, i, , 2],
              status = ifelse(
                is.finite(boot_gt_tau[b, i, , 1]) &
                  is.finite(boot_gt_tau[b, i, , 2]),
                "ok",
                "infeasible"
              )
            )
          })[valid],
          tau_assumptions = tau_levels
        ),
        alpha = alpha
      )
      row_idx <- tau_df$group == gt_res[[i]]$group & tau_df$time == gt_res[[i]]$time
      tau_df$ci_lower[row_idx] <- tau_res$lower_ci
      tau_df$ci_upper[row_idx] <- tau_res$upper_ci
      tau_df$c_crit[row_idx] <- tau_res$c_crit
      tau_df$n_boot_valid[row_idx] <- tau_res$n_boot_valid
    }
  }

  ## Aggregate CIs
  agg_zeta_ci <- matrix(NA, J, 2)

  if (sum(valid) > 10) {
    for (j in seq_len(J)) {
      agg_zeta_ci[j, ] <- stats::quantile(boot_agg_zeta[valid, j],
                                           c(alpha / 2, 1 - alpha / 2))
    }
  }

  agg_tau_point <- agg$agg_tau
  names(agg_tau_point)[names(agg_tau_point) == "tau_lb"] <- "effect_lb"
  names(agg_tau_point)[names(agg_tau_point) == "tau_ub"] <- "effect_ub"

  agg_tau_res <- .FinalizeRelativeEffects(
    relative_effects = agg_tau_point,
    re_summary = .SummariseRelativeBootstrap(
      lapply(seq_len(n_boot), function(b) {
        tibble::tibble(
          tau_assumption = tau_levels,
          effect_lb = boot_agg_tau[b, , 1],
          effect_ub = boot_agg_tau[b, , 2],
          status = ifelse(
            is.finite(boot_agg_tau[b, , 1]) & is.finite(boot_agg_tau[b, , 2]),
            "ok",
            "infeasible"
          )
        )
      })[valid],
      tau_assumptions = tau_levels
    ),
    alpha = alpha
  )
  names(agg_tau_res)[names(agg_tau_res) == "effect_lb"] <- "tau_lb"
  names(agg_tau_res)[names(agg_tau_res) == "effect_ub"] <- "tau_ub"

  agg_effects <- list(
    agg_zeta    = agg$agg_zeta,
    agg_zeta_ci = agg_zeta_ci,
    agg_tau     = agg_tau_res
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
      weights = weights,
      tau_assumption = tau_assumption
    ),
    n_boot_valid = sum(valid),
    n_boot_fail = n_fail
  )

  class(result) <- "orddid_staggered"
  result
}
