
.RunBootstrap <- function(df, cluster, outcome, post, treat) {
  # Bootstrap
  df_boot <- .SampleDf(df = df, cluster = cluster)

  # Format data
  df_format <- .orddid_extract_Y_cells(
    data = df_boot,
    outcome = outcome,
    post = post,
    treat = treat
  )

  # Estimation
  estimate <- .GetPointEstimate(df = df_format)$estimated_effects

  return(estimate)
}


#' Block bootstrap sampling: resample clusters with replacement
#' @noRd
#' @importFrom dplyr %>% distinct slice_sample inner_join across all_of mutate row_number select
.SampleDf <- function(df, cluster) {
  # Get unique cluster IDs
  cluster_ids <- df %>%
    distinct(across(all_of(cluster)))
  # Resample clusters with replacement
  sampled_clusters <- cluster_ids %>%
    slice_sample(n = nrow(cluster_ids), replace = TRUE) %>%
    mutate(.boot_id = row_number())
  # Collect all observations from sampled clusters
  df_boot <- sampled_clusters %>%
    inner_join(df, by = cluster, relationship = "many-to-many") %>%
    select(-".boot_id")
  return(df_boot)
}

#' @importFrom stats quantile sd
#' @noRd
.ComputeCI <- function(boot_estimates, alpha = 0.05) {
  # Collect bootstrap estimates for diff effects
  diff_effects <- do.call(rbind, boot_estimates$diff_effects)
  # Compute percentile CI
  lower_bound <- apply(diff_effects, 2, quantile, probs = alpha / 2)
  upper_bound <- apply(diff_effects, 2, quantile, probs = 1 - alpha / 2)

  # Collect estimates for relative effect
  relative_effects <- do.call(rbind, boot_estimates$relative_effect)
  # Compute percentile CI
  re_se <- apply(relative_effects, 2, sd)
  return(list(lower = lower_bound, upper = upper_bound, re_se = re_se))
}


#' Imbens–Manski (2004) CI for a partially identified scalar parameter
#'
#' @param lb_hat numeric, point estimate of the lower bound
#' @param ub_hat numeric, point estimate of the upper bound
#' @param se_lb  numeric, standard error of lb_hat (bootstrap or asymptotic)
#' @param se_ub  numeric, standard error of ub_hat
#' @param alpha  significance level (default 0.05)
#' @param truncate_diff logical; if TRUE, uses max(ub_hat - lb_hat, 0) in the
#'        critical-value equation (recommended for coverage)
#' @return list with c_crit, ci (length-2), delta, se_max
imbens_manski_ci <- function(lb_hat, ub_hat, se_lb, se_ub, alpha = 0.05,
                             truncate_diff = TRUE) {

  stopifnot(is.finite(lb_hat), is.finite(ub_hat),
            is.finite(se_lb),  is.finite(se_ub),
            se_lb >= 0, se_ub >= 0)

  if (se_lb == 0 && se_ub == 0) {
    warning("Both standard errors are zero; returning the estimated bounds.")
    return(list(c_crit = 0, ci = c(lb_hat, ub_hat),
                delta = NA_real_, se_max = 0))
  }

  se_max <- max(se_lb, se_ub)

  # signal-to-noise ratio in IM equation
  diff_hat <- ub_hat - lb_hat
  if (truncate_diff) diff_hat <- max(diff_hat, 0)
  delta <- diff_hat / se_max

  # Solve: Phi(c + delta) - Phi(-c) = 1 - alpha
  f <- function(c) stats::pnorm(c + delta) - stats::pnorm(-c) - (1 - alpha)

  # Bracket: for delta >= 0, c is in [z_{1-α}, z_{1-α/2}]
  lo <- stats::qnorm(1 - alpha)
  hi <- stats::qnorm(1 - alpha/2)
  f_lo <- f(lo); f_hi <- f(hi)

  if (is.na(f_lo) || is.na(f_hi) || f_lo * f_hi > 0) {
    # widen if needed (extreme inputs)
    root <- stats::uniroot(f, interval = c(0, 8))$root
  } else {
    root <- stats::uniroot(f, interval = c(lo, hi))$root
  }

  c_star <- as.numeric(root)

  ci <- c(lb_hat - c_star * se_lb,
          ub_hat + c_star * se_ub)

  list(
    bound = c(lb = lb_hat, ub = ub_hat),
    c_crit = c_star, ci = ci, delta = delta, se_max = se_max)
}
