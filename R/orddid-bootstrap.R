
.RunBootstrap <- function(df, cluster, outcome, post, treat,
                          tau_assumption = "none") {
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
  estimate <- .GetPointEstimate(
    df = df_format,
    tau_assumption = tau_assumption
  )$estimated_effects

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
  valid <- !vapply(boot_estimates$diff_effects, is.null, logical(1)) &
           !vapply(boot_estimates$relative_effect, is.null, logical(1))
  n_fail <- sum(!valid)

  if (!any(valid)) {
    stop("All bootstrap iterations failed.", call. = FALSE)
  }

  # Collect bootstrap estimates for diff effects
  diff_effects <- do.call(rbind, boot_estimates$diff_effects[valid])
  # Compute percentile CI
  lower_bound <- apply(diff_effects, 2, quantile, probs = alpha / 2)
  upper_bound <- apply(diff_effects, 2, quantile, probs = 1 - alpha / 2)

  return(list(
    lower = lower_bound,
    upper = upper_bound,
    re_summary = .SummariseRelativeBootstrap(
      boot_estimates$relative_effect[valid]
    ),
    n_fail = n_fail,
    n_valid = sum(valid)
  ))
}


.SummariseRelativeBootstrap <- function(relative_effect_list,
                                        tau_assumptions = NULL) {
  if (length(relative_effect_list) == 0) {
    return(tibble::tibble(
      tau_assumption = character(),
      se_lb = numeric(),
      se_ub = numeric(),
      n_valid = integer()
    ))
  }

  if (is.null(tau_assumptions)) {
    tau_assumptions <- unique(unlist(lapply(
      relative_effect_list,
      function(x) x$tau_assumption
    )))
  }

  out <- lapply(tau_assumptions, function(assump) {
    lb_vals <- ub_vals <- numeric()

    for (tab in relative_effect_list) {
      row <- tab[tab$tau_assumption == assump, , drop = FALSE]
      if (nrow(row) == 0L) next
      if (!is.finite(row$effect_lb[1]) || !is.finite(row$effect_ub[1])) next
      lb_vals <- c(lb_vals, row$effect_lb[1])
      ub_vals <- c(ub_vals, row$effect_ub[1])
    }

    tibble::tibble(
      tau_assumption = assump,
      se_lb = if (length(lb_vals) > 1L) stats::sd(lb_vals) else NA_real_,
      se_ub = if (length(ub_vals) > 1L) stats::sd(ub_vals) else NA_real_,
      n_valid = length(lb_vals)
    )
  })

  dplyr::bind_rows(out)
}


.FinalizeRelativeEffects <- function(relative_effects, re_summary, alpha = 0.05) {
  out <- relative_effects
  out$lower_ci <- NA_real_
  out$upper_ci <- NA_real_
  out$c_crit <- NA_real_
  out$n_boot_valid <- 0L

  for (i in seq_len(nrow(out))) {
    row_summary <- re_summary[
      re_summary$tau_assumption == out$tau_assumption[i],
      ,
      drop = FALSE
    ]
    if (nrow(row_summary) == 0L) next

    out$n_boot_valid[i] <- row_summary$n_valid[1]
    if (out$status[i] != "ok") next

    rel_ci <- tryCatch(
      imbens_manski_ci(
        lb_hat = out$effect_lb[i],
        ub_hat = out$effect_ub[i],
        se_lb = row_summary$se_lb[1],
        se_ub = row_summary$se_ub[1],
        alpha = alpha,
        truncate_diff = TRUE
      ),
      error = function(e) list(ci = c(NA_real_, NA_real_), c_crit = NA_real_)
    )

    out$lower_ci[i] <- rel_ci$ci[1]
    out$upper_ci[i] <- rel_ci$ci[2]
    out$c_crit[i] <- rel_ci$c_crit
  }

  out
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
