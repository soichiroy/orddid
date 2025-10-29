
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
  estimate <- .GetPointEstimate(df = df_format)

  return(estimate)
}


#' Sampling function
#' @noRd
#' @importFrom dplyr group_by slice_sample ungroup across all_of
.SampleDf <- function(df, cluster) {
  df_boot <- df %>%
    group_by(across(all_of(cluster))) %>%
    slice_sample(prop = 1) %>%
    ungroup()
  return(df_boot)
}

.ComputeCI <- function(boot_estimates, alpha = 0.05) {
  # Collect bootstrap estimates
  boot_diff_effects <- sapply(boot_estimates, function(x) {
    x$estimated_effects$diff_effects
  })

  # Compute percentile CI
  lower_bound <- apply(boot_diff_effects, 1, quantile, probs = alpha / 2)
  upper_bound <- apply(boot_diff_effects, 1, quantile, probs = 1 - alpha / 2)

  return(list(lower = lower_bound, upper = upper_bound))
}