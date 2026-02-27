#' Run equivalence test
#'
#' @param data A data frame containing the data.
#' @param outcome Character scalar; column name for the outcome variable.
#' @param post Character scalar; column name for the post-period indicator (logical).
#' @param treat Character scalar; column name for the treatment indicator (logical).
#' @param cluster Character scalar; column name for the cluster id.
#'   If left as \code{NULL}, bootstrap is implemented at the individual level.
#' @param n_boot The number of boostrapt iterations for estimating the variance. Default is \code{n_boot = 500}.
#' @param use_parallel Logical; if TRUE, use parallel computing for bootstrap.
#' @param v_range Numeric vector; the range of quantiles to evaluate the equivalence test.
#'
#' @importFrom dplyr sym group_by ungroup slice_sample
#' @importFrom rlang !! sym
#' @importFrom foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster detectCores
#' @export
run_equivalence <- function(
  data,
  outcome,
  post,
  treat,
  cluster,
  n_boot = 500,
  use_parallel = FALSE,
  v_range = seq(0, 1, by = 0.01)
) {
  df_format <- .orddid_extract_Y_cells(
    data = data,
    outcome = outcome,
    post = post,
    treat = treat
  )

  # Compute point estimate
  point_estimate <- .GetPointEstimatePre(
    df = df_format,
    v_range = v_range
  )

  # Bootstrap with optional parallelization using foreach/doParallel
  has_foreach <- requireNamespace("foreach", quietly = TRUE)
  has_doParallel <- requireNamespace("doParallel", quietly = TRUE)
  has_parallel <- requireNamespace("parallel", quietly = TRUE)

  if (!(has_foreach && has_doParallel && has_parallel)) {
    use_parallel <- FALSE
  }

  if (use_parallel) {
    n_cores <- max(1, parallel::detectCores() - 1)
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    on.exit(
      {
        parallel::stopCluster(cl)
        foreach::registerDoSEQ()
      },
      add = TRUE
    )
  } else {
    foreach::registerDoSEQ()
  }

  # Bootstrap iterations
  boot_estimate <- foreach::foreach(
    i = seq_len(n_boot),
    .combine = "rbind",
    .inorder = FALSE,
    .packages = c("dplyr", "rlang", "orddid")
  ) %dopar%
    {
      .RunBootstrapPre(
        df = data,
        cluster = cluster,
        outcome = outcome,
        post = post,
        treat = treat,
        v_range = v_range
      )
    }

  # Compute point-wise CI
  ci_upper_pw <- apply(boot_estimate, 2, function(x) {
    quantile(x, probs = 0.95)
  })
  ci_lower_pw <- apply(boot_estimate, 2, function(x) {
    quantile(x, probs = 0.05)
  })

  df_out <- tibble(
    v = v_range,
    point_estimate = point_estimate$r,
    lower_ci = ci_lower_pw,
    upper_ci = ci_upper_pw
  )
  return(list(
    estimates = df_out,
    M = point_estimate$M,
    inf_factor = 1 / point_estimate$M
  ))
}

#' Estimate parameters for the pre-treatment data
#' @noRd
.GetPointEstimatePre <- function(df, v_range = seq(0, 1, by = 0.01)) {
  # D = 0 group
  fit_y00 <- .EstimateLatentOutcomeParams(Ydt = df$Y00, cutoff = NULL)
  fit_y01 <- .EstimateLatentOutcomeParams(Ydt = df$Y01, cutoff = fit_y00$cutoff)
  q0 <- .ComputeQFunction(fit0 = fit_y00, fit1 = fit_y01, v_range = v_range)

  # D = 1 group (use shared cutoffs from Y00)
  fit_y10 <- .EstimateLatentOutcomeParams(Ydt = df$Y10, cutoff = fit_y00$cutoff)
  fit_y11 <- .EstimateLatentOutcomeParams(Ydt = df$Y11, cutoff = fit_y00$cutoff)
  q1 <- .ComputeQFunction(fit0 = fit_y10, fit1 = fit_y11, v_range = v_range)

  # Estimate M
  M <- .EstimateLowerConst(
    fit0 = fit_y00,
    fit1 = fit_y01,
    v_range = v_range
  )

  return(list(q0 = q0, q1 = q1, r = q1 - q0, M = M))
}

#' Compute the quantile function values for given fit parameters
#' @noRd
#' @importFrom stats pnorm qnorm
.ComputeQFunction <- function(fit0, fit1, v_range) {
  x <- pnorm((fit1$mu - fit0$mu) / fit0$sd + (fit1$sd / fit0$sd) * qnorm(v_range))
  return(x)
}

#' Run a single bootstrap iteration for pre-treatment data
#'
#' @importFrom dplyr sym group_by ungroup slice_sample
#' @importFrom rlang !! sym
#' @noRd
.RunBootstrapPre <- function(
  df,
  cluster,
  outcome,
  post,
  treat,
  v_range
) {
  # Block bootstrap: resample clusters with replacement
  df_boot <- .SampleDf(df = df, cluster = cluster)

  # Compute bootstrap estimate
  df_format <- .orddid_extract_Y_cells(
    data = df_boot,
    outcome = outcome,
    post = post,
    treat = treat
  )

  boot_estimate <- .GetPointEstimatePre(df = df_format, v_range = v_range)$r
  return(boot_estimate)
}

#' Compute the identified bounds for M
#' @noRd
.EstimateLowerConst <- function(fit0, fit1, v_range) {
  b <- fit1$sd / fit0$sd
  a <- (fit1$mu - fit0$mu) / fit0$sd
  x <- b * dnorm(a + b * qnorm(v_range)) / dnorm(qnorm(v_range))
  x <- x[is.finite(x)]
  return(min(x))
}


#' Plot equivalence test result
#' @importFrom ggplot2 ggplot aes geom_hline geom_line labs annotate theme_minimal
#' @export
plot_equivalence <- function(equivalence_result) {
  df_plot <- equivalence_result$estimates

  # Rejection threshold
  d_max <- max(abs(df_plot$upper_ci), abs(df_plot$lower_ci))

  # Worst case bias with the d_max
  worst_case_bias <- 2 * d_max / equivalence_result$M

  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = v)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') +
    ggplot2::geom_line(aes(y = upper_ci), linetype = "dashed") +
    ggplot2::geom_line(aes(y = lower_ci), linetype = "dashed") +
    ggplot2::geom_line(
      ggplot2::aes(y = point_estimate),
      color = "#006284",
      linewidth = 1.2
    ) +
    ggplot2::labs(
      title = "Equivalence Test Result",
      x = "Quantile (v)",
      y = "r(v) = q1(v) - q0(v)"
    ) +
    ggplot2::geom_hline(
      yintercept = c(-d_max, d_max),
      linetype = "dotted",
      color = "#990000"
    ) +
    ggplot2::annotate(
      "text",
      x = 0.5,
      y = d_max,
      label = paste0(
        "Rejection Threshold: ±",
        round(d_max, 3),
        " ",
        "(|Bias| ≤ ",
        round(worst_case_bias, 3),
        ")"
      ),
      color =  "#990000",
      vjust = -0.5
    ) +
    ggplot2::ylim(-1.1 * d_max, 1.1 * d_max) +
    ggplot2::theme_bw()

  return(p)
}
