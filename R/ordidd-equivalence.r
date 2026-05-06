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
#' @importFrom tibble tibble
#' @importFrom stats quantile dnorm
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
#' @param equivalence_result An output from \code{\link{run_equivalence}}.
#' @importFrom ggplot2 ggplot aes geom_hline geom_line labs annotate theme_minimal
#' @importFrom rlang .data
#' @export
plot_equivalence <- function(equivalence_result) {
  df_plot <- equivalence_result$estimates

  # Rejection threshold
  d_max <- max(abs(df_plot$upper_ci), abs(df_plot$lower_ci))

  # Worst case bias with the d_max
  worst_case_bias <- 2 * d_max / equivalence_result$M

  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = .data$v)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') +
    ggplot2::geom_line(ggplot2::aes(y = .data$upper_ci), linetype = "dashed") +
    ggplot2::geom_line(ggplot2::aes(y = .data$lower_ci), linetype = "dashed") +
    ggplot2::geom_line(
      ggplot2::aes(y = .data$point_estimate),
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
        "Rejection Threshold: +/-",
        round(d_max, 3),
        " ",
        "(|Bias| <= ",
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


#' Equivalence test based on category-specific effects
#'
#' Implements the category-specific equivalence test described in
#' the appendix section ``Equivalence test with category specific effects''.
#' Treats two pre-treatment periods as if they were the pre and post periods,
#' estimates \eqn{\widetilde{\zeta}_j = \Pr(Y_{i1} = j \mid D_i = 1) - \Pr(Y^N_{11} = j)}
#' under the parametric (location-scale) model, and returns bootstrap
#' confidence intervals at the requested level.
#'
#' @param data A data frame containing the data.
#' @param outcome Character; column name for the outcome variable.
#' @param post Character; column name for the post-period indicator (logical).
#'   In the equivalence-test setup, \code{post = TRUE} corresponds to the
#'   second pre-treatment period (e.g., 2012) and \code{post = FALSE} to the
#'   first pre-treatment period (e.g., 2010).
#' @param treat Character; column name for the (eventual) treatment indicator.
#' @param cluster Character; column name for the cluster id used for block
#'   bootstrap. If \code{NULL}, bootstrap is at the individual level.
#' @param n_boot Number of bootstrap iterations. Default \code{500}.
#' @param use_parallel Logical; if \code{TRUE}, run the bootstrap in parallel.
#' @param alpha Significance level for the point-wise CI.
#'   Default \code{0.10} (i.e., a 90\% CI).
#' @param link Character; \code{"probit"} (default) or \code{"logit"}.
#' @return A list with
#'   \item{estimates}{A tibble with one row per category containing the point
#'     estimate, bootstrap SE, and lower/upper CI bounds.}
#'   \item{boot_estimates}{The matrix of bootstrap zeta estimates
#'     (\code{n_boot} rows by \code{J} columns).}
#'   \item{alpha}{The significance level used to construct the CI.}
#' @importFrom dplyr sym
#' @importFrom rlang !! sym
#' @importFrom foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom tibble tibble
#' @importFrom stats quantile sd
#' @export
run_equivalence_zeta <- function(
  data,
  outcome,
  post,
  treat,
  cluster,
  n_boot = 500,
  use_parallel = FALSE,
  alpha = 0.10,
  link = c("probit", "logit")
) {
  link <- match.arg(link)

  df_format <- .orddid_extract_Y_cells(
    data = data,
    outcome = outcome,
    post = post,
    treat = treat
  )

  # Number of outcome categories J. Used downstream to drop bootstrap
  # iterations whose resampled cells happen to omit one of the J categories
  # (those iterations would otherwise mis-specify the cutoffs and contaminate
  # the percentile CI). See PR #8 for the matching change in the parametric
  # path.
  n_categories <- length(unique(c(
    df_format$Y00, df_format$Y01, df_format$Y10, df_format$Y11
  )))

  # Point estimate of zeta_j
  point_estimate <- .GetPointEstimateZeta(df = df_format, link = link)

  # Bootstrap setup
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

  boot_estimate <- foreach::foreach(
    i = seq_len(n_boot),
    .combine = "rbind",
    .inorder = FALSE,
    .packages = c("dplyr", "rlang", "orddid")
  ) %dopar%
    {
      tryCatch(
        .RunBootstrapZeta(
          df = data,
          cluster = cluster,
          outcome = outcome,
          post = post,
          treat = treat,
          link = link,
          n_categories = n_categories
        ),
        error = function(e) rep(NA_real_, length(point_estimate))
      )
    }

  # Drop failed bootstrap rows
  valid_rows <- stats::complete.cases(boot_estimate)
  if (sum(valid_rows) == 0L) {
    stop("All bootstrap iterations failed.", call. = FALSE)
  }
  boot_clean <- boot_estimate[valid_rows, , drop = FALSE]

  J <- length(point_estimate)
  ci_lower <- apply(
    boot_clean, 2, stats::quantile, probs = alpha / 2, na.rm = TRUE
  )
  ci_upper <- apply(
    boot_clean, 2, stats::quantile, probs = 1 - alpha / 2, na.rm = TRUE
  )
  se_zeta <- apply(boot_clean, 2, stats::sd, na.rm = TRUE)

  df_out <- tibble::tibble(
    category = seq_len(J),
    estimate = point_estimate,
    se = se_zeta,
    lower_ci = ci_lower,
    upper_ci = ci_upper
  )

  return(list(
    estimates = df_out,
    boot_estimates = boot_clean,
    alpha = alpha,
    n_boot_valid = sum(valid_rows),
    n_boot_fail = sum(!valid_rows)
  ))
}


#' Compute zeta_j point estimates for given pre-period cells
#'
#' Same logic as \code{.GetPointEstimate}, but returns only the
#' category-specific differences (\code{zeta}) so the bootstrap loop
#' can stack them into a matrix.
#' @noRd
.GetPointEstimateZeta <- function(df, link = c("probit", "logit")) {
  link <- match.arg(link)
  y1_obs <- .EstimateObservedDist(Y11 = df$Y11, link = link)
  y1_cf <- .EstimateCounterfactualDist(
    Y10 = df$Y10,
    Y01 = df$Y01,
    Y00 = df$Y00,
    link = link
  )
  zeta <- .ComputeTreatmentEffects(
    y1_prop = y1_obs$prob,
    y0_prop = y1_cf$prob
  )
  return(zeta)
}


#' Single bootstrap iteration for the zeta-based equivalence test
#' @noRd
.RunBootstrapZeta <- function(
  df,
  cluster,
  outcome,
  post,
  treat,
  link = c("probit", "logit"),
  n_categories = NULL
) {
  link <- match.arg(link)
  df_boot <- .SampleDf(df = df, cluster = cluster)
  df_format <- .orddid_extract_Y_cells(
    data = df_boot,
    outcome = outcome,
    post = post,
    treat = treat
  )
  if (!is.null(n_categories)) {
    if (any(vapply(df_format, function(y) length(unique(y)),
                   FUN.VALUE = integer(1)) < n_categories)) {
      stop("bootstrap sample is missing one or more outcome categories",
           call. = FALSE)
    }
  }
  .GetPointEstimateZeta(df = df_format, link = link)
}


#' Plot the category-specific equivalence test result
#'
#' Plots the point estimates of \eqn{\widetilde{\zeta}_j} with their
#' bootstrap confidence intervals and overlays two horizontal equivalence
#' thresholds derived from the post-treatment effect estimates.
#'
#' @param equivalence_result Output from \code{\link{run_equivalence_zeta}}.
#' @param post_zeta Numeric vector of post-treatment effect estimates
#'   \eqn{\widehat{\zeta}_j}. Used to compute the equivalence thresholds
#'   \eqn{\epsilon_{\max} = \max_j |\widehat{\zeta}_j|} and
#'   \eqn{\epsilon_{\min} = \min_j |\widehat{\zeta}_j|}.
#' @param category_labels Optional character vector of length \code{J}
#'   to label the x-axis. Defaults to category indices.
#' @param title Optional plot title.
#' @param point_color Color for the zeta point estimates and CI bars.
#'   Defaults to the package's standard 25-mile color.
#' @return A \code{ggplot2} object.
#' @importFrom ggplot2 ggplot aes geom_hline geom_pointrange labs theme_bw
#'   scale_x_continuous scale_x_discrete annotate geom_point geom_errorbar
#' @importFrom rlang .data
#' @export
plot_equivalence_zeta <- function(
  equivalence_result,
  post_zeta = NULL,
  category_labels = NULL,
  title = NULL,
  point_color = "#006284"
) {
  df_plot <- equivalence_result$estimates
  alpha <- equivalence_result$alpha
  ci_pct <- round(100 * (1 - alpha))

  J <- nrow(df_plot)
  if (is.null(category_labels)) {
    category_labels <- as.character(df_plot$category - 1L)
  }
  df_plot$category_label <- factor(
    df_plot$category,
    levels = df_plot$category,
    labels = category_labels
  )

  p <- ggplot2::ggplot(
    df_plot,
    ggplot2::aes(x = .data$category_label, y = .data$estimate)
  ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey60")

  if (!is.null(post_zeta)) {
    eps_max <- max(abs(post_zeta), na.rm = TRUE)
    eps_min <- min(abs(post_zeta), na.rm = TRUE)
    # Draw thresholds underneath the points so the CI bars stay on top.
    p <- p +
      ggplot2::geom_hline(
        yintercept = c(-eps_max, eps_max),
        linetype = "solid",
        color = "grey50",
        linewidth = 0.9
      ) +
      ggplot2::geom_hline(
        yintercept = c(-eps_min, eps_min),
        linetype = "dashed",
        color = "grey50",
        linewidth = 0.9
      )
  }

  p +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$lower_ci, ymax = .data$upper_ci),
      width = 0.15,
      linewidth = 0.6,
      color = point_color
    ) +
    ggplot2::geom_point(size = 2.5, color = point_color) +
    ggplot2::labs(
      x = "Outcome category",
      y = expression(tilde(zeta)[j]),
      title = title
    ) +
    ggplot2::theme_bw()
}
