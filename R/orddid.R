#' Ordinal Difference-in-Differences for Panel Data
#'
#' \code{ord_did()} implements the difference-in-differences for the ordinal outcome.
#'
#' @param data A data frame containing the variables.
#' @param outcome Character scalar; column name for the outcome variable.
#' @param post Character scalar; column name for the post-period indicator.
#' @param treat Character scalar; column name for the treatment indicator.
#' @param cluster Character scalar; column name for the cluster id.
#' @param n_boot The number of bootstrap iterations. Default is \code{n_boot = 500}.
#' @param use_parallel If \code{TRUE}, use parallel computing for bootstrap.
#' @param method Character; \code{"parametric"} (default) for ordered probit
#'   without covariates, \code{"semiparametric"} for NPMLE-based semiparametric
#'   estimation with covariates, or \code{"parametric_covariates"} for ordered
#'   probit with covariates (uses normal CDF instead of semiparametric F).
#' @param covariates Character vector of covariate column names (required for
#'   \code{method = "semiparametric"} and \code{method = "parametric_covariates"}).
#'   The first covariate should be a continuous variable with positive Lebesgue
#'   density (no intercept is added; the cutpoints handle the baseline level).
#' @param fixed_index Integer; which column of the covariate matrix to normalize
#'   in the semiparametric Step 1 (default 1, i.e., the first covariate).
#' @param n_cores Integer; number of cores for parallel computing. Default
#'   \code{NULL} uses all available cores minus one.
#' @param tau_assumption Character vector controlling how the relative-effect
#'   bounds are computed. \code{"none"} uses the unconstrained sharp bounds of
#'   Lu, Zhang, and Ding (2020). \code{"stochastic_dominance"} imposes positive
#'   regression dependence on the joint distribution of the post-treatment
#'   potential outcomes and computes the tighter feasible bounds by linear
#'   programming. By default, both are returned.
#' @param semiparametric_regularization Character; regularization used for the
#'   first-stage error CDF when \code{method = "semiparametric"}. The default
#'   \code{"ispline"} uses monotone I-spline smoothing; \code{"ridge"} uses the
#'   vanishing-ridge regularization. Ignored for other methods.
#' @return A list with components:
#' \item{inputs}{Function arguments for reproducibility.}
#' \item{estimate_effects}{Data frame of category-specific effects with CIs.}
#' \item{relative_effects}{Data frame of relative treatment effect bounds with CIs.}
#' @examples
#'\donttest{
#' ## load packages
#' library(orddid)
#' library(dplyr)
#'
#' ## load example data
#' data("gun_twowave")
#'
#' ## run
#' ## fit the ordinal DID
#' set.seed(1234)
#' treated_id <- gun_twowave %>%
#'   filter(treat_25mi == 1) %>%
#'   select(caseid, treated_group = treat_25mi)
#' dat_use <- gun_twowave %>%
#'   left_join(treated_id, by = "caseid") %>%
#'   mutate(
#'     post = ifelse(year == 2012, TRUE, FALSE),
#'     treated = ifelse(!is.na(treated_group), TRUE, FALSE)
#'   )
#'
#' dat_use %>%
#'   group_by(pid3) %>%
#'   group_map(
#'     ~ ord_did(
#'       data = .x,
#'       outcome = "guns",
#'       post = "post",
#'       treat = "treated",
#'       cluster = "caseid"
#'     )
#'   )
#'
#' ## view summary of the output
#' ## non-cumulative effects
#' summary(fit, cumulative = FALSE)
#'
#' ## cumulative effects
#' summary(fit)
#' }
#' @importFrom tibble tibble
#' @export
ord_did <- function(
  data,
  outcome,
  post,
  treat,
  cluster,
  n_boot = 500,
  use_parallel = FALSE,
  method = c("parametric", "semiparametric", "parametric_covariates"),
  covariates = NULL,
  fixed_index = 1,
  n_cores = NULL,
  tau_assumption = c("none", "stochastic_dominance"),
  semiparametric_regularization = c("ispline", "ridge")
) {
  method <- match.arg(method)
  semiparametric_regularization <- match.arg(
    semiparametric_regularization
  )
  tau_assumption <- .orddid_match_tau_assumptions(tau_assumption)

  if (method == "semiparametric" || method == "parametric_covariates") {
    return(.ord_did_with_covariates(
      data = data,
      outcome = outcome,
      post = post,
      treat = treat,
      cluster = cluster,
      covariates = covariates,
      n_boot = n_boot,
      use_parallel = use_parallel,
      fixed_index = fixed_index,
      n_cores = n_cores,
      method = method,
      tau_assumption = tau_assumption,
      semiparametric_regularization = semiparametric_regularization
    ))
  }

  ## ---- Parametric (existing behavior) ----
  df_format <- .orddid_extract_Y_cells(
    data = data,
    outcome = outcome,
    post = post,
    treat = treat
  )

  # Compute point estimate
  point_estimate <- .GetPointEstimate(
    df = df_format,
    tau_assumption = tau_assumption
  )

  if (use_parallel) {
    if (is.null(n_cores)) n_cores <- max(1, parallel::detectCores() - 1)
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

  # Compute bootstrap estimates
  boot_estimate <- foreach::foreach(
    i = seq_len(n_boot),
    .combine = dplyr::bind_rows,
    .inorder = FALSE,
    .packages = c("dplyr", "rlang", "orddid")
  ) %dopar%
    {
      est <- tryCatch(
        .RunBootstrap(
          df = data,
          cluster = cluster,
          outcome = outcome,
          post = post,
          treat = treat,
          tau_assumption = tau_assumption
        ),
        error = function(e) NULL
      )
      if (is.null(est)) {
        return(tibble::tibble(
          iteration = i,
          diff_effects = list(NULL),
          relative_effect = list(NULL)
        ))
      }
      tibble::tibble(
        iteration = i,
        diff_effects = list(est$diff_effects),
        relative_effect = list(est$relative_effect)
      )
    }

  ci_diff <- .ComputeCI(boot_estimates = boot_estimate, alpha = 0.05)

  # Summarize results
  diff_res <- tibble(
    category = 1:length(point_estimate$estimated_effects$diff_effects),
    effect = point_estimate$estimated_effects$diff_effects,
    lower_ci = ci_diff$lower,
    upper_ci = ci_diff$upper
  )

  relative_res <- .FinalizeRelativeEffects(
    relative_effects = point_estimate$estimated_effects$relative_effect,
    re_summary = ci_diff$re_summary,
    alpha = 0.05
  )
  return(list(
    inputs = list(
      outcome = outcome,
      post = post,
      treat = treat,
      cluster = cluster,
      n_boot = n_boot,
      method = "parametric",
      tau_assumption = tau_assumption
    ),
    estimate_effects = diff_res,
    relative_effects = relative_res,
    n_boot_fail = ci_diff$n_fail
  ))
}

## ---- Semiparametric path ----

.ord_did_with_covariates <- function(data, outcome, post, treat, cluster,
                                     covariates, n_boot, use_parallel,
                                     fixed_index, n_cores = NULL,
                                     method = "semiparametric",
                                     tau_assumption = "none",
                                     semiparametric_regularization = c(
                                       "ispline", "ridge"
                                     )) {
  tau_assumption <- .orddid_match_tau_assumptions(tau_assumption)
  semiparametric_regularization <- match.arg(
    semiparametric_regularization
  )

  if (is.null(covariates)) {
    stop("covariates must be specified for method = '", method, "'.")
  }
  if (!is.character(covariates) || length(covariates) < 1L ||
      anyNA(covariates) || any(!nzchar(covariates)) ||
      anyDuplicated(covariates)) {
    stop("covariates must be a non-empty character vector of unique column names.")
  }

  required_columns <- unique(c(outcome, post, treat, cluster, covariates))
  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns) > 0L) {
    stop(
      "Missing required columns: ",
      paste(missing_columns, collapse = ", "),
      "."
    )
  }

  if (length(n_boot) != 1L || !is.numeric(n_boot) || is.na(n_boot) ||
      !is.finite(n_boot) || n_boot < 0 || n_boot != floor(n_boot)) {
    stop("n_boot must be one non-negative integer.")
  }
  n_boot <- as.integer(n_boot)

  if (length(fixed_index) != 1L || !is.numeric(fixed_index) ||
      is.na(fixed_index) || !is.finite(fixed_index) ||
      fixed_index != floor(fixed_index) || fixed_index < 1L ||
      fixed_index > length(covariates)) {
    stop("fixed_index must identify one of the supplied covariates.")
  }
  fixed_index <- as.integer(fixed_index)

  ## Use one complete-case sample for both the point estimate and bootstrap.
  ## In particular, do not let split() silently discard missing cluster IDs.
  keep <- stats::complete.cases(data[, required_columns, drop = FALSE])
  data_use <- data[keep, , drop = FALSE]
  if (nrow(data_use) == 0L) {
    stop("No complete observations remain for covariate estimation.")
  }

  x_cols <- covariates
  X_mat <- as.matrix(data_use[, x_cols, drop = FALSE])
  if (!is.numeric(X_mat) || any(!is.finite(X_mat))) {
    stop("All covariates must be numeric and finite on the estimation sample.")
  }

  y_vec <- suppressWarnings(as.integer(data_use[[outcome]]))
  if (anyNA(y_vec)) {
    stop("outcome must be coercible to ordered integer categories.")
  }
  y_vec <- y_vec - min(y_vec)
  if (!identical(sort(unique(y_vec)), 0:2)) {
    stop(
      "Covariate estimators currently require exactly three pooled outcome ",
      "categories, coded consecutively."
    )
  }

  p_vec <- as.logical(data_use[[post]])
  t_vec <- as.logical(data_use[[treat]])
  if (anyNA(p_vec) || anyNA(t_vec)) {
    stop("post and treat must be coercible to non-missing logical indicators.")
  }

  ## Cell membership: 1=Y00, 2=Y01, 3=Y10, 4=Y11.
  cell_id <- 1L + as.integer(p_vec) + 2L * as.integer(t_vec)
  cell_sizes <- tabulate(cell_id, nbins = 4L)
  if (any(cell_sizes == 0L)) {
    stop("All four treatment-by-period cells must contain observations.")
  }

  i00 <- cell_id == 1L
  i01 <- cell_id == 2L
  i10 <- cell_id == 3L
  i11 <- cell_id == 4L

  if (!identical(sort(unique(y_vec[i00])), 0:2)) {
    stop("The control-pre cell (Y00) must contain all three outcome categories 0:2.")
  }
  if (method == "semiparametric" &&
      length(unique(X_mat[i00, fixed_index])) < 2L) {
    stop("The covariate selected by fixed_index must vary in the control-pre cell.")
  }

  cells <- list(
    Y00 = y_vec[i00], X00 = X_mat[i00, , drop = FALSE],
    Y01 = y_vec[i01], X01 = X_mat[i01, , drop = FALSE],
    Y10 = y_vec[i10], X10 = X_mat[i10, , drop = FALSE],
    Y11 = y_vec[i11], X11 = X_mat[i11, , drop = FALSE]
  )

  ## Select estimation function based on method
  est_fn <- if (method == "semiparametric") {
    .semiparametric_did_once
  } else {
    .parametric_covariates_did_once
  }

  ## Pass the CDF regularization only to the semiparametric workhorse so the
  ## ordered-probit-with-covariates interface remains unchanged.
  fit_once <- function(Y00, X00, Y01, X01, Y10, X10, Y11, X11) {
    fit_args <- list(
      Y00 = Y00, X00 = X00,
      Y01 = Y01, X01 = X01,
      Y10 = Y10, X10 = X10,
      Y11 = Y11, X11 = X11,
      fixed_index = fixed_index,
      verbose = FALSE,
      tau_assumption = tau_assumption
    )
    if (method == "semiparametric") {
      fit_args$semiparametric_regularization <-
        semiparametric_regularization
    }
    do.call(est_fn, fit_args)
  }

  ## Point estimate
  point <- fit_once(
    Y00 = cells$Y00, X00 = cells$X00,
    Y01 = cells$Y01, X01 = cells$X01,
    Y10 = cells$Y10, X10 = cells$X10,
    Y11 = cells$Y11, X11 = cells$X11
  )

  J <- length(point$zeta)
  if (J != 3L || any(!is.finite(point$zeta))) {
    stop("The point estimate did not return three finite category effects.")
  }

  ## Bootstrap — precompute vectors and indices for fast resampling
  cluster_vec <- data_use[[cluster]]

  ## Cluster index: map each cluster to its row positions (integer vectors)
  cluster_idx <- split(seq_along(y_vec), cluster_vec)
  n_clusters  <- length(cluster_idx)

  ## Draw the complete bootstrap plan in the calling R process. Workers only
  ## consume these indices, so results do not depend on scheduling or n_cores.
  bootstrap_cluster_samples <- lapply(seq_len(n_boot), function(b) {
    sample.int(n_clusters, n_clusters, replace = TRUE)
  })

  ## Set up parallel backend
  if (use_parallel && is.null(n_cores)) {
    detected_cores <- parallel::detectCores()
    n_cores <- if (is.na(detected_cores)) 1L else max(1L, detected_cores - 1L)
  }
  if (use_parallel &&
      (length(n_cores) != 1L || !is.numeric(n_cores) || is.na(n_cores) ||
       !is.finite(n_cores) || n_cores < 1L || n_cores != floor(n_cores))) {
    stop("n_cores must be one positive integer when use_parallel = TRUE.")
  }
  if (use_parallel && n_boot > 0L) {
    n_cores <- as.integer(n_cores)
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    on.exit({
      parallel::stopCluster(cl)
      foreach::registerDoSEQ()
    }, add = TRUE)
  } else {
    foreach::registerDoSEQ()
  }

  ## Capture internal function for export to parallel workers
  ## (needed when using devtools::load_all instead of installed package)
  .boot_covariate_fn <- fit_once

  ## Run cluster bootstrap in parallel
  boot_results <- if (n_boot == 0L) {
    list()
  } else {
    foreach::foreach(
      b = seq_len(n_boot),
      .inorder = TRUE,
      .packages = c("orddid")
    ) %dopar% {
      sampled <- bootstrap_cluster_samples[[b]]
      row_idx <- unlist(cluster_idx[sampled], use.names = FALSE)

      y_b    <- y_vec[row_idx]
      X_b    <- X_mat[row_idx, , drop = FALSE]
      cell_b <- cell_id[row_idx]

      b00 <- cell_b == 1L
      b01 <- cell_b == 2L
      b10 <- cell_b == 3L
      b11 <- cell_b == 4L

      if (min(sum(b00), sum(b01), sum(b10), sum(b11)) < 10L ||
          !identical(sort(unique(y_b[b00])), 0:2)) {
        return(NULL)
      }

      fit_b <- tryCatch(
        .boot_covariate_fn(
          Y00 = y_b[b00], X00 = X_b[b00, , drop = FALSE],
          Y01 = y_b[b01], X01 = X_b[b01, , drop = FALSE],
          Y10 = y_b[b10], X10 = X_b[b10, , drop = FALSE],
          Y11 = y_b[b11], X11 = X_b[b11, , drop = FALSE]
        ),
        error = function(e) NULL
      )

      if (is.null(fit_b) || length(fit_b$zeta) != J ||
          any(!is.finite(fit_b$zeta))) {
        return(NULL)
      }
      fit_b
    }
  }

  ## Collect results
  boot_zeta <- matrix(NA, n_boot, J)
  boot_tau <- vector("list", n_boot)
  n_fail <- 0

  for (b in seq_len(n_boot)) {
    if (!is.null(boot_results[[b]])) {
      boot_zeta[b, ] <- boot_results[[b]]$zeta
      boot_tau[[b]] <- boot_results[[b]]$tau
    } else {
      n_fail <- n_fail + 1
    }
  }

  ## Compute CIs
  valid <- if (n_boot == 0L) {
    logical(0)
  } else {
    apply(boot_zeta, 1L, function(z) all(is.finite(z)))
  }

  zeta_lower <- zeta_upper <- rep(NA, J)
  if (sum(valid) > 10) {
    for (j in seq_len(J)) {
      qq <- stats::quantile(boot_zeta[valid, j], c(0.025, 0.975))
      zeta_lower[j] <- qq[1]
      zeta_upper[j] <- qq[2]
    }
  }

  diff_res <- tibble::tibble(
    category = seq_len(J),
    effect   = point$zeta,
    lower_ci = zeta_lower,
    upper_ci = zeta_upper
  )

  relative_res <- .FinalizeRelativeEffects(
    relative_effects = point$tau,
    re_summary = .SummariseRelativeBootstrap(
      boot_tau[!vapply(boot_tau, is.null, logical(1))],
      tau_assumptions = point$tau$tau_assumption
    ),
    alpha = 0.05
  )

  list(
    inputs = list(
      outcome = outcome,
      post = post,
      treat = treat,
      cluster = cluster,
      covariates = covariates,
      n_boot = n_boot,
      method = method,
      tau_assumption = tau_assumption,
      fixed_index = fixed_index,
      semiparametric_regularization = semiparametric_regularization,
      use_parallel = use_parallel,
      n_cores = n_cores
    ),
    estimate_effects = diff_res,
    relative_effects = relative_res,
    n_boot_fail = n_fail,
    bootstrap_effects = boot_zeta
  )
}

#' Extract cell data for semiparametric estimation
#' @keywords internal
.extract_semiparametric_cells <- function(data, outcome, post, treat, x_cols) {
  y <- data[[outcome]]
  p <- as.logical(data[[post]])
  t <- as.logical(data[[treat]])

  keep <- !(is.na(y) | is.na(p) | is.na(t))
  y <- y[keep]; p <- p[keep]; t <- t[keep]
  data_k <- data[keep, , drop = FALSE]

  ## Recode Y to 0-indexed (0, 1, ..., J-1) if it starts from a higher value
  y <- as.integer(y) - min(as.integer(y))

  ## Build design matrix (no intercept; cutpoints handle the baseline level,
  ## and in the semiparametric model F(0) is estimated nonparametrically)
  X <- as.matrix(data_k[, x_cols, drop = FALSE])

  i00 <- !p & !t
  i01 <-  p & !t
  i10 <- !p &  t
  i11 <-  p &  t

  list(
    Y00 = y[i00], X00 = X[i00, , drop = FALSE],
    Y01 = y[i01], X01 = X[i01, , drop = FALSE],
    Y10 = y[i10], X10 = X[i10, , drop = FALSE],
    Y11 = y[i11], X11 = X[i11, , drop = FALSE]
  )
}

.GetPointEstimate <- function(df,
                              tau_assumption = "none") {
  tau_assumption <- .orddid_match_tau_assumptions(tau_assumption)
  y1_obs <- .EstimateObservedDist(Y11 = df$Y11)
  y1_cf <- .EstimateCounterfactualDist(
    Y10 = df$Y10,
    Y01 = df$Y01,
    Y00 = df$Y00
  )

  # Compute effects
  diff_effects <- .ComputeTreatmentEffects(
    y1_prop = y1_obs$prob,
    y0_prop = y1_cf$prob
  )
  relative <- .ComputeRelativeEffect(
    p1 = y1_obs$prob,
    p0 = y1_cf$prob,
    tau_assumption = tau_assumption
  )
  return(
    list(
      estimated_props = list(
        y1_obs = y1_obs$prob,
        y1_cf = y1_cf$prob
      ),
      estimated_effects = list(
        diff_effects = diff_effects,
        relative_effect = relative
      )
    )
  )
}



.EstimateCounterfactualDist <- function(Y10, Y01, Y00,
                                        link = c("probit", "logit")) {
  link <- match.arg(link)
  # Estimate parameters for the observed outcome
  fit_obs <- .EstimateControlParams(
    Y00 = Y00,
    Y01 = Y01,
    Y10 = Y10,
    link = link
  )

  # Estimate parameters for the counterfactual outcome
  fit_cf <- .EstimateCounterfactualParams(fit = fit_obs)

  # Estimate the outcome proportions
  prob_y11 <- .EstimateOutcomeProportions(
    mu = fit_cf$mu,
    sd = fit_cf$sd,
    cutoff = fit_cf$cutoff,
    link = link
  )

  return(list(
    params = fit_cf,
    prob = prob_y11
  ))
}

.EstimateObservedDist <- function(Y11, link = c("probit", "logit")) {
  link <- match.arg(link)
  # Estimate parameters for the observed outcome
  fit_obs <- .EstimateLatentOutcomeParams(
    Ydt = Y11,
    cutoff = NULL,
    link = link
  )

  # Estimate the outcome proportions
  prob_y11 <- .EstimateOutcomeProportions(
    mu = fit_obs$mu,
    sd = fit_obs$sd,
    cutoff = fit_obs$cutoff,
    link = link
  )

  return(list(
    params = fit_obs,
    prob = prob_y11
  ))
}
