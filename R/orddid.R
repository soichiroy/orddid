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
  n_cores = NULL
) {
  method <- match.arg(method)

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
      method = method
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
  point_estimate <- .GetPointEstimate(df = df_format)

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
      est <- .RunBootstrap(
        df = data,
        cluster = cluster,
        outcome = outcome,
        post = post,
        treat = treat
      )
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

  # Relative effect results
  relative_ci <- imbens_manski_ci(
    lb_hat = point_estimate$estimated_effects$relative_effect[1],
    ub_hat = point_estimate$estimated_effects$relative_effect[2],
    se_lb = ci_diff$re_se[1],
    se_ub = ci_diff$re_se[2],
    alpha = 0.05,
    truncate_diff = TRUE
  )

  relative_res <- tibble(
    effect_lb = point_estimate$estimated_effects$relative_effect[1],
    effect_ub = point_estimate$estimated_effects$relative_effect[2],
    lower_ci = relative_ci$ci[1],
    upper_ci = relative_ci$ci[2],
    c_crit = relative_ci$c_crit
  )
  return(list(
    inputs = list(
      outcome = outcome,
      post = post,
      treat = treat,
      cluster = cluster,
      n_boot = n_boot,
      method = "parametric"
    ),
    estimate_effects = diff_res,
    relative_effects = relative_res
  ))
}

## ---- Semiparametric path ----

.ord_did_with_covariates <- function(data, outcome, post, treat, cluster,
                                     covariates, n_boot, use_parallel,
                                     fixed_index, n_cores = NULL,
                                     method = "semiparametric") {

  if (is.null(covariates)) {
    stop("covariates must be specified for method = '", method, "'.")
  }

  ## Select estimation function based on method
  est_fn <- if (method == "semiparametric") {
    .semiparametric_did_once
  } else {
    .parametric_covariates_did_once
  }

  x_cols <- covariates

  ## Build cell data
  cells <- .extract_semiparametric_cells(data, outcome, post, treat, x_cols)

  ## Point estimate
  point <- est_fn(
    Y00 = cells$Y00, X00 = cells$X00,
    Y01 = cells$Y01, X01 = cells$X01,
    Y10 = cells$Y10, X10 = cells$X10,
    Y11 = cells$Y11, X11 = cells$X11,
    fixed_index = fixed_index,
    verbose = FALSE
  )

  J <- length(point$zeta)

  ## Bootstrap — precompute vectors and indices for fast resampling

  ## Precompute: extract outcome, post, treat, and design matrix as raw vectors/matrices
  y_all <- data[[outcome]]
  p_all <- as.logical(data[[post]])
  t_all <- as.logical(data[[treat]])
  keep  <- !(is.na(y_all) | is.na(p_all) | is.na(t_all))

  y_vec       <- as.integer(y_all[keep]) - min(as.integer(y_all[keep]))
  X_mat       <- as.matrix(data[keep, x_cols, drop = FALSE])
  cluster_vec <- data[[cluster]][keep]

  ## Cell membership: 1=Y00, 2=Y01, 3=Y10, 4=Y11
  cell_id <- 1L + as.integer(as.logical(data[[post]][keep])) +
             2L * as.integer(as.logical(data[[treat]][keep]))

  ## Cluster index: map each cluster to its row positions (integer vectors)
  cluster_idx <- split(seq_along(y_vec), cluster_vec)
  n_clusters  <- length(cluster_idx)

  ## Set up parallel backend
  if (use_parallel) {
    if (is.null(n_cores)) n_cores <- max(1, parallel::detectCores() - 1)
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
  .boot_semi_fn <- est_fn

  ## Run cluster bootstrap in parallel
  boot_results <- foreach::foreach(
    b = seq_len(n_boot),
    .inorder = FALSE,
    .packages = c("orddid")
  ) %dopar% {
    ## Resample clusters with replacement, gather row indices
    sampled <- sample.int(n_clusters, n_clusters, replace = TRUE)
    row_idx <- unlist(cluster_idx[sampled], use.names = FALSE)

    ## Subset via integer indexing (fast for vectors and matrices)
    y_b    <- y_vec[row_idx]
    X_b    <- X_mat[row_idx, , drop = FALSE]
    cell_b <- cell_id[row_idx]

    i00 <- cell_b == 1L
    i01 <- cell_b == 2L
    i10 <- cell_b == 3L
    i11 <- cell_b == 4L

    if (min(sum(i00), sum(i01), sum(i10), sum(i11)) < 10) return(NULL)

    tryCatch(
      .boot_semi_fn(
        Y00 = y_b[i00], X00 = X_b[i00, , drop = FALSE],
        Y01 = y_b[i01], X01 = X_b[i01, , drop = FALSE],
        Y10 = y_b[i10], X10 = X_b[i10, , drop = FALSE],
        Y11 = y_b[i11], X11 = X_b[i11, , drop = FALSE],
        fixed_index = fixed_index,
        verbose = FALSE
      ),
      error = function(e) NULL
    )
  }

  ## Collect results
  boot_zeta <- matrix(NA, n_boot, J)
  boot_tau  <- matrix(NA, n_boot, 2)
  colnames(boot_tau) <- c("lower", "upper")
  n_fail <- 0

  for (b in seq_len(n_boot)) {
    if (!is.null(boot_results[[b]])) {
      boot_zeta[b, ] <- boot_results[[b]]$zeta
      boot_tau[b, ]  <- boot_results[[b]]$tau
    } else {
      n_fail <- n_fail + 1
    }
  }

  ## Compute CIs
  valid <- !is.na(boot_zeta[, 1])

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

  ## Relative effect
  re_se <- if (sum(valid) > 10) {
    apply(boot_tau[valid, , drop = FALSE], 2, stats::sd)
  } else {
    c(NA, NA)
  }

  relative_ci <- tryCatch(
    imbens_manski_ci(
      lb_hat = point$tau[1],
      ub_hat = point$tau[2],
      se_lb = re_se[1],
      se_ub = re_se[2],
      alpha = 0.05,
      truncate_diff = TRUE
    ),
    error = function(e) list(ci = c(NA, NA), c_crit = NA)
  )

  relative_res <- tibble::tibble(
    effect_lb = point$tau[1],
    effect_ub = point$tau[2],
    lower_ci  = relative_ci$ci[1],
    upper_ci  = relative_ci$ci[2],
    c_crit    = relative_ci$c_crit
  )

  list(
    inputs = list(
      outcome = outcome,
      post = post,
      treat = treat,
      cluster = cluster,
      covariates = covariates,
      n_boot = n_boot,
      method = method
    ),
    estimate_effects = diff_res,
    relative_effects = relative_res,
    n_boot_fail = n_fail
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

.GetPointEstimate <- function(df) {
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
    p0 = y1_cf$prob
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

.EstimateCounterfactualDist <- function(Y10, Y01, Y00) {
  # Estimate parameters for the observed outcome
  fit_obs <- .EstimateControlParams(
    Y00 = Y00,
    Y01 = Y01,
    Y10 = Y10
  )

  # Estimate parameters for the counterfactual outcome
  fit_cf <- .EstimateCounterfactualParams(fit = fit_obs)

  # Estimate the outcome proportions
  prob_y11 <- .EstimateOutcomeProportions(
    mu = fit_cf$mu,
    sd = fit_cf$sd,
    cutoff = fit_cf$cutoff
  )

  return(list(
    params = fit_cf,
    prob = prob_y11
  ))
}

.EstimateObservedDist <- function(Y11) {
  # Estimate parameters for the observed outcome
  fit_obs <- .EstimateLatentOutcomeParams(
    Ydt = Y11,
    cutoff = NULL
  )

  # Estimate the outcome proportions
  prob_y11 <- .EstimateOutcomeProportions(
    mu = fit_obs$mu,
    sd = fit_obs$sd,
    cutoff = fit_obs$cutoff
  )

  return(list(
    params = fit_obs,
    prob = prob_y11
  ))
}
