#' Ordinal Difference-in-Differences for Panel Data
#'
#' \code{ord_did()} implements the difference-in-differences for the ordinal outcome.
#'
#' @param Ynew A numeric vector of ordinal outcome for the post-treatment period.
#' @param Yold A numeric vector of ordinal outcome for the pre-treatment period.
#' @param treat A numeric vector of treatment indicator.
#'  The treatment group should take 1 and the control group should take 0.
#' @param id_cluster A vector of cluster id.
#'   If left as \code{NULL}, bootstrap is implemented at the individual level.
#' @param cut A vector of cutoffs. Two numeric values should be specified. Default is \code{cut = c(0, 1)}.
#' @param n_boot The number of boostrapt iterations for estimating the variance. Default is \code{n_boot = 500}.
#' @param pre A boolean argument used to indicate if the data comes entirely from pre-treatment periods.
#'  This should be \code{TRUE} when the output is supplied to \code{\link{equivalence_test}}.
#' @param verbose If \code{TRUE}, print the progress of bootstrap iterations.
#' @return \code{ord_did()} returns a list of class `orddid' containing the following components:
#' \item{fit}{A list with the output of the ordinal DID estimators,
#'            which contains parameter estimates and predicted probabilities for each category.}
#' \item{boot}{A list with the output of bootstraps,
#'             which contains parameter estimates and predicted probabilities for each category.}
#' \item{boot_params}{A list with all objects generated during the bootstrap step.}
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
  use_parallel = FALSE
) {
  df_format <- .orddid_extract_Y_cells(
    data = data,
    outcome = outcome,
    post = post,
    treat = treat
  )

  # Compute point estimate
  point_estimate <- .GetPointEstimate(df = df_format)

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
      n_boot = n_boot
    ),
    estimate_effects = diff_res,
    relative_effects = relative_res
  ))
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
