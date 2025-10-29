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
#' fit <- ord_did(
#'   Ynew = gun_twowave %>% filter(year == 2012) %>% pull(guns),
#'   Yold = gun_twowave %>% filter(year == 2010) %>% pull(guns),
#'   treat = gun_twowave %>% filter(year == 2012) %>% pull(treat_100mi),
#'   id_cluster = gun_twowave %>% filter(year == 2010) %>% pull(reszip),
#'   n_boot = 10,
#'   pre = FALSE,
#'   verbose = FALSE
#' )
#'
#' ## view summary of the output
#' ## non-cumulative effects
#' summary(fit, cumulative = FALSE)
#'
#' ## cumulative effects
#' summary(fit)
#' }
#' @export
ord_did <- function(data, outcome, post, treat, cluster) {
  df_format <- .orddid_extract_Y_cells(
    data = data,
    outcome = outcome,
    post = post,
    treat = treat
  )

  # Compute point estimate
  point_estimate <- .GetPointEstimate(df = df_format)

  # Bootstrap
  boot_estimate <- lapply(1:50, function(i) {
    .RunBootstrap(
      df = data,
      cluster = cluster,
      outcome = outcome,
      post = post,
      treat = treat
    )
  })

  ci_diff <- .ComputeCI(boot_estimates = boot_estimate, alpha = 0.05)
  # Summarize results
  res <- tibble(
    category = 1:length(point_estimate$estimated_effects$diff_effects),
    effect = point_estimate$estimated_effects$diff_effects,
    lower_ci = ci_diff$lower,
    upper_ci = ci_diff$upper
  )
  return(list(
    estimate_effects = res
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

# ord_did <- function(Ynew, Yold, treat, id_cluster = NULL, cut = c(0, 1),
#                     n_boot = 500, pre = FALSE, verbose = FALSE) {

#   ## quick input check
#   ord_did_check_input(Ynew, Yold, treat, id_cluster, n_boot)

#   ## data summary
#   Yc <- c(Ynew, Yold)
#   J  <- length(unique(Yc))
#   n  <- length(treat)
#   n1 <- sum(treat)

#   ## fit on the obs data
#   fit <- ord_did_run(Ynew, Yold, treat, cut, pre)

#   ## bootstrap (assuming a panel)
#   boot <- ord_did_boot(Ynew, Yold, treat, cut, id_cluster, J, n_boot, verbose, pre)

#   ## return objects
#   return_list <- list(
#     'fit' = fit,
#     'boot' = boot$boot_save,
#     'boot_params' = boot$boot_params)

#   ## add attributes to the returning object
#   attr(return_list, "input")    <- list(n_boot = n_boot, pre = pre, cut = cut)
#   attr(return_list, "n_choice") <- J
#   attr(return_list, "n")        <- n
#   attr(return_list, "n1")       <- n1

#   class(return_list) <- c("orddid", "orddid.fit")
#   return(return_list)
# }
