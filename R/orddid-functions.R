#' Estimate parameters for the observed outcome
#' 
#' @noRd
#' @keywords internal
.EstimateControlParams <- function(Y00, Y01, Y10) {
  # This ordering of estimation assumes that we fix sd_00 = 1.
  # Y00: Estimate mu_00, and the cutoffs
  fit_y00 <- .EstimateLatentOutcomeParams(Ydt = Y00, cutoff = NULL)

  # Y01: Estimate mu_01 and sd_01 with cutoff fixed
  fit_y01 <- .EstimateLatentOutcomeParams(Ydt = Y01, cutoff = fit_y00$cutoff)

  # Y10: Estimate mu_10 and sd_10 with cutoff fixed
  fit_y10 <- .EstimateLatentOutcomeParams(Ydt = Y10, cutoff = fit_y00$cutoff)

  return(list(y00 = fit_y00, y01 = fit_y01, y10 = fit_y10))
}

#' Estimate parameters for the latent outcome (mu_dt, sd_dt)
#'
#' @param Y00 A vector of pre-treatment control observations.
#' @param cutoff A value for the first cutoff value to be fixed. Default is 0.
#'
#' @noRd
.EstimateLatentOutcomeParams <- function(Ydt, cutoff) {
  ## Implement the ordinal probit regression only with the intercept.
  fit <- .FitOrderedProbit(Y = Ydt, cutoff_val = cutoff)
  return(list(mu = fit$mu, sd = fit$sd, cutoff = fit$cutoff))
}

#' Estimate parameters for the counterfactual outcome.
#'
#' Use the identification formula to estimate mu_11 and sd_11
#' @noRd
.EstimateCounterfactualParams <- function(fit) {
  mu11 <- with(fit, y10$mu + (y01$mu - y00$mu) * (y10$sd / y00$sd))
  sd11 <- with(fit, y10$sd * y01$sd / y00$sd)

  return(list(mu = mu11, sd = sd11, cutoff = fit$y00$cutoff))
}

.EstimateOutcomeProportions <- function(mu, sd, cutoff) {
  # calculate the cutoff values
  cut_new <- .PadInf(cutoff)
  prob_yj <- rep(NA, (length(cutoff) + 1))

  for (j in 1:(length(cutoff) + 1)) {
    prob_yj[j] <- pnorm(cut_new[j + 1], mean = mu, sd = sd) -
      pnorm(cut_new[j], mean = mu, sd = sd)
  }

  return(prob_yj)
}

.PadInf <- function(x) {
  return(c(-Inf, x, Inf))
}



.ComputeTreatmentEffects <- function(y1_prop, y0_prop) {
  effect_cumulative <- .CumulativeEffects(y1_prop, y0_prop)
  return(list(
    effect_cumulative = effect_cumulative,
    y1_prop = y1_prop,
    y0_prop = y0_prop
  ))
}

.CumulativeEffects <- function(y1_prop, y0_prop) {
  effect_est <- map2(y1_prop, y0_prop, function(y1, y0) {
    # Calculate cumulative effects
    cum_y1 <- cumsum(y1)
    cum_y0 <- cumsum(y0)

    # Calculate treatment effects
    return(cum_y1 - cum_y0)
  })

  # Compute the quantile based CI
  effect_mat <- bind_rows(effect_est)
  effect_ci <- apply(effect_mat, 2, function(x) {
    quantile(x, probs = c(0.025, 0.975))
  })

  return(effect_ci)
}

.RelativeEffect <- function() {}

