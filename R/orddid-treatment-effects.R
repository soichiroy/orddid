

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
