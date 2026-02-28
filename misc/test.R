## Example
data(gun_twowave_sub)


# library(tidyverse)

# # Format data
# Y00 <- gun_twowave_sub %>%
#   filter(treat10_12_2 == 0) %>%
#   pull(guns10)
# table(Y00)

# Y01 <- gun_twowave_sub %>%
#   filter(treat10_12_2 == 0) %>%
#   pull(guns12)
# table(Y01)

# Y10 <- gun_twowave_sub %>%
#   filter(treat10_12_2 == 1) %>%
#   pull(guns10)
# fit <- .EstimateControlParams(
#   Y00 = Y00,
#   Y01 = Y01,
#   Y10 = Y10
# )

# .EstimateCounterfactualParams(fit)

# prY_cf <- .EstimateOutcomeProportions(
#   mu = fit$y00$mu,
#   sd = fit$y00$sd,
#   cutoff = fit$y00$cutoff
# )

# Y11 <- gun_twowave_sub %>%
#   filter(treat10_12_2 == 1) %>%
#   pull(guns12)
# fit_y11 <- .EstimateLatentOutcomeParams(Ydt = Y11, cutoff = NULL)
# fit_y11
# prY11 <- .EstimateOutcomeProportions(
#   mu = fit_y11$mu,
#   sd = fit_y11$sd,
#   cutoff = fit$y00$cutoff
# )
# prY11 <- table(Y11)/length(Y11)

# .ComputeTreatmentEffects(
#   y1_prop = prY11,
#   y0_prop  = prY_cf
# )

# compute_relative_effect_bounds(
#   p1 = prY11,
#   p0 = prY_cf
# )