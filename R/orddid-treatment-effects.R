

.ComputeTreatmentEffects <- function(y1_prop, y0_prop) {
  effect_cumulative <- .CumulativeEffects(y1_prop, y0_prop)
  return(list(
    effect_cumulative = effect_cumulative,
    y1_prop = y1_prop,
    y0_prop = y0_prop
  ))
}

.CumulativeEffects <- function(y1_prop, y0_prop) {
  effect_est <- map2_dbl(y1_prop, y0_prop, function(y1, y0) {
    # Calculate cumulative effects
    cum_y1 <- cumsum(y1)
    cum_y0 <- cumsum(y0)

    # Calculate treatment effects
    return(cum_y1 - cum_y0)
  })

  # Compute the quantile based CI
  # effect_mat <- bind_rows(effect_est)
  # effect_ci <- apply(effect_mat, 2, function(x) {
  #   quantile(x, probs = c(0.025, 0.975))
  # })

  return(effect_est)
}
#' Compute sharp bounds on the relative treatment effect for ordinal outcomes
#'
#' Given the marginal distributions of the potential outcomes under treatment
#' (``p1``) and control (``p0``), this function returns the sharp upper and
#' lower bounds on the relative treatment effect γ as derived by Lu et al.
#' (2020) in Equations (3) and (4):contentReference[oaicite:0]{index=0}:contentReference[oaicite:1]{index=1}.
#' The bounds are expressed in closed form using auxiliary quantities δ⁠_j⁠m
#' and ξ⁠_j⁠m (Equations (2) and (5)) that depend only on the marginals
#':contentReference[oaicite:2]{index=2}:contentReference[oaicite:3]{index=3}.
#'
#' @param p1 A numeric vector of non‑negative probabilities giving the
#'   marginal distribution of the treated potential outcome ``Y(1)``. Its
#'   length ``J`` determines the number of ordinal categories (assumed to be
#'   labelled 0, …, J − 1). The elements should sum to one.
#' @param p0 A numeric vector of non‑negative probabilities giving the
#'   marginal distribution of the control potential outcome ``Y(0)``. It must
#'   be the same length as ``p1`` and should also sum to one.
#'
#' @return A list with components:
#'
#' * `lower` – the sharp lower bound γ_L on the relative treatment effect.
#' * `upper` – the sharp upper bound γ_U on the relative treatment effect.
#'
#' For binary outcomes (length of ``p1`` equal to 2), the bounds coincide with
#' the usual average treatment effect, i.e. the difference in means
#' ``sum((0:(J-1)) * p1) - sum((0:(J-1)) * p0)``:contentReference[oaicite:4]{index=4}.  If ``p1`` and
#' ``p0`` contain only a single category, the relative treatment effect is
#' identically zero.
#'
#' @examples
#' # Example with three ordinal categories
#' p1 <- c(0.2, 0.3, 0.5)
#' p0 <- c(0.5, 0.3, 0.2)
#' bounds <- compute_relative_effect_bounds(p1, p0)
#' bounds$lower  # lower bound
#' bounds$upper  # upper bound
#'
#' @export 
compute_relative_effect_bounds <- function(p1, p0) {
  # Check that p1 and p0 are numeric and of the same length
  if (!is.numeric(p1) || !is.numeric(p0)) {
    stop("p1 and p0 must be numeric vectors of probabilities.")
  }
  if (length(p1) != length(p0)) {
    stop("p1 and p0 must have the same length.")
  }
  J <- length(p1)
  # If there is only one category, the relative treatment effect is zero
  if (J == 1) {
    return(list(lower = 0, upper = 0))
  }
  # For a binary outcome, γ reduces to the difference in means
  if (J == 2) {
    # categories are labelled 0 and 1; the expectation of Y(1) is p1[2] and of
    # Y(0) is p0[2]
    gamma <- sum((0:(J-1)) * p1) - sum((0:(J-1)) * p0)
    return(list(lower = gamma, upper = gamma))
  }
  # Ensure the probabilities sum to one (within numerical tolerance).  A warning
  # is issued rather than a hard stop to allow slight floating point error.
  tol <- 1e-8
  if (abs(sum(p1) - 1) > tol) {
    warning("p1 does not sum to 1 (within tolerance)")
  }
  if (abs(sum(p0) - 1) > tol) {
    warning("p0 does not sum to 1 (within tolerance)")
  }
  # Compute δ_jm and ξ_jm for j = 1,…,J−1 and m = 1,…,J−j.  The formulae
  # correspond exactly to Equations (2) and (5) in Lu et al. (2020) after
  # translating the 0‑based indices to R's 1‑based indexing.
  delta_vals <- c()
  xi_vals <- c()
  for (j in 1:(J - 1)) {
    for (m in 1:(J - j)) {
      ## δ_jm (Equation (2))
      # Sum over k = j,…,J−1 (0‑based).  In R this becomes indices (j+1)…J.
      start1 <- j + 1
      sum1 <- if (start1 <= J) sum(p1[start1:J]) else 0
      # Sum over k = j + m,…,J−1 (0‑based).  Start index (j + m) + 1 = j + m + 1
      start2 <- j + m + 1
      sum2 <- if (start2 <= J) sum(p1[start2:J]) else 0
      # Sum over l = 0,…,j − 2 (0‑based).  Corresponds to R indices 1…(j−1).
      sum3 <- if (j >= 2) sum(p0[1:(j - 1)]) else 0
      # Sum over l = j + m − 1,…,J − 1 (0‑based).  Start index (j + m − 1) + 1 = j + m
      start4 <- j + m
      sum4 <- if (start4 <= J) sum(p0[start4:J]) else 0
      delta <- sum1 + sum2 + sum3 - sum4
      delta_vals <- c(delta_vals, delta)
      ## ξ_jm (Equation (5))
      # Sum over k = j + m − 1,…,J − 1 (0‑based).  Start index (j + m − 1) + 1 = j + m
      start_xi1 <- j + m
      sum_xi1 <- if (start_xi1 <= J) sum(p1[start_xi1:J]) else 0
      # Sum over l = j,…,J − 1 (0‑based).  Corresponds to R indices (j + 1)…J.
      start_xi2 <- j + 1
      sum_xi2 <- if (start_xi2 <= J) sum(p0[start_xi2:J]) else 0
      # Sum over l = j + m,…,J − 1 (0‑based).  Start index (j + m) + 1 = j + m + 1
      start_xi3 <- j + m + 1
      sum_xi3 <- if (start_xi3 <= J) sum(p0[start_xi3:J]) else 0
      # Sum over k = 0,…,j − 2 (0‑based).  R indices 1…(j − 1).
      sum_xi4 <- if (j >= 2) sum(p1[1:(j - 1)]) else 0
      xi <- sum_xi1 - sum_xi2 - sum_xi3 - sum_xi4
      xi_vals <- c(xi_vals, xi)
    }
  }
  # The sharp upper bound is the minimum δ_jm and the sharp lower bound is the
  # maximum ξ_jm:contentReference[oaicite:5]{index=5}:contentReference[oaicite:6]{index=6}.
  gamma_upper <- min(delta_vals)
  gamma_lower <- max(xi_vals)
  return(list(lower = gamma_lower, upper = gamma_upper))
}

#' Unit tests for `compute_relative_effect_bounds()`
#'
#' This helper function runs a suite of small tests to verify that
#' ``compute_relative_effect_bounds()`` accurately reproduces the analytical
#' formulas derived in Lu et al. (2020).  It checks the implementation
#' against hand‑computed expressions for the three‑category case (J = 3)
#':contentReference[oaicite:7]{index=7}, ensures the binary case behaves correctly and
#' exercises a few edge cases.  The function invisibly returns ``TRUE`` if all
#' tests pass and throws an error otherwise.
#'
#' @return ``TRUE`` if all checks succeed.  Errors are thrown on failure.
#' @examples
#' test_compute_bounds()  # run the tests
#' @noRd 
#
# test_compute_bounds <- function() {
#   # Helper to compute expected δ and ξ for J = 3 by explicit formulas.  The
#   # analytic expressions for δ⁠_11, δ⁠_12, δ⁠_21 and ξ⁠_11, ξ⁠_12, ξ⁠_21 are
#   # given in Remark 1 of Lu et al. (2020):contentReference[oaicite:8]{index=8}.
#   expected_bounds_J3 <- function(p1, p0) {
#     # p1 and p0 of length 3
#     # δ11 = p1[2] + 2*p1[3] − p0[2] − p0[3]
#     delta11 <- p1[2] + 2 * p1[3] - p0[2] - p0[3]
#     # δ12 = p1[2] + p1[3] − p0[3]
#     delta12 <- p1[2] + p1[3] - p0[3]
#     # δ21 = p1[3] − p0[3] + p0[1]
#     delta21 <- p1[3] - p0[3] + p0[1]
#     gammaU <- min(c(delta11, delta12, delta21))
#     # ξ11 = p1[2] + p1[3] − p0[2] − 2*p0[3]
#     xi11 <- p1[2] + p1[3] - p0[2] - 2 * p0[3]
#     # ξ12 = p1[3] − p0[2] − p0[3]
#     xi12 <- p1[3] - p0[2] - p0[3]
#     # ξ21 = p1[3] − p0[3] - p1[1]
#     xi21 <- p1[3] - p0[3] - p1[1]
#     gammaL <- max(c(xi11, xi12, xi21))
#     return(list(lower = gammaL, upper = gammaU))
#   }
#   # Test case: three categories with arbitrary probabilities
#   p1 <- c(0.2, 0.3, 0.5)
#   p0 <- c(0.5, 0.3, 0.2)
#   result <- compute_relative_effect_bounds(p1, p0)
#   expected <- expected_bounds_J3(p1, p0)
#   if (abs(result$lower - expected$lower) > 1e-8 ||
#       abs(result$upper - expected$upper) > 1e-8) {
#     stop("Test failed for J = 3 example 1.")
#   }
#   # Second three‑category example
#   p1 <- c(0.1, 0.2, 0.7)
#   p0 <- c(0.4, 0.4, 0.2)
#   result <- compute_relative_effect_bounds(p1, p0)
#   expected <- expected_bounds_J3(p1, p0)
#   if (abs(result$lower - expected$lower) > 1e-8 ||
#       abs(result$upper - expected$upper) > 1e-8) {
#     stop("Test failed for J = 3 example 2.")
#   }
#   # Binary outcome should give the difference in means
#   p1 <- c(0.3, 0.7)
#   p0 <- c(0.5, 0.5)
#   result <- compute_relative_effect_bounds(p1, p0)
#   gamma <- sum((0:1) * p1) - sum((0:1) * p0)
#   if (abs(result$lower - gamma) > 1e-8 || abs(result$upper - gamma) > 1e-8) {
#     stop("Test failed for binary outcome.")
#   }
#   # Edge case with degenerate distributions (all mass in one category)
#   p1 <- c(0, 1, 0)
#   p0 <- c(0, 0, 1)
#   result <- compute_relative_effect_bounds(p1, p0)
#   # For this configuration γ is necessarily ≤ 1 and ≥ −1
#   if (result$upper > 1 + 1e-8 || result$lower < -1 - 1e-8) {
#     stop("Test failed for degenerate distributions.")
#   }
#   # Random test: generate several random 4‑category distributions and check
#   # that the computed bounds satisfy gamma_lower ≤ gamma_upper and are within
#   # [−1, 1].
#   set.seed(123)
#   for (i in 1:5) {
#     x <- runif(4)
#     y <- runif(4)
#     p1 <- x / sum(x)
#     p0 <- y / sum(y)
#     result <- compute_relative_effect_bounds(p1, p0)
#     if (result$lower - 1e-8 > result$upper) {
#       stop("Random test failed: lower bound exceeds upper bound.")
#     }
#     if (result$upper > 1 + 1e-8 || result$lower < -1 - 1e-8) {
#       stop("Random test failed: bounds outside [−1,1].")
#     }
#   }
#   invisible(TRUE)
# }
