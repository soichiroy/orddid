ordDID_semiparametric <- function(
  Y00, X00,
  Y01, X01,
  Y10, X10,
  Y11, X11
) {
  # Data clearing
  # Fit Y00
  fit_Y00 <- ord_two_stage_fit(Y = Y00, X = X00)

  # Fit Y01
  fit_Y01 <- ord_locscale_fit(
    Y = Y01,
    X = X01,
    Fhat = fit_Y00$Fhat,
    alpha = fit_Y00$alpha
  )

  fit_Y10 <- ord_locscale_fit(
    Y = Y10,
    X = X10,
    Fhat = fit_Y00$Fhat,
    alpha = fit_Y00$alpha
  )

  # Get parameters for Y11
  mu10 <- X11 %*% fit_Y10$beta
  sd10 <- exp(X11 %*% fit_Y10$gamma)
  mu01 <- X11 %*% fit_Y01$beta
  sd01 <- exp(X11 %*% fit_Y01$gamma)
  mu00 <- X11 %*% fit_Y00$beta

  mu11 <- mu10 + (mu01 - mu00) * sd10
  sd11 <- sd01 * sd10
  kappa2 <- fit_Y00$alpha

  # Predicted probabilities
  F_S <- fit_Y00$Fhat$eval(mu11 / sd11)
  F_S_a <- fit_Y00$Fhat$eval((mu11 + kappa2) / sd11)
  prY11 <- cbind(
    p1 = F_S,
    p2 = pmax(0, F_S_a - F_S),
    p3 = pmax(0, 1 - F_S_a)
  )
  return(prY11)
}
