# =========================================================
# Location-scale model: simultaneous mean/scale estimation
# Yamauchi (2024) Section 4.1, Appendix A.6
# =========================================================


ord_locscale_fit <- function(
  Y,
  X,
  Fhat,
  alpha,
  fixed_index_mu = 1,
  w = NULL,
  beta_init = NULL,
  gamma_init = NULL,
  method = c("spectral", "BB"),
  maxit = 300,
  tol = 1e-6,
  bt_factor = 0.5,
  rel_improve = 0.9,
  verbose = TRUE
) {
  method <- match.arg(method)
  X_mu <- as.matrix(X)
  W_sig <- as.matrix(X)
  n <- nrow(X_mu)
  p_mu <- ncol(X_mu)
  p_gm <- ncol(W_sig)

  y1 <- as.integer(Y == min(Y))
  y3 <- as.integer(Y == max(Y))
  ww <- if (is.null(w)) rep(1, n) else w
  sw <- sum(ww)


  # Moment conditions
  moment_map <- function(theta) {
    pa <- unpack(theta)
    beta <- pa$b
    mu <- as.numeric(X_mu %*% beta)
    sigma <- as.numeric(exp(W_sig %*% pa$gamma))
    S <- mu / sigma

    g1 <- y1 - Fhat$predict(S)
    g2 <- (1 - y3) - Fhat$predict(S + alpha / sigma)

    m1 <- as.numeric(colSums(X_mu * (g1 * ww)) / sw)
    m2 <- as.numeric(colSums(W_sig * (g2 * ww)) / sw)

    c(m1, m2)
  }

  # Initialize parameters
  if (is.null(beta_init)) {
    beta_init <- .init_beta_probit(Y, X_mu, fixed_index_mu)
  }
  if (is.null(gamma_init)) {
    gamma_init <- rep(0, p_gm)
  }

  theta <- c(beta_init, gamma_init)

  unpack <- function(theta) {
    list(
      b = theta[seq_len(p_mu)],
      gamma = theta[(p_mu + 1):length(theta)]
    )
  }

  # Solve moment equations using selected method
  if (method == "BB" && requireNamespace("BB", quietly = TRUE)) {
    if (verbose) {
      message("Using BB::dfsane …")
    }
    result <- .solve_BB_locscale(moment_map, theta, maxit, tol, verbose)
  } else {
    if (verbose) {
      message("Using spectral residual (Barzilai–Borwein) with backtracking …")
    }
    result <- .solve_spectral(
      moment_map,
      theta,
      maxit,
      tol,
      bt_factor,
      rel_improve,
      verbose
    )
  }

  # Extract results
  theta <- result$par
  conv <- result$converged
  iter <- result$iter
  score <- result$score
  score_norm <- result$score_norm

  # Extract parameter estimates
  pa <- unpack(theta)
  beta_hat <- pa$b
  gamma_hat <- pa$gamma
  mu_hat <- as.numeric(X_mu %*% beta_hat)
  sigma_hat <- as.numeric(exp(W_sig %*% gamma_hat))
  S_hat <- mu_hat / sigma_hat

  # Prediction function
  predict_probs <- function(Xnew_mu, Wnew_sig) {
    Xnew_mu <- as.matrix(Xnew_mu)
    Wnew_sig <- as.matrix(Wnew_sig)
    mu <- as.numeric(Xnew_mu %*% beta_hat)
    sig <- as.numeric(exp(Wnew_sig %*% gamma_hat))
    S <- mu / sig
    F_S <- Fhat$predict(S)
    F_S_a <- Fhat$predict(S + alpha / sig)
    cbind(
      p1 = F_S,
      p2 = pmax(0, F_S_a - F_S),
      p3 = pmax(0, 1 - F_S_a)
    )
  }

  list(
    beta = beta_hat,
    gamma = gamma_hat,
    mu = mu_hat,
    sigma = sigma_hat,
    S = S_hat,
    converged = conv,
    iter = iter,
    score = score,
    score_norm = score_norm,
    predict_probs = predict_probs
  )
}

# Convenience wrapper for predictions
predict_probs_locscale <- function(fit, X_mu, W_sig, Fhat, alpha) {
  fit$predict_probs(X_mu, W_sig)
}


#' BB::dfsane solver wrapper for location-scale model
#' @keywords internal
.solve_BB_locscale <- function(moment_map, theta_init, maxit, tol, verbose) {
  root <- BB::dfsane(
    par = theta_init,
    fn = moment_map,
    control = list(maxit = maxit, tol = tol, trace = if (verbose) 1 else 0)
  )

  list(
    par = root$par,
    converged = (root$convergence == 0),
    iter = root$feval,
    score = moment_map(root$par),
    score_norm = max(abs(moment_map(root$par)))
  )
}
