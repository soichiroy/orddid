## ============================================================
## orddid-semiparametric-profile.R
## Profile likelihood for Step 2: estimate (beta_dt, xi_dt)
##
## Convention:
##   P(Y = j | X) = Ftilde((kappa_{j+1} - X'beta) / exp(X'xi))
##                - Ftilde((kappa_j   - X'beta) / exp(X'xi))
## ============================================================

#' Profile likelihood estimator for a single (d,t) cell
#'
#' @param Y       integer vector of outcomes (0, 1, ..., J-1)
#' @param X       n x p design matrix (with intercept column)
#' @param Ftilde  smoothed CDF function from .smooth_Fhat()$Ftilde
#' @param kappa   cutoff vector c(-Inf, 0, alpha, Inf) from Step 1
#' @param beta_init optional initial value for beta (length p)
#' @param xi_init   optional initial value for xi (length p)
#' @param ftilde  smoothed density function (if provided, analytical gradient
#'   is used for BFGS, avoiding finite-difference approximation)
#' @return list with beta, xi, loglik, convergence
#' @keywords internal
.step2_profile_likelihood <- function(Y, X, Ftilde, kappa,
                                      beta_init = NULL,
                                      xi_init = NULL,
                                      ftilde = NULL) {

  n <- length(Y)
  p <- ncol(X)
  J <- length(kappa) - 1

  ind <- matrix(0, n, J)
  for (j in seq_len(J)) {
    ind[, j] <- as.integer(Y == (j - 1))
  }

  if (is.null(beta_init)) {
    y1 <- as.integer(Y == 0)
    beta_init <- tryCatch({
      c_hat <- stats::coef(stats::glm(y1 ~ X - 1,
                                       family = stats::binomial(link = "probit")))
      -as.numeric(c_hat)
    }, error = function(e) {
      rep(0, p)
    })
  }

  if (is.null(xi_init)) {
    xi_init <- rep(0, p)
  }

  theta_init <- c(beta_init, xi_init)

  ## Precompute which kappa values are finite (avoids per-call ifelse)
  kappa_finite <- is.finite(kappa)

  ## Cache for combined fn+gradient evaluation (BFGS calls fn and gr
  ## at the same theta, so we avoid computing twice)
  .cache <- new.env(parent = emptyenv())
  .cache$theta <- NULL
  .cache$val   <- NULL
  .cache$grad  <- NULL

  .eval_at <- function(theta) {
    if (identical(theta, .cache$theta)) return(invisible(NULL))

    beta  <- theta[1:p]
    xi    <- theta[(p + 1):(2 * p)]
    mu    <- as.numeric(X %*% beta)
    sigma <- exp(as.numeric(X %*% xi))

    ll <- 0
    dll_dbeta <- rep(0, p)
    dll_dxi   <- rep(0, p)
    has_grad  <- !is.null(ftilde)

    for (j in seq_len(J)) {
      ## Lower bound
      if (kappa_finite[j]) {
        u_lo <- (kappa[j] - mu) / sigma
        F_lo <- Ftilde(u_lo)
        f_lo <- if (has_grad) ftilde(u_lo) else NULL
      } else {
        u_lo <- rep(if (kappa[j] < 0) -Inf else Inf, n)
        F_lo <- rep(if (kappa[j] < 0) 0 else 1, n)
        f_lo <- if (has_grad) rep(0, n) else NULL
      }

      ## Upper bound
      if (kappa_finite[j + 1]) {
        u_hi <- (kappa[j + 1] - mu) / sigma
        F_hi <- Ftilde(u_hi)
        f_hi <- if (has_grad) ftilde(u_hi) else NULL
      } else {
        u_hi <- rep(if (kappa[j + 1] < 0) -Inf else Inf, n)
        F_hi <- rep(if (kappa[j + 1] < 0) 0 else 1, n)
        f_hi <- if (has_grad) rep(0, n) else NULL
      }

      prob_j <- pmax(F_hi - F_lo, 1e-12)
      ll <- ll + sum(ind[, j] * log(prob_j))

      if (has_grad) {
        w_j <- ind[, j] / prob_j

        ## d P_j / d beta_k = (f_lo - f_hi) * x_k / sigma
        dll_dbeta <- dll_dbeta +
          as.numeric(crossprod(X, w_j * (f_lo - f_hi) / sigma))

        ## d P_j / d xi_k = (f_lo*u_lo - f_hi*u_hi) * x_k
        ## Handle 0 * Inf = 0 for boundary terms
        fu_lo <- ifelse(f_lo == 0, 0, f_lo * u_lo)
        fu_hi <- ifelse(f_hi == 0, 0, f_hi * u_hi)
        dll_dxi <- dll_dxi +
          as.numeric(crossprod(X, w_j * (fu_lo - fu_hi)))
      }
    }

    .cache$theta <- theta
    .cache$val   <- -ll
    .cache$grad  <- if (has_grad) -c(dll_dbeta, dll_dxi) else NULL
  }

  negloglik <- function(theta) {
    .eval_at(theta)
    .cache$val
  }

  grad_negloglik <- if (!is.null(ftilde)) {
    function(theta) {
      .eval_at(theta)
      .cache$grad
    }
  } else {
    NULL
  }

  fit <- stats::optim(theta_init, negloglik,
                      gr = grad_negloglik,
                      method = "BFGS",
                      control = list(maxit = 500, reltol = 1e-8))

  list(
    beta        = fit$par[1:p],
    xi          = fit$par[(p + 1):(2 * p)],
    loglik      = -fit$value,
    convergence = fit$convergence,
    message     = fit$message
  )
}
