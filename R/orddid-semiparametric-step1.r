
# =========================================================
# One-call wrapper: Stage 1(ii) -> Stage 2 + reusable \hat F
# Returns β̂, α̂, \hat F, index g, and a predict_probs() closure
# =========================================================

#' Two-Stage Semiparametric Estimation for Ordered Response Models
#'
#' Estimates a 3-category ordered response model using the two-stage semiparametric
#' method from Yamauchi (2024). Stage 1 solves moment equations to estimate the
#' index coefficients (beta) and CDF via isotonic regression. Stage 2 estimates
#' the threshold parameter (alpha) using monotone zero-crossing.
#'
#' @param Y Numeric vector of ordered responses (must have exactly 3 categories).
#'   Will be recoded internally to 1, 2, 3 based on sorted unique values.
#' @param X Numeric matrix of covariates (n x K). Must include a constant/intercept
#'   column if desired (not added automatically).
#' @param fixed_index Integer indicating which coefficient to normalize to 1 for
#'   identification. Default is 1 (typically the intercept). This coefficient will
#'   be fixed at 1 in the estimation.
#' @param w Optional numeric vector of observation weights (length n). If NULL,
#'   all observations are weighted equally.
#' @param beta_init Optional numeric vector of initial values for beta (length K).
#'   If NULL, initialized via ordered probit MLE (with OLS fallback).
#' @param method Character string specifying the solver: "spectral" (default) uses
#'   spectral gradient with Barzilai-Borwein step and backtracking line search;
#'   "BB" uses BB::dfsane if the BB package is available.
#' @param verbose Logical. If TRUE, prints iteration progress and convergence info.
#'
#' @return A list with components:
#' \describe{
#'   \item{beta}{Numeric vector (length K) of estimated index coefficients with
#'     beta[fixed_index] = 1.}
#'   \item{alpha}{Numeric scalar for the estimated threshold parameter separating
#'     categories 2 and 3.}
#'   \item{Fhat}{FhatIsotonic object containing the estimated CDF. Has methods
#'     $predict() and $quantile() for evaluation.}
#'   \item{index}{Numeric vector (length n) of fitted index values X'beta.}
#'   \item{predict_probs}{Function closure taking new covariate matrix Xnew and
#'     returning predicted probabilities (n_new x 3 matrix with columns p1, p2, p3).}
#'   \item{stage1}{Full output list from stage1ii_solve including convergence info.}
#'   \item{stage2}{Full output list from EstimateCutoff including bracket and method.}
#' }
#'
#' @details
#' This function implements the complete two-stage estimation procedure:
#' \itemize{
#'   \item \strong{Stage 1(i):} Isotonic regression (PAVA) to estimate F(X'beta | beta)
#'   \item \strong{Stage 1(ii):} Solve moment condition E[X_{-fixed} * (1{Y=1} - F(X'beta))] = 0
#'   \item \strong{Stage 2:} Find alpha satisfying E[1 - 1{Y=3} - F(X'beta + alpha)] = 0
#' }
#'
#' The model assumes: P(Y = j | X) = F(X'beta + alpha_{j-1}) - F(X'beta + alpha_{j-2})
#' where alpha_0 = -Inf, alpha_1 = 0 (normalized), alpha_2 = alpha (estimated), alpha_3 = Inf.
#'
#' @references
#' Yamauchi, S. (2024). "Simple Semiparametric Estimation of Ordered Response Models."
#'
#' @examples
#' \dontrun{
#' # Generate 3-category ordered response data
#' set.seed(123)
#' n <- 500
#' X <- cbind(1, rnorm(n), rnorm(n))
#' beta_true <- c(1, 0.5, -0.3)
#' alpha_true <- 1.2
#' latent <- X %*% beta_true + rnorm(n)
#' Y <- cut(latent, breaks = c(-Inf, 0, alpha_true, Inf), labels = FALSE)
#'
#' # Estimate model
#' fit <- ord_two_stage_fit(Y, X, verbose = TRUE)
#'
#' # View results
#' print(fit$beta)
#' print(fit$alpha)
#' print(fit$stage1$converged)
#' print(fit$stage2$converged)
#'
#' # Predict probabilities for new data
#' X_new <- cbind(1, rnorm(10), rnorm(10))
#' probs <- fit$predict_probs(X_new)
#' }
#'
#' @seealso \code{\link{stage1ii_solve}}, \code{\link{EstimateCutoff}}, \code{\link{stage1i_isotonic}}
#' @export
ord_two_stage_fit <- function(
  Y,
  X,
  fixed_index = 1,
  w = NULL,
  beta_init = NULL,
  method = c("spectral", "BB"),
  verbose = TRUE
) {
  method <- match.arg(method)
  # Stage 1(ii)
  s1 <- stage1ii_solve(
    Y,
    X,
    fixed_index,
    w,
    beta_init,
    method,
    verbose = verbose
  )
  # Stage 2
  s2 <- EstimateCutoff(
    Y,
    X,
    beta_hat = s1$beta,
    Fhat = s1$Fhat,
    w = w,
    verbose = verbose
  )
  # Convenience
  predict_probs <- function(Xnew) {
    predict_probs_3cat(Xnew, s1$beta, s2$alpha, s1$Fhat)
  }
  list(
    beta = s1$beta,
    alpha = s2$alpha,
    Fhat = s1$Fhat,
    index = as.numeric(X %*% s1$beta),
    predict_probs = predict_probs,
    stage1 = s1,
    stage2 = s2
  )
}

# =========================================================
# Stage 1(i): isotonic PAVA for Fhat(·; beta)  (eq. (2.1))
# Returns a reusable object with $predict (CDF) and $quantile (left-inverse)
# =========================================================

#' Weighted Pool Adjacent Violators Algorithm (PAVA)
#' @keywords internal
.pava_weighted <- function(mu, wt) {
  len <- rep(1L, length(mu))
  i <- 1L
  while (i < length(mu)) {
    if (mu[i] > mu[i + 1] + 1e-15) {
      new_w <- wt[i] + wt[i + 1]
      new_mu <- (wt[i] * mu[i] + wt[i + 1] * mu[i + 1]) / new_w
      mu[i] <- new_mu
      wt[i] <- new_w
      len[i] <- len[i] + len[i + 1]
      mu <- mu[-(i + 1)]
      wt <- wt[-(i + 1)]
      len <- len[-(i + 1)]
      if (i > 1L) i <- i - 1L
    } else {
      i <- i + 1L
    }
  }
  list(mu = mu, wt = wt, len = len)
}

#' Compress tied observations
#' @keywords internal
.compress_ties <- function(u_sorted, y_sorted, w_sorted, tol = 1e-12) {
  grp_break <- c(TRUE, diff(u_sorted) > tol)
  gid <- cumsum(grp_break)
  list(
    u = as.numeric(tapply(u_sorted, gid, function(z) z[1])),
    sw = as.numeric(tapply(w_sorted, gid, sum)),
    sy = as.numeric(tapply(w_sorted * y_sorted, gid, sum))
  )
}

stage1i_isotonic <- function(Y, X, beta, w = NULL, tol = 1e-12) {
  n <- nrow(X)
  y1 <- as.integer(Y == min(Y))
  u <- as.numeric(X %*% beta)
  w <- if (is.null(w)) rep(1, n) else w

  # Sort by index
  ord <- order(u, method = "radix")
  uo <- u[ord]
  yo <- y1[ord]
  wo <- w[ord]

  # Compress ties and compute weighted means
  comp <- .compress_ties(uo, yo, wo, tol)
  mu <- comp$sy / comp$sw

  # PAVA isotonic regression
  pava <- .pava_weighted(mu, comp$sw)

  # Build step function
  block_end_idx <- cumsum(pava$len)
  knots <- comp$u[block_end_idx]
  values <- pava$mu

  # CDF evaluator
  Fhat_u <- function(ut) {
    j <- findInterval(ut, knots, rightmost.closed = TRUE)
    j[j == 0] <- 1L
    values[j]
  }

  # Quantile function (left-inverse)
  quantile_fun <- function(p, extended = FALSE) {
    vapply(
      p,
      function(pi) {
        if (!is.finite(pi) || pi < 0 || pi > 1) {
          return(NA_real_)
        }
        if (pi <= values[1]) {
          return(if (extended) -Inf else knots[1])
        }
        if (pi > values[length(values)]) {
          return(if (extended) Inf else knots[length(knots)])
        }
        knots[which(values >= pi)[1]]
      },
      numeric(1)
    )
  }

  # Predict function (flexible input)
  predict_fun <- function(Xnew, beta_for_X = beta) {
    if (is.numeric(Xnew) && is.null(dim(Xnew))) {
      Fhat_u(as.numeric(Xnew))
    } else {
      Fhat_u(as.numeric(as.matrix(Xnew) %*% beta_for_X))
    }
  }

  fitted <- numeric(n)
  fitted[ord] <- Fhat_u(uo)

  structure(
    list(
      knots = knots,
      values = values,
      predict = predict_fun,
      quantile = quantile_fun,
      export = function() list(knots = knots, values = values),
      u = u,
      beta = beta,
      fitted = fitted,
      eval = Fhat_u
    ),
    class = "FhatIsotonic"
  )
}

# Helpers to save/load Fhat across periods/scripts
save_Fhat <- function(Fhat, file) {
  saveRDS(Fhat$export(), file = file)
}

load_Fhat <- function(file_or_list, beta_for_X = NULL) {
  spec <- if (is.character(file_or_list)) {
    readRDS(file_or_list)
  } else {
    file_or_list
  }

  knots <- spec$knots
  values <- spec$values

  Fhat_u <- function(ut) {
    j <- findInterval(ut, knots, rightmost.closed = TRUE)
    j[j == 0] <- 1L
    values[j]
  }

  predict_fun <- function(Xnew, beta_for_X = beta_for_X) {
    if (is.numeric(Xnew) && is.null(dim(Xnew))) {
      Fhat_u(as.numeric(Xnew))
    } else {
      Fhat_u(as.numeric(as.matrix(Xnew) %*% beta_for_X))
    }
  }

  quantile_fun <- function(p, extended = FALSE) {
    vapply(
      p,
      function(pi) {
        if (!is.finite(pi) || pi < 0 || pi > 1) {
          return(NA_real_)
        }
        if (pi <= values[1]) {
          return(if (extended) -Inf else knots[1])
        }
        if (pi > values[length(values)]) {
          return(if (extended) Inf else knots[length(knots)])
        }
        knots[which(values >= pi)[1]]
      },
      numeric(1)
    )
  }

  structure(
    list(
      knots = knots,
      values = values,
      predict = predict_fun,
      quantile = quantile_fun,
      export = function() list(knots = knots, values = values),
      beta = beta_for_X
    ),
    class = "FhatIsotonic"
  )
}

# =========================================================
# Stage 1(ii): solve Ψ_n(β)=0 (eq. (2.2))
# Spectral method with backtracking or BB::dfsane
# =========================================================

#' Initialize beta via probit regression
#' @keywords internal
.init_beta_probit <- function(Y, X, fixed_index) {
  y1 <- as.integer(Y == min(Y))
  K <- ncol(X)

  # Try probit GLM
  dfX <- as.data.frame(X)
  names(dfX) <- paste0("x", seq_len(K))
  form <- as.formula(paste("y1 ~", paste(names(dfX), collapse = " + "), "- 1"))

  bfull <- try(
    coef(glm(
      form,
      data = cbind(dfX, y1 = y1),
      family = binomial(link = "probit")
    )),
    silent = TRUE
  )

  # Fallback to OLS if probit fails
  if (inherits(bfull, "try-error") || !is.finite(bfull[fixed_index])) {
    bfull <- try(solve(crossprod(X), crossprod(X, y1)), silent = TRUE)
  }

  if (inherits(bfull, "try-error") || !is.finite(bfull[fixed_index])) {
    stop("Failed to construct beta_init; please provide beta_init.")
  }

  as.numeric(bfull / bfull[fixed_index])
}

#' Spectral method (Barzilai-Borwein with backtracking)
#' @keywords internal
.solve_spectral <- function(
  moment_map,
  b_init,
  maxit,
  tol,
  bt_factor,
  rel_improve,
  verbose
) {
  b <- b_init
  g <- moment_map(b)
  gnorm <- max(abs(g))

  if (verbose) {
    message(sprintf("iter %3d |score|_inf = %.3e", 0L, gnorm))
  }

  prev_b <- b
  prev_g <- g
  step <- 1
  iter <- 0L

  while (gnorm > tol && iter < maxit) {
    # Compute spectral step size
    if (iter > 0) {
      s <- b - prev_b
      y <- g - prev_g
      denom <- sum(s * y)
      step <- if (is.finite(denom) && abs(denom) > 1e-12) {
        sum(s * s) / denom
      } else {
        1
      }
      if (!is.finite(step) || step <= 0) step <- 1
    }

    # Line search with backtracking
    cand_b <- b - step * g
    new_g <- moment_map(cand_b)
    new_norm <- max(abs(new_g))
    bt <- 0L

    while (new_norm > rel_improve * gnorm && bt < 20L) {
      step <- step * bt_factor
      cand_b <- b - step * g
      new_g <- moment_map(cand_b)
      new_norm <- max(abs(new_g))
      bt <- bt + 1L
    }

    prev_b <- b
    prev_g <- g
    b <- cand_b
    g <- new_g
    gnorm <- new_norm
    iter <- iter + 1L

    if (verbose) {
      message(sprintf(
        "iter %3d |score|_inf = %.3e (step=%.2e, bt=%d)",
        iter,
        gnorm,
        step,
        bt
      ))
    }
  }

  list(
    par = b,
    converged = (gnorm <= tol),
    iter = iter,
    score = g,
    score_norm = gnorm
  )
}

#' BB::dfsane solver wrapper
#' @keywords internal
.solve_BB <- function(moment_map, b_init, maxit, tol, verbose) {
  root <- BB::dfsane(
    par = b_init,
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


# stage1ii_solve
stage1ii_solve <- function(
  Y,
  X,
  fixed_index = 1,
  w = NULL,
  beta_init = NULL,
  method = c("spectral", "BB"),
  maxit = 200,
  tol = 1e-6,
  bt_factor = 0.5,
  rel_improve = 0.9,
  verbose = TRUE
) {
  # enforce normalization
  .compose_beta <- function(b) {
    bf <- numeric(K)
    bf[fixed_index] <- 1
    bf[-fixed_index] <- b
    bf
  }

  # Moment condition: E[X_{-fixed} * (1{Y=1} - Fhat(X'beta))] = 0
  .moment_map <- function(b) {
    beta_full <- .compose_beta(b)
    u <- as.numeric(X %*% beta_full)
    Fhat <- stage1i_isotonic(Y, X, beta_full, w = ww)
    r <- y1 - Fhat$predict(u)
    as.numeric(colSums(X_minus * (r * ww)) / sum(ww))
  }

  # Organize inputs
  method <- match.arg(method)
  X <- as.matrix(X)
  n <- nrow(X)
  K <- ncol(X)

  y1 <- as.integer(Y == min(Y))
  X_minus <- X[, -fixed_index, drop = FALSE]
  ww <- if (is.null(w)) rep(1, n) else w

  # Initialize via probit if not provided
  if (is.null(beta_init)) {
    beta_init <- .init_beta_probit(Y, X, fixed_index)
  }
  b <- beta_init[-fixed_index]

  # Solve moment equations using selected method
  if (method == "BB" && requireNamespace("BB", quietly = TRUE)) {
    if (verbose) {
      message("Using BB::dfsane ...")
    }
    result <- .solve_BB(.moment_map, b, maxit, tol, verbose)
  } else {
    if (verbose) {
      message(
        "Using spectral residual (Barzilai–Borwein) with backtracking ..."
      )
    }
    result <- .solve_spectral(
      .moment_map,
      b,
      maxit,
      tol,
      bt_factor,
      rel_improve,
      verbose
    )
  }

  # Extract results
  bhat <- result$par
  converged <- result$converged
  iter <- result$iter
  score <- result$score
  score_norm <- result$score_norm

  # Compose the full beta vector
  beta_hat <- .compose_beta(bhat)

  # get the final estimate of F given the estimated beta
  Fhat <- stage1i_isotonic(Y, X, beta_hat, w = ww)

  list(
    beta = beta_hat,
    converged = converged,
    iter = iter,
    score = score,
    score_norm = score_norm,
    Fhat = Fhat,
    method = method
  )
}

# =========================================================
# Stage 2: threshold α via monotone zero-crossing (eq. (2.3))
# =========================================================
EstimateCutoff <- function(
  Y,
  X,
  beta_hat,
  Fhat = NULL,
  w = NULL,
  alpha_lower = 0,
  alpha_upper = NULL,
  tol = 1e-6,
  max_expand = 40L,
  grid_points = 2000L,
  verbose = TRUE
) {
  X <- as.matrix(X)
  y3 <- as.integer(Y == max(Y))
  u <- as.numeric(X %*% beta_hat)
  ww <- if (is.null(w)) rep(1, nrow(X)) else w
  denom <- sum(ww)

  # Moment condition for alpha
  psi <- function(a) sum(ww * (1 - y3 - Fhat$predict(u + a))) / denom

  # Initialize bracket
  lo <- max(alpha_lower, 0)
  psi_lo <- psi(lo)

  if (psi_lo <= 0) {
    if (verbose) {
      message(sprintf("Stage 2: Ψ(lo)=%.3e ≤ 0; α̂ = %.6g.", psi_lo, lo))
    }
    return(list(
      alpha = lo,
      psi_at = psi_lo,
      converged = TRUE,
      bracket = c(lo, lo),
      method = "boundary",
      Fhat = Fhat
    ))
  }

  # Find upper bracket
  if (is.null(alpha_upper)) {
    span <- diff(range(u))
    hi <- if (is.finite(span) && span > 0) 2 * span else 2 * max(sd(u), 1)
  } else {
    hi <- alpha_upper
  }

  psi_hi <- psi(hi)
  n_expand <- 0L
  while (psi_hi > 0 && n_expand < max_expand) {
    hi <- hi * 2
    psi_hi <- psi(hi)
    n_expand <- n_expand + 1L
  }

  # Solve for zero
  if (psi_hi <= 0) {
    root <- try(uniroot(psi, lower = lo, upper = hi, tol = tol), silent = TRUE)
    if (!inherits(root, "try-error")) {
      return(list(
        alpha = root$root,
        psi_at = psi(root$root),
        converged = TRUE,
        bracket = c(lo, hi),
        method = "uniroot",
        Fhat = Fhat
      ))
    }
  }

  # Fallback: grid search
  grid <- seq(lo, hi, length.out = grid_points)
  vals <- vapply(grid, psi, numeric(1))
  idx <- which(vals <= 0)[1]

  if (is.na(idx)) {
    alpha <- grid[which.min(abs(vals))]
    method <- "grid-min"
    crossed <- FALSE
  } else {
    alpha <- grid[idx]
    method <- "grid-cross"
    crossed <- TRUE
  }

  list(
    alpha = alpha,
    psi_at = psi(alpha),
    converged = crossed,
    bracket = c(lo, hi),
    method = method,
    Fhat = Fhat
  )
}

# =========================================================
# Predicted probabilities (3-category ordered model; eq. (1.2))
# =========================================================
predict_probs_3cat <- function(X, beta_hat, alpha_hat, Fhat) {
  g <- as.numeric(X %*% beta_hat)
  Fg <- Fhat$predict(g)
  Fga <- Fhat$predict(g + alpha_hat)
  p1 <- Fg
  p2 <- pmax(0, Fga - Fg)
  p3 <- pmax(0, 1 - Fga)
  cbind(p1 = p1, p2 = p2, p3 = p3)
}
