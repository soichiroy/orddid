## ============================================================
## orddid-semiparametric-profile.R
## Profile likelihood for the control-post and treated-pre cells
## ============================================================

#' Select a profile candidate, prioritizing successful optimization
#' @param fits list of profile optimization results
#' @return integer index of the selected fit
#' @keywords internal
.select_profile_fit <- function(fits) {
  if (!length(fits)) stop("No profile-likelihood candidates were supplied.")
  values <- vapply(fits, function(fit) fit$value, numeric(1))
  converged <- which(vapply(
    fits,
    function(fit) isTRUE(fit$convergence == 0L),
    logical(1)
  ))
  eligible <- if (length(converged)) converged else seq_along(fits)
  eligible[which.min(values[eligible])]
}

#' Profile likelihood estimator for one treatment-by-period cell
#'
#' Fits
#' `P(Y=j|X) = F((kappa[j+1]-X beta)/exp(X xi)) -
#' F((kappa[j]-X beta)/exp(X xi))`.
#'
#' Multiple deterministic starts are evaluated. The largest finite likelihood
#' among converged candidates is retained; nonconverged candidates are used
#' only when no candidate converged. With no analytic density, BFGS uses finite
#' differences. `optimizer = "both"` evaluates BFGS and Nelder-Mead from every
#' start, which is useful for regularized CDFs with finite jumps.
#'
#' @param Y integer outcomes coded 0, ..., J-1
#' @param X design matrix
#' @param Ftilde regularized CDF
#' @param kappa cutoff vector
#' @param beta_init initial location coefficients
#' @param xi_init initial log-scale coefficients
#' @param ftilde density for analytic BFGS gradients, or NULL
#' @param optimizer optimization method
#' @param multistart whether to use deterministic alternative starts
#' @return estimates, selected likelihood, and per-start diagnostics
#' @keywords internal
.step2_profile_likelihood <- function(
    Y, X, Ftilde, kappa, beta_init = NULL, xi_init = NULL,
    ftilde = NULL, optimizer = c("BFGS", "Nelder-Mead", "both"),
    multistart = TRUE) {
  optimizer <- match.arg(optimizer)
  Y <- as.integer(Y)
  X <- as.matrix(X)
  n <- length(Y)
  p <- ncol(X)
  J <- length(kappa) - 1L
  if (nrow(X) != n || any(!is.finite(X)) || anyNA(Y) ||
      any(Y < 0L | Y >= J) || length(unique(Y)) < 2L) {
    stop("Y, X, and kappa do not define a valid profile-likelihood cell.")
  }

  indicators <- matrix(0, n, J)
  for (j in seq_len(J)) indicators[, j] <- as.integer(Y == (j - 1L))

  if (is.null(beta_init)) {
    y0 <- as.integer(Y == 0L)
    beta_init <- tryCatch(
      -as.numeric(stats::coef(stats::glm(
        y0 ~ X - 1,
        family = stats::binomial(link = "probit")
      ))),
      error = function(e) rep(0, p)
    )
  }
  if (is.null(xi_init)) xi_init <- rep(0, p)
  if (length(beta_init) != p || length(xi_init) != p ||
      any(!is.finite(c(beta_init, xi_init)))) {
    stop("beta_init and xi_init must be finite vectors with ncol(X) entries.")
  }

  kappa_finite <- is.finite(kappa)
  cache <- new.env(parent = emptyenv())
  cache$theta <- NULL

  evaluate <- function(theta) {
    if (identical(theta, cache$theta)) return(invisible(NULL))
    beta <- theta[seq_len(p)]
    xi <- theta[p + seq_len(p)]
    mu <- as.numeric(X %*% beta)
    eta <- as.numeric(X %*% xi)
    if (any(!is.finite(c(mu, eta))) || any(abs(eta) > 30)) {
      cache$theta <- theta
      cache$value <- 1e100
      cache$gradient <- rep(0, 2L * p)
      return(invisible(NULL))
    }
    sigma <- exp(eta)
    has_gradient <- !is.null(ftilde) && optimizer != "Nelder-Mead"
    log_likelihood <- 0
    score_beta <- rep(0, p)
    score_xi <- rep(0, p)

    for (j in seq_len(J)) {
      if (kappa_finite[j]) {
        u_lower <- (kappa[j] - mu) / sigma
        F_lower <- Ftilde(u_lower)
        f_lower <- if (has_gradient) ftilde(u_lower) else NULL
      } else {
        u_lower <- rep(-Inf, n)
        F_lower <- rep(0, n)
        f_lower <- if (has_gradient) rep(0, n) else NULL
      }
      if (kappa_finite[j + 1L]) {
        u_upper <- (kappa[j + 1L] - mu) / sigma
        F_upper <- Ftilde(u_upper)
        f_upper <- if (has_gradient) ftilde(u_upper) else NULL
      } else {
        u_upper <- rep(Inf, n)
        F_upper <- rep(1, n)
        f_upper <- if (has_gradient) rep(0, n) else NULL
      }
      if (any(!is.finite(c(F_lower, F_upper)))) {
        log_likelihood <- -1e100
        break
      }

      probability <- pmax(F_upper - F_lower, 1e-12)
      log_likelihood <- log_likelihood +
        sum(indicators[, j] * log(probability))

      if (has_gradient) {
        weight <- indicators[, j] / probability
        score_beta <- score_beta + as.numeric(crossprod(
          X, weight * (f_lower - f_upper) / sigma
        ))
        f_u_lower <- ifelse(f_lower == 0, 0, f_lower * u_lower)
        f_u_upper <- ifelse(f_upper == 0, 0, f_upper * u_upper)
        score_xi <- score_xi + as.numeric(crossprod(
          X, weight * (f_u_lower - f_u_upper)
        ))
      }
    }

    cache$theta <- theta
    cache$value <- -log_likelihood
    cache$gradient <- if (has_gradient) {
      -c(score_beta, score_xi)
    } else {
      NULL
    }
  }
  objective <- function(theta) {
    evaluate(theta)
    cache$value
  }
  gradient <- if (!is.null(ftilde) && optimizer != "Nelder-Mead") {
    function(theta) {
      evaluate(theta)
      cache$gradient
    }
  } else {
    NULL
  }

  base_start <- c(beta_init, xi_init)
  starts <- list(base_start)
  if (isTRUE(multistart)) {
    cell_beta <- tryCatch(
      .ordered_probit_mle(Y, X)$beta,
      error = function(e) beta_init
    )
    if (length(cell_beta) != p || any(!is.finite(cell_beta))) {
      cell_beta <- beta_init
    }
    starts[[2L]] <- c(cell_beta, rep(0, p))

    location_direction <- tryCatch(
      solve(
        crossprod(X) + diag(1e-8, p),
        as.numeric(crossprod(X, rep(1, n)))
      ),
      error = function(e) rep(0, p)
    )
    direction_scale <- sqrt(mean(as.numeric(X %*% location_direction)^2))
    if (!is.finite(direction_scale) || direction_scale <= 1e-10) {
      location_direction <- rep(0, p)
    } else {
      location_direction <- location_direction / direction_scale
    }
    scale_shift <- log(1.5) * location_direction
    starts[[3L]] <- c(beta_init, xi_init + scale_shift)
    starts[[4L]] <- c(beta_init, xi_init - scale_shift)
  }

  start_matrix <- do.call(rbind, starts)
  start_key <- apply(round(start_matrix, 12L), 1L, paste, collapse = ":")
  starts <- starts[!duplicated(start_key)]
  run_optimizer <- function(start, start_id, method) {
    fit <- tryCatch(
      if (method == "BFGS") {
        stats::optim(
          start, objective, gr = gradient, method = "BFGS",
          control = list(maxit = 500L, reltol = 1e-8)
        )
      } else {
        stats::optim(
          start, objective, method = "Nelder-Mead",
          control = list(maxit = 1500L, reltol = 1e-8)
        )
      },
      error = function(e) NULL
    )
    if (is.null(fit) || !is.finite(fit$value) || fit$value >= 1e99 ||
        any(!is.finite(fit$par))) {
      return(NULL)
    }
    fit$start_id <- start_id
    fit$start <- start
    fit$optimizer <- method
    fit$analytic_gradient <- method == "BFGS" && !is.null(gradient)
    fit
  }

  methods <- if (optimizer == "both") {
    c("BFGS", "Nelder-Mead")
  } else {
    optimizer
  }
  fits <- unlist(lapply(seq_along(starts), function(start_id) {
    lapply(methods, function(method) {
      run_optimizer(starts[[start_id]], start_id, method)
    })
  }), recursive = FALSE)
  fits <- Filter(Negate(is.null), fits)
  if (!length(fits)) stop("All profile-likelihood optimization starts failed.")

  negative_log_likelihoods <- vapply(fits, function(fit) fit$value, numeric(1))
  selected <- fits[[.select_profile_fit(fits)]]
  candidate_fits <- data.frame(
    start_id = vapply(fits, function(fit) fit$start_id, integer(1)),
    loglik = -negative_log_likelihoods,
    convergence = vapply(fits, function(fit) fit$convergence, integer(1)),
    optimizer = vapply(fits, function(fit) fit$optimizer, character(1)),
    analytic_gradient = vapply(
      fits,
      function(fit) fit$analytic_gradient,
      logical(1)
    ),
    stringsAsFactors = FALSE
  )

  list(
    beta = selected$par[seq_len(p)],
    xi = selected$par[p + seq_len(p)],
    loglik = -selected$value,
    convergence = selected$convergence,
    message = selected$message,
    optimizer = selected$optimizer,
    analytic_gradient = selected$analytic_gradient,
    selected_start = selected$start_id,
    candidate_fits = candidate_fits,
    candidate_parameters = lapply(fits, function(fit) fit$par)
  )
}
