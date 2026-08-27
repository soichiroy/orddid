## ============================================================
## orddid-semiparametric-first-stage.R
## Paper-aligned first stage for the three-category model
## ============================================================

#' Initialize the semiparametric first stage by ordered probit
#' @keywords internal
.initialize_semiparametric_j3 <- function(Y, X, fixed_index = 1L,
                                          scale_tol = 1e-8) {
  Y <- as.integer(Y)
  X <- as.matrix(X)
  if (anyNA(Y) || !identical(sort(unique(Y)), 0:2)) {
    stop("Y00 must contain exactly the three categories 0, 1, and 2.")
  }
  if (nrow(X) != length(Y) || ncol(X) < 1L || any(!is.finite(X))) {
    stop("X00 must be a finite numeric matrix with one row per outcome.")
  }
  if (length(fixed_index) != 1L || fixed_index < 1L ||
      fixed_index > ncol(X)) {
    stop("fixed_index must identify a column of X00.")
  }

  fit <- .ordered_probit_mle(Y, X)
  if (length(fit$beta) != ncol(X) || any(!is.finite(fit$beta)) ||
      length(fit$kappa) != 4L || !is.finite(fit$alpha) || fit$alpha <= 0) {
    stop("Ordered-probit initialization did not return finite J=3 parameters.")
  }

  b_raw <- -as.numeric(fit$beta)
  scale <- b_raw[fixed_index]
  if (!is.finite(scale) || scale <= scale_tol) {
    stop(
      "The normalized covariate has the wrong or near-zero orientation: ",
      "the ordered-probit estimate requires -beta[fixed_index] > 0."
    )
  }

  b <- b_raw / scale
  b[fixed_index] <- 1
  alpha <- fit$alpha / scale
  if (!is.finite(alpha) || alpha <= 0) {
    stop("The normalized ordered-probit cutoff is not positive and finite.")
  }

  list(
    b = b,
    beta = -b,
    alpha = alpha,
    scale = scale,
    ordered_probit = fit
  )
}

#' Identify columns whose location cannot be updated with an unrestricted CDF
#' @keywords internal
.semiparametric_constant_columns <- function(X, tol = sqrt(.Machine$double.eps)) {
  X <- as.matrix(X)
  spans <- apply(X, 2L, function(x) diff(range(x)))
  scales <- apply(X, 2L, function(x) max(1, max(abs(x))))
  which(spans <= tol * scales)
}

#' Evaluate the lowest-category coefficient moment
#' @keywords internal
.semiparametric_beta_moment <- function(Y, X, b, Ftilde, free_index) {
  if (!length(free_index)) return(numeric(0))
  residual <- as.integer(Y == 0L) - pmin(pmax(Ftilde(as.numeric(X %*% b)), 0), 1)
  as.numeric(colMeans(X[, free_index, drop = FALSE] * residual))
}

#' Update identified first-stage coefficients for a fixed regularized CDF
#'
#' The normalized coefficient and all constant columns are held fixed. Holding
#' constant columns is an additional location normalization required when a
#' simulation design explicitly includes an intercept with an unrestricted F.
#'
#' @keywords internal
.update_beta_npmle_moment <- function(Y, X, b_start, Ftilde, ftilde = NULL,
                                      fixed_index = 1L, maxit = 500L,
                                      tol = 1e-5, verbose = FALSE) {
  Y <- as.integer(Y)
  X <- as.matrix(X)
  p <- ncol(X)
  if (anyNA(Y) || nrow(X) != length(Y) || any(!is.finite(X)) ||
      length(b_start) != p || any(!is.finite(b_start)) ||
      fixed_index < 1L || fixed_index > p) {
    stop("b_start and fixed_index are incompatible with X.")
  }

  constant_index <- .semiparametric_constant_columns(X)
  if (fixed_index %in% constant_index) {
    stop("The covariate selected by fixed_index must not be constant.")
  }
  free_index <- setdiff(seq_len(p), c(fixed_index, constant_index))
  fixed_value <- b_start[fixed_index]

  compose_b <- function(free) {
    b <- b_start
    b[fixed_index] <- fixed_value
    b[free_index] <- free
    b
  }
  moment <- function(free) {
    .semiparametric_beta_moment(
      Y, X, compose_b(free), Ftilde, free_index
    )
  }
  objective <- function(free) sum(moment(free)^2)
  gradient <- if (!is.null(ftilde)) {
    function(free) {
      b <- compose_b(free)
      density <- pmax(ftilde(as.numeric(X %*% b)), 0)
      X_free <- X[, free_index, drop = FALSE]
      jacobian <- -crossprod(X_free, X_free * density) / nrow(X)
      as.numeric(2 * crossprod(jacobian, moment(free)))
    }
  } else {
    NULL
  }

  start_free <- b_start[free_index]
  start_score <- moment(start_free)
  start_objective <- sum(start_score^2)
  if (!length(free_index)) {
    return(list(
      b = compose_b(numeric(0)), converged = TRUE, updated = FALSE,
      score = numeric(0), score_norm = 0, start_score_norm = 0,
      objective = 0, start_objective = 0, convergence = 0L,
      message = "no identified free coefficients",
      free_index = free_index, constant_index = constant_index
    ))
  }

  fit <- tryCatch(
    stats::optim(
      par = start_free,
      fn = objective,
      gr = gradient,
      method = "BFGS",
      control = list(maxit = maxit, reltol = 1e-12)
    ),
    error = function(e) NULL
  )
  candidates <- list(list(
    par = start_free, value = start_objective, convergence = NA_integer_,
    message = "starting value"
  ))
  if (!is.null(fit) && is.finite(fit$value)) candidates[[2L]] <- fit

  if (is.null(fit) || !is.finite(fit$value) ||
      max(abs(moment(fit$par))) > tol) {
    fallback <- tryCatch(
      if (length(start_free) == 1L) {
        radius <- 5 * (1 + abs(start_free))
        one_dimensional <- stats::optimize(
          objective,
          interval = c(start_free - radius, start_free + radius),
          tol = 1e-10
        )
        list(
          par = one_dimensional$minimum,
          value = one_dimensional$objective,
          convergence = 0L,
          message = "bounded one-dimensional fallback"
        )
      } else {
        stats::optim(
          par = start_free,
          fn = objective,
          method = "Nelder-Mead",
          control = list(maxit = maxit, reltol = 1e-12)
        )
      },
      error = function(e) NULL
    )
    if (!is.null(fallback) && is.finite(fallback$value)) {
      candidates[[length(candidates) + 1L]] <- fallback
    }
  }

  values <- vapply(candidates, function(x) x$value, numeric(1))
  selected <- candidates[[which.min(values)]]
  final_score <- moment(selected$par)
  score_norm <- max(abs(final_score))
  b <- compose_b(selected$par)
  b[fixed_index] <- fixed_value

  if (verbose) {
    message(sprintf(
      "Step 1 beta moment: %.3e -> %.3e",
      if (length(start_score)) max(abs(start_score)) else 0,
      score_norm
    ))
  }

  list(
    b = b,
    converged = is.finite(score_norm) && score_norm <= tol,
    updated = selected$value < start_objective - 1e-12,
    score = final_score,
    score_norm = score_norm,
    start_score_norm = max(abs(start_score)),
    objective = selected$value,
    start_objective = start_objective,
    convergence = selected$convergence,
    message = selected$message,
    free_index = free_index,
    constant_index = constant_index
  )
}

#' Update the positive second cutoff for a fixed regularized CDF
#' @keywords internal
.estimate_alpha_smooth <- function(Y00, X00, b_00, Ftilde,
                                   alpha_start = NULL, tol = 1e-8) {
  u <- as.numeric(as.matrix(X00) %*% b_00)
  target <- as.integer(Y00 <= 1L)
  score <- function(alpha) mean(target - Ftilde(u + alpha))

  lower <- 0
  lower_score <- score(lower)
  if (!is.finite(lower_score) || lower_score <= 0) return(NA_real_)

  scale <- stats::sd(u)
  if (!is.finite(scale) || scale <= 0) scale <- 1
  upper <- max(1e-3, scale, if (is.null(alpha_start)) 0 else 1.5 * alpha_start)
  upper_score <- score(upper)
  expansions <- 0L
  while ((!is.finite(upper_score) || upper_score > 0) && expansions < 40L) {
    upper <- 2 * upper
    upper_score <- score(upper)
    expansions <- expansions + 1L
  }
  if (!is.finite(upper_score) || upper_score > 0) return(NA_real_)

  tryCatch(
    stats::uniroot(score, lower = lower, upper = upper, tol = tol)$root,
    error = function(e) NA_real_
  )
}

#' Run the alternating NPMLE, coefficient, and cutoff first stage
#'
#' A direct I-spline is used only to make the NPMLE moment map numerically
#' smooth. The final reporting regularization is selected after this common
#' first stage.
#' @keywords internal
.run_semiparametric_first_stage <- function(
    Y00, X00, fixed_index = 1L, min_iter = 2L, max_iter = 20L,
    parameter_tol = 5e-3, moment_tol = 1e-3, stable_passes = 2L,
    verbose = FALSE) {
  init <- .initialize_semiparametric_j3(Y00, X00, fixed_index)
  b <- init$b
  alpha <- init$alpha
  n0 <- length(Y00)
  trace <- vector("list", max_iter)
  stable_count <- 0L
  converged <- FALSE
  cycle_detected <- FALSE
  stop_reason <- "max_iter"

  for (iteration in seq_len(max_iter)) {
    npmle <- .compute_joint_npmle(Y00, X00, b, alpha)
    if (is.null(npmle)) {
      stop("The interval-censoring NPMLE failed during the first-stage update.")
    }
    update_smoother <- .smooth_Fhat(npmle, n0 = n0)
    beta_update <- .update_beta_npmle_moment(
      Y = Y00,
      X = X00,
      b_start = b,
      Ftilde = update_smoother$Ftilde,
      ## The direct spline is a proper monotone approximation on its fitted
      ## support, but the imposed 0/1 tails create finite endpoint jumps. Its
      ## ordinary density therefore is not a valid global derivative of the
      ## moment objective; use BFGS finite differences for this update.
      ftilde = NULL,
      fixed_index = fixed_index,
      tol = moment_tol,
      verbose = verbose
    )
    b_new <- beta_update$b
    alpha_new <- .estimate_alpha_smooth(
      Y00, X00, b_new, update_smoother$Ftilde,
      alpha_start = alpha
    )
    if (!is.finite(alpha_new) || alpha_new <= 0) {
      stop("The cutoff moment did not have a positive finite root.")
    }

    identified <- unique(c(
      beta_update$free_index,
      fixed_index
    ))
    change <- max(c(
      abs(b_new[identified] - b[identified]) / (1 + abs(b[identified])),
      abs(alpha_new - alpha) / (1 + abs(alpha))
    ))
    alpha_score <- mean(
      as.integer(Y00 <= 1L) -
        update_smoother$Ftilde(as.numeric(X00 %*% b_new) + alpha_new)
    )
    trace[[iteration]] <- list(
      iteration = iteration,
      b = b_new,
      alpha = alpha_new,
      parameter_change = change,
      beta_score_norm = beta_update$score_norm,
      alpha_score = alpha_score,
      beta_update = beta_update
    )
    b <- b_new
    alpha <- alpha_new

    is_stable <- change <= parameter_tol &&
      beta_update$score_norm <= moment_tol && abs(alpha_score) <= 1e-6
    stable_count <- if (is_stable) stable_count + 1L else 0L
    if (iteration >= min_iter && stable_count >= stable_passes) {
      converged <- TRUE
      stop_reason <- "stable"
      break
    }

    ## A Turnbull update is a discrete map and can enter a short numerical
    ## cycle even when every moment is small. Treat recurrence within the same
    ## practical tolerance as a diagnosed stopping condition rather than
    ## silently exhausting max_iter.
    if (iteration >= min_iter + 2L &&
        beta_update$score_norm <= moment_tol && abs(alpha_score) <= 1e-6) {
      candidate_lags <- 2L:min(5L, iteration - 1L)
      cycle_change <- vapply(candidate_lags, function(lag) {
        prior <- trace[[iteration - lag]]
        max(c(
          abs(b[identified] - prior$b[identified]) /
            (1 + abs(prior$b[identified])),
          abs(alpha - prior$alpha) / (1 + abs(prior$alpha))
        ))
      }, numeric(1))
      if (any(cycle_change <= parameter_tol)) {
        cycle_lag <- candidate_lags[which.min(cycle_change)]
        cycle_detected <- TRUE
        stop_reason <- paste0("cycle_", cycle_lag)
        break
      }
    }
  }
  trace <- trace[!vapply(trace, is.null, logical(1))]

  final_npmle <- .compute_joint_npmle(Y00, X00, b, alpha)
  if (is.null(final_npmle)) {
    stop("The final NPMLE refit failed at the updated coefficients and cutoff.")
  }
  final_smoother <- .smooth_Fhat(final_npmle, n0 = n0)
  constant_index <- .semiparametric_constant_columns(X00)
  free_index <- setdiff(seq_len(ncol(X00)), c(fixed_index, constant_index))
  final_beta_score <- .semiparametric_beta_moment(
    Y00, X00, b, final_smoother$Ftilde, free_index
  )
  final_alpha_score <- mean(
    as.integer(Y00 <= 1L) -
      final_smoother$Ftilde(as.numeric(X00 %*% b) + alpha)
  )
  final_scores_practical <-
    (if (length(final_beta_score)) max(abs(final_beta_score)) else 0) <=
      parameter_tol && abs(final_alpha_score) <= parameter_tol
  converged <- converged && final_scores_practical
  if (stop_reason == "stable" && !final_scores_practical) {
    stop_reason <- "stable_update_but_final_score_large"
  }

  list(
    b = b,
    beta = -b,
    alpha = alpha,
    npmle = final_npmle,
    iterations = length(trace),
    completed = TRUE,
    converged = converged,
    fixed_point_stable = converged,
    cycle_detected = cycle_detected,
    stop_reason = stop_reason,
    trace = trace,
    initializer = init,
    fixed_index = fixed_index,
    constant_index = constant_index,
    final_beta_score = final_beta_score,
    final_alpha_score = final_alpha_score
  )
}
