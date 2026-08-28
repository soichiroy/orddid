

.ComputeTreatmentEffects <- function(y1_prop, y0_prop) {
  return(y1_prop - y0_prop)
}

#' Sharp bounds on the relative treatment effect for ordinal outcomes
#'
#' @description
#' Computes the sharp lower and upper bounds \eqn{(\gamma_L,\gamma_U)} on the
#' relative treatment effect \eqn{\gamma=\Pr\{Y(1)>Y(0)\}-\Pr\{Y(1)<Y(0)\}}
#' for ordinal outcomes, given the marginal probability vectors
#' \eqn{p_1=(p_{0+},\ldots,p_{J-1,+})} and \eqn{p_0=(p_{+0},\ldots,p_{+,J-1})}.
#' The function implements the closed-form solutions in Lu, Zhang and Ding (2020),
#' using the quantities \eqn{\delta_{jm}} and \eqn{\xi_{jm}} defined for
#' \eqn{j=1,\ldots,J-1} and \eqn{m=1,\ldots,J-j}.
#'
#' @param p1 Numeric vector of length \eqn{J} with \eqn{p_{k+}=\Pr\{Y(1)=k\}}, ordered from the worst category (0) to the best (J-1).
#' @param p0 Numeric vector of length \eqn{J} with \eqn{p_{+l}=\Pr\{Y(0)=l\}}, ordered from the worst category (0) to the best (J-1).
#'
#' @details
#' Let \eqn{J=\mathrm{length}(p1)=\mathrm{length}(p0)\ge 3}. Define tail sums
#' \eqn{T_1(a)=\sum_{k=a}^{J-1}p_{k+}} and \eqn{T_0(b)=\sum_{l=b}^{J-1}p_{+l}}
#' (empty sums are zero), and head sums \eqn{H_1(c)=\sum_{k=0}^{c}p_{k+}},
#' \eqn{H_0(c)=\sum_{l=0}^{c}p_{+l}}. Then, for \eqn{j=1,\ldots,J-1} and \eqn{m=1,\ldots,J-j},
#' \deqn{\delta_{jm}=T_1(j)+T_1(j+m)+H_0(j-2)-T_0(j+m-1),}
#' \deqn{\xi_{jm}=T_1(j+m-1)-T_0(j)-T_0(j+m)-H_1(j-2).}
#' The sharp bounds are \eqn{\gamma_U=\min_{j,m}\delta_{jm}} and
#' \eqn{\gamma_L=\max_{j,m}\xi_{jm}}.
#'
#' @return A named numeric vector of length 2 with elements:
#' \describe{
#'   \item{lower}{\eqn{\gamma_L} — sharp lower bound}
#'   \item{upper}{\eqn{\gamma_U} — sharp upper bound}
#' }
#'
#' @references
#' Lu, J., Zhang, Y., & Ding, P. (2020). Sharp bounds on the relative treatment effect for ordinal outcomes.
#' \emph{Biometrics}, 76(2), 664–669. \doi{10.1111/biom.13148}
#'
#' @examples
#' ## J = 3 example:
#' p1 <- c(0.2, 0.3, 0.5)  # Pr{Y(1)=0,1,2}
#' p0 <- c(0.4, 0.4, 0.2)  # Pr{Y(0)=0,1,2}
#' bound_gamma_ordinal(p1, p0)
#'
#' @param tau_assumption Character vector. \code{"none"} uses the closed-form
#'   sharp bounds from Lu, Zhang, and Ding (2020). \code{"stochastic_dominance"}
#'   imposes positive regression dependence on the joint distribution of
#'   \eqn{(Y(1), Y(0))} and computes the resulting bounds by linear programming.
#'
#' @noRd
#' @keywords internal
.ComputeRelativeEffect <- function(p1, p0,
                                   tau_assumption = "none") {
  tau_assumption <- .orddid_match_tau_assumptions(tau_assumption)
  probs <- .orddid_validate_relative_effect_probs(p1 = p1, p0 = p0)

  out <- lapply(tau_assumption, function(assump) {
    bounds <- tryCatch(
      {
        if (assump == "none") {
          .ComputeRelativeEffectClosedForm(p1 = probs$p1, p0 = probs$p0)
        } else {
          .ComputeRelativeEffectStochasticDominance(p1 = probs$p1, p0 = probs$p0)
        }
      },
      error = function(e) c(lower = NA_real_, upper = NA_real_)
    )

    tibble::tibble(
      tau_assumption = assump,
      effect_lb = bounds[["lower"]],
      effect_ub = bounds[["upper"]],
      status = if (all(is.finite(bounds))) "ok" else "infeasible"
    )
  })

  dplyr::bind_rows(out)
}


.orddid_match_tau_assumptions <- function(
  tau_assumption = "none"
) {
  unique(match.arg(
    tau_assumption,
    choices = c("none", "stochastic_dominance"),
    several.ok = TRUE
  ))
}


.orddid_validate_relative_effect_probs <- function(p1, p0, tol = 1e-10) {
  if (!is.numeric(p1) || !is.numeric(p0)) {
    stop("p1 and p0 must be numeric probability vectors.", call. = FALSE)
  }
  if (length(p1) != length(p0)) {
    stop("p1 and p0 must have the same length.", call. = FALSE)
  }
  if (length(p1) < 2) {
    stop("Relative-effect bounds require at least two outcome categories.",
         call. = FALSE)
  }
  if (any(!is.finite(p1)) || any(!is.finite(p0))) {
    stop("p1 and p0 must contain only finite values.", call. = FALSE)
  }
  if (any(p1 < -tol) || any(p0 < -tol)) {
    stop("p1 and p0 must be nonnegative probability vectors.", call. = FALSE)
  }

  p1[p1 < 0] <- 0
  p0[p0 < 0] <- 0

  s1 <- sum(p1)
  s0 <- sum(p0)
  if (s1 <= tol || s0 <= tol) {
    stop("p1 and p0 must each sum to a positive value.", call. = FALSE)
  }

  list(p1 = p1 / s1, p0 = p0 / s0)
}


.ComputeRelativeEffectClosedForm <- function(p1, p0) {
  J <- length(p1)

  # Tail sums with zero padding at the end: tail[i] = sum_{i..J}
  tail1 <- c(rev(cumsum(rev(p1))), 0)
  tail0 <- c(rev(cumsum(rev(p0))), 0)

  # Head sums with zero padding at the beginning: head_pre[i] = sum_{1..(i-1)}
  head1_pre <- c(0, cumsum(p1))  # index i gives sum_{1..(i-1)}
  head0_pre <- c(0, cumsum(p0))

  deltas <- c()
  xis    <- c()

  for (j in 1:(J - 1)) {
    for (m in 1:(J - j)) {
      # Map formula indices (0-based categories) to 1-based vector indices via the helpers above
      # δ_{jm} = T1(j) + T1(j+m) + H0(j-2) - T0(j+m-1)
      delta_jm <- tail1[j + 1] +                 # T1(j)
                  tail1[j + m + 1] +             # T1(j+m)
                  head0_pre[j] -                 # H0(j-2) because head0_pre[j] = sum_{1..(j-1)}
                  tail0[j + m]                   # T0(j+m-1) since tail0[i] = sum_{i..J}, with i=j+m

      # ξ_{jm} = T1(j+m-1) - T0(j) - T0(j+m) - H1(j-2)
      xi_jm <- tail1[j + m] -                    # T1(j+m-1)
               tail0[j + 1] -                    # T0(j)
               tail0[j + m + 1] -                # T0(j+m)
               head1_pre[j]                      # H1(j-2) because head1_pre[j] = sum_{1..(j-1)}

      deltas <- c(deltas, delta_jm)
      xis    <- c(xis, xi_jm)
    }
  }

  c(lower = max(xis), upper = min(deltas))
}


.ComputeRelativeEffectStochasticDominance <- function(p1, p0) {
  J <- length(p1)
  idx_mat <- matrix(seq_len(J * J), nrow = J, ncol = J)
  obj <- as.vector(sign(row(idx_mat) - col(idx_mat)))

  const_rows <- list()
  const_dir <- character()
  const_rhs <- numeric()

  add_constraint <- function(coeffs, dir, rhs) {
    const_rows[[length(const_rows) + 1L]] <<- coeffs
    const_dir <<- c(const_dir, dir)
    const_rhs <<- c(const_rhs, rhs)
  }

  for (j in seq_len(J)) {
    coeffs <- numeric(J * J)
    coeffs[idx_mat[j, ]] <- 1
    add_constraint(coeffs, "==", p1[j])
  }

  for (k in seq_len(J)) {
    coeffs <- numeric(J * J)
    coeffs[idx_mat[, k]] <- 1
    add_constraint(coeffs, "==", p0[k])
  }

  for (k_low in seq_len(J - 1L)) {
    k_high <- k_low + 1L
    if (p0[k_low] <= 0 || p0[k_high] <= 0) next

    for (threshold in seq_len(J - 1L)) {
      coeffs <- numeric(J * J)
      coeffs[idx_mat[seq_len(threshold), k_low]] <- p0[k_high]
      coeffs[idx_mat[seq_len(threshold), k_high]] <- -p0[k_low]
      add_constraint(coeffs, ">=", 0)
    }
  }

  for (j_low in seq_len(J - 1L)) {
    j_high <- j_low + 1L
    if (p1[j_low] <= 0 || p1[j_high] <= 0) next

    for (threshold in seq_len(J - 1L)) {
      coeffs <- numeric(J * J)
      coeffs[idx_mat[j_low, seq_len(threshold)]] <- p1[j_high]
      coeffs[idx_mat[j_high, seq_len(threshold)]] <- -p1[j_low]
      add_constraint(coeffs, ">=", 0)
    }
  }

  const_mat <- do.call(rbind, c(const_rows, list(diag(J * J))))
  const_dir <- c(const_dir, rep(">=", J * J))
  const_rhs <- c(const_rhs, rep(0, J * J))

  lp_upper <- lpSolve::lp(
    direction = "max",
    objective.in = obj,
    const.mat = const_mat,
    const.dir = const_dir,
    const.rhs = const_rhs
  )
  lp_lower <- lpSolve::lp(
    direction = "min",
    objective.in = obj,
    const.mat = const_mat,
    const.dir = const_dir,
    const.rhs = const_rhs
  )

  if (lp_upper$status != 0 || lp_lower$status != 0) {
    stop(
      "The stochastic_dominance tau assumption is infeasible for the supplied marginals.",
      call. = FALSE
    )
  }

  c(lower = lp_lower$objval, upper = lp_upper$objval)
}

#' Count pooled outcome categories across DID cells
#'
#' @param cells Named list of outcome vectors.
#' @return Integer number of observed pooled outcome categories.
#' @keywords internal
#' @noRd
.orddid_pooled_category_count <- function(cells) {
  if (!is.list(cells) || length(cells) == 0L) {
    stop("cells must be a non-empty list of outcome vectors.", call. = FALSE)
  }

  pooled <- unlist(cells, use.names = FALSE)
  pooled <- pooled[!is.na(pooled)]
  n_categories <- length(unique(pooled))
  if (n_categories < 2L) {
    stop("Ordinal DID requires at least two pooled outcome categories.",
         call. = FALSE)
  }
  n_categories
}


#' Count observed outcome categories in required DID cells
#'
#' @param cells Named list of outcome vectors.
#' @param required_cells Character vector naming cells to inspect.
#' @return Named integer vector of observed category counts.
#' @keywords internal
#' @noRd
.orddid_category_support_counts <- function(
    cells, required_cells = names(cells)) {
  if (!is.list(cells) || is.null(names(cells))) {
    stop("cells must be a named list of outcome vectors.", call. = FALSE)
  }
  if (!is.character(required_cells) || length(required_cells) == 0L ||
      anyNA(required_cells) || any(!required_cells %in% names(cells))) {
    stop("required_cells must identify one or more supplied cells.",
         call. = FALSE)
  }

  vapply(cells[required_cells], function(y) {
    length(unique(y[!is.na(y)]))
  }, integer(1))
}


#' Check fixed ordinal support in required DID cells
#'
#' Bootstrap resampling can remove a rare outcome category. Estimators that
#' infer cutoffs from the observed categories must reject such draws rather
#' than silently relabel the remaining categories.
#'
#' @param cells Named list of outcome vectors.
#' @param n_categories Number of categories in the original estimation sample.
#' @param required_cells Character vector naming cells that must retain support.
#' @return Logical scalar.
#' @keywords internal
#' @noRd
.orddid_has_full_category_support <- function(
    cells, n_categories, required_cells = names(cells)) {
  counts <- .orddid_category_support_counts(cells, required_cells)
  all(counts == as.integer(n_categories))
}


#' Require fixed ordinal support in specified DID cells
#'
#' @inheritParams .orddid_has_full_category_support
#' @param context Short description used in the error message.
#' @return The category counts, invisibly.
#' @keywords internal
#' @noRd
.orddid_require_full_category_support <- function(
    cells, n_categories, required_cells = names(cells),
    context = "The sample") {
  if (length(n_categories) != 1L || !is.numeric(n_categories) ||
      is.na(n_categories) || !is.finite(n_categories) ||
      n_categories < 2L || n_categories != floor(n_categories)) {
    stop("n_categories must be one integer greater than or equal to two.",
         call. = FALSE)
  }

  counts <- .orddid_category_support_counts(cells, required_cells)
  deficient <- counts != as.integer(n_categories)
  if (any(deficient)) {
    detail <- paste0(
      names(counts)[deficient], " (", counts[deficient], "/",
      as.integer(n_categories), ")",
      collapse = ", "
    )
    stop(
      context, " must retain all ", as.integer(n_categories),
      " pooled outcome categories in the required DID cells; deficient ",
      "cells: ", detail, ".",
      call. = FALSE
    )
  }

  invisible(counts)
}


#' Extract outcome vectors by post/treat cells
#'
#' Given a data frame and the column names for the outcome, post indicator, and
#' treatment indicator (all logical/boolean), this helper returns the four
#' vectors commonly used in DID with two periods and two groups:
#' - Y11: post & treated
#' - Y10: pre  & treated
#' - Y01: post & control
#' - Y00: pre  & control
#'
#' This function performs only minimal handling of missing values (it removes
#' rows with missing outcome or missing indicators when `remove_na = TRUE`).
#' It assumes the inputs are already validated elsewhere (same convention as
#' other internal helpers in this package).
#'
#' @param data A data.frame-like object containing the variables
#' @param outcome Character scalar; column name for the outcome variable
#' @param post Character scalar; column name for the post-period indicator (logical)
#' @param treat Character scalar; column name for the treatment indicator (logical)
#' @param remove_na Logical; if TRUE (default), drop rows with NA in any of
#'   outcome, post, or treat before subsetting
#'
#' @return A named list with elements `Y11`, `Y10`, `Y01`, and `Y00`, each a
#'   vector of the outcome in the corresponding cell.
#'
#' @keywords internal
#' @export 
.orddid_extract_Y_cells <- function(data, outcome, post, treat, remove_na = TRUE) {
  y <- data[[outcome]]
  p <- data[[post]]
  t <- data[[treat]]

  if (remove_na) {
    keep <- !(is.na(y) | is.na(p) | is.na(t))
    y <- y[keep]
    p <- p[keep]
    t <- t[keep]
  }

  Y11 <- y[p & t]
  Y10 <- y[(!p) & t]
  Y01 <- y[p & (!t)]
  Y00 <- y[(!p) & (!t)]

  return(list(Y11 = Y11, Y10 = Y10, Y01 = Y01, Y00 = Y00))
}
