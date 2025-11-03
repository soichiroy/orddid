

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
#' @noRd
#' @keywords internal
.ComputeRelativeEffect <- function(p1, p0) {
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
