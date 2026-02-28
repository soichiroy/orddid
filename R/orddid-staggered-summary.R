## ============================================================
## orddid-staggered-summary.R
## S3 methods for orddid_staggered objects
## ============================================================

#' @export
summary.orddid_staggered <- function(object, ...) {
  cat("Ordinal DID: Staggered Adoption Design\n")
  cat("=======================================\n\n")

  ## Group-time effects
  gt <- object$group_time_effects
  groups <- unique(gt$group)
  times  <- unique(gt$time)

  cat("Group-Time Category-Specific Effects (zeta):\n\n")
  for (gg in groups) {
    for (tt in unique(gt$time[gt$group == gg])) {
      sub <- gt[gt$group == gg & gt$time == tt, ]
      cat(sprintf("  (G=%s, T=%s):\n", gg, tt))
      for (r in seq_len(nrow(sub))) {
        cat(sprintf("    Category %d: %7.4f  [%7.4f, %7.4f]\n",
                    sub$category[r], sub$effect[r],
                    sub$lower_ci[r], sub$upper_ci[r]))
      }
    }
  }

  ## Group-time tau
  tau <- object$group_time_tau
  cat("\nGroup-Time Relative Effects (tau bounds):\n\n")
  for (r in seq_len(nrow(tau))) {
    cat(sprintf("  (G=%s, T=%s): [%7.4f, %7.4f]  CI: [%7.4f, %7.4f]\n",
                tau$group[r], tau$time[r],
                tau$tau_lb[r], tau$tau_ub[r],
                tau$ci_lower[r], tau$ci_upper[r]))
  }

  ## Aggregate
  agg <- object$aggregate_effects
  if (!is.null(agg$agg_zeta)) {
    cat("\nAggregate Effects (uniform weights):\n\n")
    J <- length(agg$agg_zeta)
    for (j in seq_len(J)) {
      cat(sprintf("  zeta_%d: %7.4f  [%7.4f, %7.4f]\n",
                  j, agg$agg_zeta[j],
                  agg$agg_zeta_ci[j, 1], agg$agg_zeta_ci[j, 2]))
    }
    cat(sprintf("  tau:    [%7.4f, %7.4f]  CI: [%7.4f, %7.4f]\n",
                agg$agg_tau[1], agg$agg_tau[2],
                agg$agg_tau_ci[1], agg$agg_tau_ci[2]))
  }

  invisible(object)
}

#' @export
print.orddid_staggered <- function(x, ...) {
  n_gt <- nrow(x$group_time_tau)
  cat(sprintf("Ordinal DID (Staggered): %d group-time pairs, %d bootstrap iterations\n",
              n_gt, x$inputs$n_boot))
  invisible(x)
}
