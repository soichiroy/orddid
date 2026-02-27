#' Visualizing the result of equivalence_test()
#'
#' \code{plot.orddid.test} generates a plot that visualize the result of \code{\link{equivalence_test}}.
#'
#' @param x An output from \code{equivalence_test()}.
#' @param ylim Range of y-axis of the plot.
#' @param fill A boolean argument; if \code{TRUE}, confidence band is filled by a polygon.
#' @param ... Other arguments passed to \code{plot()} function.
#' @importFrom graphics polygon lines abline text legend
#' @export
plot.orddid.test <- function(x, ylim, fill = TRUE, ...) {
  ## extract infor
  v_range <- c(0, x$v_range, 1)
  tv      <- c(0, x$tv, 0)
  Uv      <- c(0, x$Uv, 0)
  Lv      <- c(0, x$Lv, 0)
  ULmax   <- max(abs(x$Umax), abs(x$Lmin))
  ethres  <- attr(x, "threshold")
  tmax    <- round(x$tmax, 3)

  ## plot
  plot(v_range, tv, type = 'l',
        xlab = "Quantile (v)",
        ylab = expression(hat(t)~"(v) ="~tilde(q)[1]~"(v)"~-~tilde(q)[0]~"(v)"),
        main = "Test Statistic \n (Pre-Treatment Outcome)",
        ylim = ylim, lwd = 1.5, ...)
  if (isTRUE(fill)) {
    # fill the area covered by the CI
    polygon(c(v_range, rev(v_range)), c(Uv, rev(Lv)), col = '#d3d4d2',
            border = NA)
    lines(v_range, tv, col = 'black', lwd = 1.5)
  }
  lines(v_range, Uv, lty = 2, lwd = 1.4, col = 'gray30')
  lines(v_range, Lv, lty = 2, lwd = 1.4, col = 'gray30')

  abline(h = c(-ethres, ethres), col = 'red')
  abline(h = c(-ULmax, ULmax), col = '#1E88A8')
  text(0.4, ethres, paste("equivalence threshold =", round(ethres, 3)),
        col = 'red', pos = 3, cex = 0.8)
  text(0.4, -ULmax, paste("minimum threshold =", round(ULmax, 3)),
        col = '#1E88A8', pos = 1, cex = 0.8)
  legend("topleft", legend = bquote(widehat(t)["max"]~"="~.(tmax)),
         cex = 0.9, bty = 'n')
}
