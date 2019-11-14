# Plot Equivalence Test
#'
#' @export
plot.orddid.test <- function(obj, ylim, ...) {
  ## extract infor
  v_range <- c(0, obj$v_range, 1)
  tv      <- c(0, obj$tv, 0)
  Uv      <- c(0, obj$Uv, 0)
  Lv      <- c(0, obj$Lv, 0)
  ULmax   <- max(abs(obj$Umax), abs(obj$Lmin))
  ethres  <- attr(obj, "threshold")
  tmax    <- round(obj$tmax, 3)

  ## plot
  plot(v_range, tv, type = 'l',
        xlab = "v",
        ylab = expression(hat(t)~"(v) ="~tilde(q)[1]~"(v)"~-~tilde(q)[0]~"(v)"),
        main = "Test Statistic \n (Pre-Treatment Outcome)",
        ylim = ylim, lwd = 1.5)
  lines(v_range, Uv, lty = 2, lwd = 1.3, col = 'gray40')
  lines(v_range, Lv, lty = 2, lwd = 1.3, col = 'gray40')
  abline(h = c(-ethres, ethres), col = 'red')
  abline(h = c(-ULmax, ULmax), col = '#33A6B8')
  text(0.4, ethres, paste("equivalence threshold =", round(ethres, 3)),
        col = 'red', pos = 3, cex = 0.8)
  text(0.4, -ULmax, paste("minimum threshold =", round(ULmax, 3)),
        col = '#33A6B8', pos = 1, cex = 0.8)
  legend("topleft", legend = bquote(widehat(t)["max"]~"="~.(tmax)),
         cex = 0.9, bty = 'n')
}
