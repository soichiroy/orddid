#' Create a group indicator
#' The same-length condition is already checked by \code{ord_did_check_input}.
#' @keywords internal
#' @noRd
create_group_dummy <- function(Ynew, Yold, treat, pre = FALSE) {
  ## control
  Y00 <- Yold[treat == 0]
  Y01 <- Ynew[treat == 0]
  ## treat
  Y10 <- Yold[treat == 1]
  Y11 <- Ynew[treat == 1]

  ## concatenate outcome
  if (isTRUE(pre)) {
    Yvec <- c(Y00, Y01, Y10, Y11)
    id_group <- c(
      rep(1, length(Y00)),
      rep(2, length(Y01)),
      rep(3, length(Y10)),
      rep(4, length(Y11))
    )
  } else {
    Yvec <- c(Y00, Y01, Y10)
    id_group <- c(rep(1, length(Y00)), rep(2, length(Y01)), rep(3, length(Y10)))
  }
  return(list(Y = Yvec, id_group = id_group))
}
