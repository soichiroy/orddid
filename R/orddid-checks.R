#' Input check for \link{\code{fit_ord_probit}}.
#'
#' This function check if there is missing categories. For example, if Ymin = 0 and Ymax=J-1.
#'  the function checks there are J unique values in Y.
#' @param Y a vector of ategorical outcome.
#' @keywords internal
ord_probit_check_cat <- function(Y) {
  j_min <- min(Y)
  j_max <- max(Y)
  j_uq  <- length(unique(Y))
  if (j_uq != ((j_max - j_min) + 1)) {
    stop("Outcome has a missing categories.")
  }
}


#' Input check for \link{\code{ord_did_run}}.
#'
#' @param Ynew Ynew from \code{ord_did_run()}.
#' @param Yold Yold from \code{ord_did_run()}.
#' @param treat treat from \code{ord_did_run()}.
#' @keywords internal
ord_did_check_input <- function(Ynew, Yold, treat) {

  ##
  ## check length
  ##
  n1 <- length(Ynew)
  n0 <- length(Yold)
  nT <- length(treat)

  if (!identical(n1, n0)) {
    stop("Length of Ynew and Yold do not match.")
  }

  if (!identical(n1, nT)) {
    stop("Length of Ynew and treat do not match.")
  }

  if (!identical(n0, nT)) {
    stop("Length of Yold and treat do not match.")
  }


  ##
  ## check na
  ##
  na1 <- sum(is.na(Ynew))
  na0 <- sum(is.na(Yold))
  naT <- sum(is.na(treat))

  if (na1 > 0) {
    stop("NA is not allowed in Ynew.")
  }
  if (na0 > 0) {
    stop("NA is not allowed in Yold")
  }
  if (naT > 0) {
    stop("NA is not allowed in treat")
  }

  ##
  ## input has to be numeric
  ##
  nm1 <- is.numeric(Ynew)
  nm0 <- is.numeric(Yold)
  nmT <- is.numeric(treat)
  if (isFALSE(nm1)) {
    stop('Ynew has to be numeric.')
  }
  if (isFALSE(nm0)) {
    stop("Yold has to be numeric.")
  }
  if (isFALSE(nmT)) {
    stop("treat has to be numeric.")
  }
}
