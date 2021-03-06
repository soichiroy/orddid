#' Input check for \code{fit_ord_probit}.
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


#' Input check for \code{ord_did}.
#'
#' @param Ynew Ynew from \code{ord_did()}.
#' @param Yold Yold from \code{ord_did()}.
#' @param treat treat from \code{ord_did()}.
#' @param id_cluster id_cluster from \code{ord_did()}.
#' @param n_boot n_boot from \code{ord_did()}.
#' @keywords internal
ord_did_check_input <- function(Ynew, Yold, treat, id_cluster, n_boot)  {
  if (!(n_boot >= 2)) {
    stop("n_boot has to be greater than 2.")
  }

  if (is.null(id_cluster)) {
    n1 <- length(na.omit(Ynew))
    n0 <- length(na.omit(Yold))
    nT <- length(na.omit(treat))
    if (!(length(unique(c(n1, n0, nT))) == 1)) {
      stop("Length of inputs does not match.")
    }
  } else {
    n1 <- length(na.omit(Ynew))
    n0 <- length(na.omit(Yold))
    nT <- length(na.omit(treat))
    nI <- length(na.omit(id_cluster))
    if (!(length(unique(c(n1, n0, nT, nI))) == 1)) {
      stop("Length of inputs does not match.")
    }
  }
}

#' Input check for \code{ord_did_run}.
#'
#' @param Ynew Ynew from \code{ord_did_run()}.
#' @param Yold Yold from \code{ord_did_run()}.
#' @param treat treat from \code{ord_did_run()}.
#' @keywords internal
ord_did_run_check_input <- function(Ynew, Yold, treat) {

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
