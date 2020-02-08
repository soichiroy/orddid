#' input check for \link{\code{fit_ord_probit}}.
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
