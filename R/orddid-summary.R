
#' Summary function 
#' @export
summary.orddid <- function(obj) {

  if ("orddid.fit" %in% class(obj)) {
    res <- summarise_ord_did(obj)
    ses <- sqrt(do.call("rbind", res$Var))
    cis <- do.call("rbind", res$CI)
    del <- unlist(res$Delta)
    J   <- attr(obj, "n_choice")
    
    # create a table for output 
    tab <- cbind(del, ses, cis)
    colnames(tab) <- c("Effect", "SE", 
      "90% Lower", "90% Upper", "95% Lower", "95% Upper")
    rownames(tab) <- c(paste("Delta[", 2:(J-1), "-", J,"]", sep = ""),
                      paste("Delta[", J, "]", sep = ""))
  }
  
  class(tab) <- 'summary.orddid'
  return(tab)
}

#' Print function
#' @importFrom cli cat_rule
#' @export
print.summary.orddid <- function(obj) {
  class(obj) <- NULL
  cat_rule(left = crayon::bold("Effect Estimates"))
  print.default(obj, quote = FALSE, right = TRUE, digits = 3)
  
  invisible(obj)
}

#' Summarise function for ord_did()
#' @keywords internal 
summarise_ord_did <- function(obj) {
  if ("orddid.fit" %in% class(obj)) {
    # extract the number of chices 
    J      <- attr(obj, "n_choice")
    n_boot <- attr(obj, "input")$n_boot
    
    # loop over all possible effects 
    CIs <- Delta <- Var <- list()
    for (j in 1:(J-1)) {
      # compute CI and variance from the bootstrap resutls 
      boot_dist   <- sapply(1:n_boot, function(i) {
        sum(obj$boot[[1]][i,(j+1):J]) - sum(obj$boot[[2]][i,(j+1):J])
      })    
      CIs[[j]] <- quantile(boot_dist, prob = c(0.05, 0.95, 0.025, 0.975))
      Var[[j]] <- var(boot_dist)
      
      # compute the effect 
      Delta[[j]] <- sum(obj$fit$Y1[(j+1):J]) - sum(obj$fit$Y0[(j+1):J])
    }
  } else {
    cat("Not a supported input!")
  }
  
  return(list("CI" = CIs, "Var" = Var, "Delta" = Delta))
}
