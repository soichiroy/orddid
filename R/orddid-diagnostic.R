#
# title: functions for diagnostic test 
#
# author: Soichiro Yamauchi 
# 
# note: 
#   steps for the test 
#    1. Estimate all parameters using ordered probit 
#    2. Obtain variance of parameters using bootstrap
#    3. Compute t(v) using qd() function 
#    4. Compute Var(t(v))
#




## --------------------------------------------------------------- ##
##                     Auxiliary Functions                         ##
## --------------------------------------------------------------- ##

#' Compute qd(v)
#'
#' @importFrom VGAM erf
qd <- function(u, mu1, mu0, s1, s0) {
  B <- erf(2 * u - 1, inverse = TRUE) / (s0 / s1)
  A <- erf((mu1 - mu0) / (s0 * sqrt(2)) + B)
  return((1 + A) / 2)
}


#' Compute tvmax 
#'
#' Compute max_v | q1(v) - q0(v) | 
#' @param theta a list of parameters 
tv_max <- function(theta, grid = 0.01) {
  v_range <- seq(0.001, 0.999, by = grid)
  tv <- rep(NA, length(v_range))
  for (v in 1:length(v_range)) {
    tv[v] <- qd(v_range[v], mu1 = theta$mu11, mu0 = theta$mu10,
               s1 = theta$sd11, s0 = theta$sd10) - 
            qd(v_range[v], mu1 = theta$mu01, mu0 = theta$mu00,
                       s1 = theta$sd01, s0 = theta$sd00)       
  }
  
  return(max(tv))
}


#' Compute zd(v)
#'
#' @importFrom VGAM erf
zd <- function(v, mu1, mu0, s1, s0) {
  A <- (mu1 - mu0) / (s0 * sqrt(2)) 
  B <- erf(2 * v - 1, inverse = TRUE) / (s0 / s1)
  return(A + B)
}



#' Compute gradient of t(v; theta)
#' @param theta a list of parameters
#' @importFrom VGAM erf
tv_gradient <- function(v, theta) {
  ## compute intermediate quantities 
  z0 <- zd(v, mu1 = theta$mu01, mu0 = theta$mu00, s1 = theta$sd01, s0 = theta$sd00)
  z1 <- zd(v, mu1 = theta$mu11, mu0 = theta$mu10, s1 = theta$sd11, s0 = theta$sd10)
  exp_z0 <- exp(-z0^2)
  exp_z1 <- exp(-z1^2)
  erf_1  <- erf(2 * v - 1, inverse = TRUE)
  sqrt_pi <- sqrt(pi)
  sqrt_2pi <- sqrt(2 * pi)
  
  ## compute the gradient 
  grad <- c(
    ## derivative wrt theta 0
    exp_z0 / (sqrt_2pi * theta$sd00),
    exp_z0 * z0  / (sqrt_pi * theta$sd00),
    -exp_z0 / (sqrt_2pi * theta$sd00),
    -exp_z0 * erf_1 / (sqrt_pi * theta$sd00),
    ## wrt theta 1
    -exp_z1 / (sqrt_2pi * theta$sd10),
    -exp_z1 * z1  / (sqrt_pi * theta$sd10),
    exp_z1 / (sqrt_2pi * theta$sd10),
    exp_z1 * erf_1 / (sqrt_pi * theta$sd10)
  )
  
  return(grad)
}
