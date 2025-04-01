#' Adapt_rho
#'
#' R-function implementing adaptive step-size for ADMM
#'
#' @param rho [numeric]: (Fixed) ADMM step size
#' @param r_norm [numeric]: L2-norm of primal residuals
#' @param s_norm [numeric]: L2-norm of dual residuals
#' @param tau [numeric]: Scalar
#' @param eta [numeric]: Scalar
#'
#' @returns rho [numeric]: Adaptive ADMM step size
#' @export
#'
#' @examples

Adapt_rho <- function(rho, r_norm, s_norm, tau = 2, eta = 10){
  rho <- ifelse(r_norm > eta * s_norm, tau*rho,
                ifelse(r_norm < eta * s_norm, rho / tau, rho))
  return(rho)
}
