## Adapt_rho - R-function implementing adaptive step-size for ADMM

## Input: rho    [numeric]: (Fixed) ADMM step size
##        r_norm [numeric]: L2-norm of primal residuals
##        s_norm [numeric]: L2-norm of dual residuals
##           tau [numeric]: Scalar
##           eta [numeric]: Scalar

## Output:   rho [numeric]: Adaptive ADMM step size

Adapt_rho <- function(rho, r_norm, s_norm, tau = 2, eta = 10){
  rho <- ifelse(r_norm > eta * s_norm, tau*rho,
                ifelse(r_norm < eta * s_norm, rho / tau, rho))
  return(rho)
}
