# Gradient of augmented Lagrangian of partial log-likelihood:

## cox_ll_lagr_gradient - R-function calculating the gradient of the augmented Lagrangian
##                        form of partial log-likelihood
##
## Input: X             [matrix]: Regression matrix of dimension n_obs x p_vars
##        d         [data frame]: Data set with variables Tstart, Tstop, trans and status
##        beta          [vector]: Regression parameter
##        Riskset         [list]: Risk set list
##        rho          [numeric]: Augmented Lagrangian parameter (step size; default: 1)
##        theta        [numeric]: ADMM parameter theta (dimension M x 1)
##        nu           [numeric]: ADMM parameter nu (dimension M x 1)
##
## Output: scorevector [numeric]: Gradient at beta

cox_ll_lagr_gradient <- function(X, d, beta, Riskset, rho = 1, theta, nu){

  X <- as.matrix(X)
  event <- d$status

  # Linear predictor
  lp <- X %*% beta
  ef <- exp(lp)

  # Initializing risk matrix
  n <- length(event)
  p <- length(beta)
  riskmatrix <- matrix(0, nrow = n, ncol = p)

  # Calculating risk matrix
  for (i in 1:n) {
    riskset <- Riskset[[i]]
    ef.riskset <- ef[riskset]
    currentrisk <- sum(ef.riskset)
    X.i <- X[riskset, ] / currentrisk
    riskmatrix[i, ] <- t(ef.riskset) %*% X.i
  }

  # Score vector calculation
  scorevector <- as.numeric(event %*% (X - riskmatrix)) + rho * (-theta - nu)

  return(scorevector)
}
