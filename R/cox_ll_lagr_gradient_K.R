# Gradient of augmented Lagrangian of full/partial log-likelihood with penalty matrix K:

## cox_ll_lagr_gradient_K - R-function calculating the gradient of the augmented Lagrangian
##                          form of full/partial log-likelihood with penalty matrix K
##
## Input: X             [matrix]: Design matrix of dimension n_long x p_vars
##        d         [data frame]: Data set with variables Tstart, Tstop, trans and status
##        K             [matrix]: Penalty matrix of dimension M x p (M=p+m+g)
##        beta          [vector]: Regression parameter
##        Riskset         [list]: Risk set list
##        rho          [numeric]: Augmented Lagrangian parameter (ADMM step size; default: 1)
##        theta        [numeric]: ADMM parameter theta (dimension M x 1)
##        nu           [numeric]: ADMM parameter nu (dimension M x 1)
##
## Output: scorevector [numeric]: Score function at beta

cox_ll_lagr_gradient_K <- function(X, d, K, beta, Riskset, rho = 1, theta, nu){

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
  scorevector <- t(as.numeric(event %*% (X - riskmatrix)) + (rho * (t(beta) %*% t(K) - t(theta)) - t(nu)) %*% K)

  return(scorevector)
}
