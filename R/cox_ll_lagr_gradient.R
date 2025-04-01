#' cox_ll_lagr_gradient
#'
#' R-function calculating the gradient of the augmented Lagrangian form of full/partial log-likelihood
#'
#' @param X [matrix]: Design matrix of dimension n_long x p_vars
#' @param d [data frame]: Data set with variables Tstart, Tstop, trans and status
#' @param beta [vector]: Regression parameter
#' @param Riskset [list]: Risk set list
#' @param rho [numeric]: Augmented Lagrangian parameter (ADMM step size; default: 1)
#' @param theta [vector]: ADMM parameter theta (dimension M x 1)
#' @param nu [vector]: ADMM parameter nu (dimension M x 1)
#'
#' @returns scorevector: Gradient at beta
#' @export
#'
#' @examples

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
