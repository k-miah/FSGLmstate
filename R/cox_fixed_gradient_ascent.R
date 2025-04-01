## cox_fixed_gradient_ascent - R-function implementing fixed gradient ascent
##                             for estimation in the Cox PH model
##
## Input: X          [matrix]: Regression matrix of dimension n_obs x p_vars
##        d      [data frame]: Data set with variables Tstart, Tstop, trans and status
##        K          [matrix]: Penalty matrix of dimension M x p (M=p+m+g)

##        eps       [numeric]: Step size in [0,1] (default: .01)
##        beta.init  [vector]: Initial value of beta (default: 0)
##        Riskset      [list]: Risk set list
##        tolerance [numeric]: Tolerance for stopping criterion (default: 1e-6)
##        max_iter  [numeric]: Maximum number of iterations (default: 1000)
##        rho       [numeric]: Augmented Lagrangian parameter (step size; default: 1)
##        theta     [numeric]: ADMM parameter theta = beta
##        nu        [numeric]: ADMM parameter nu = theta - beta
##
## Output: res         [list]: Beta estimation for Cox model at stopping iteration 'iter'

cox_fixed_gradient_ascent <- function(X, d, K, eps = 0.01, beta.init = NULL, Riskset,
                                      tolerance = 1e-6, max_iter = 1000,
                                      rho = 1, theta, nu){

  # Initialize coefficient beta
  if(is.null(beta.init)){
    beta <- rep(0, ncol(X))
  } else{
    beta <- beta.init
  }

  # Iterate until convergence or maximum iterations reached
  it <- 1
  tol <- 1
  ll <- 0

  while(tol > tolerance && it < max_iter){

    # Update coefficients: Fixed gradient ascent step
    if(is.null(K)){
      gradient <- cox_ll_lagr_gradient(X, d = d, beta = beta, Riskset = Riskset,
                                       rho = rho, theta = theta, nu = nu)
    } else{
      gradient <- cox_ll_lagr_gradient_K(X, d = d, K = K, beta = beta, Riskset = Riskset,
                                         rho = rho, theta = theta, nu = nu)
    }
    beta_new <- beta + eps * gradient

    # Stopping criterion: partial log-likelihood
    ll_old <- ll
    ll <- partial_ll(beta = beta_new, X, d, risksetlist = Riskset)
    tol <- abs(ll - ll_old)

    # Check step size for overshooting
    if(it > 1 & ll_old > ll){
      # reduce step size (step-halving)
      eps <- eps/2
      beta_new <- beta
      ll <- ll_old
    }

    it <- it + 1
    beta <- beta_new
  }
  # If max_iter is reached, print warning and exit loop
  if(it == max_iter){
    warning(paste("Gradient ascent did not converge after", max_iter, "iterations\n"))
  }

  # Return estimated coefficients
  return(list(beta = beta, partial_loglik = ll, iter = it))
}
