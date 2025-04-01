## cox_newton_raphson - R-function implementing Newton-Raphson for estimation
##                      in the Cox PH model
##
## Input: X          [matrix]: Regression matrix of dimension n_obs x p_vars
##        d      [data frame]: Data set with variables Tstart, Tstop, trans and status
##        K          [matrix]: Penalty matrix of dimension M x p (M=p+m+g)

##        beta.init  [vector]: Initial value of beta (default: 0)
##        Riskset      [list]: Risk set list
##        max_iter  [numeric]: Maximum number of iterations (default: 1000)
##        tolerance [numeric]: Tolerance for stopping criterion (default: 1e-6)
##        eps       [numeric]: Step size in [0,1] (default: .01)
##        rho       [numeric]: Augmented Lagrangian parameter (step size; default: 1)
##        theta     [numeric]: ADMM parameter theta = beta
##        nu        [numeric]: ADMM parameter nu = theta - beta
##
## Output:         res [list]: Beta estimation for Cox model at stopping iteration 'iter'

cox_newton_raphson <- function(X, d, K, beta.init = NULL, Riskset,
                               max_iter = 1000, tolerance = 1e-6, eps = 0.01,
                               rho = 1, theta, nu){

  # Initialize coefficients
  if(is.null(beta.init)){
    beta <- rep(0, ncol(X))
  } else{
    beta <- beta.init
  }

  # Update until convergence or maximum iterations reached
  it <- 1
  tol <- 1
  ll <- 0

  while(tol > tolerance && it < max_iter){

    if(it == 1){
      if(is.null(K)){
        gradient <- cox_ll_lagr_gradient(X, d = d, beta = beta, Riskset = Riskset,
                                         rho = rho, theta = theta, nu = nu)
      } else{
        gradient <- cox_ll_lagr_gradient_K(X, d = d, K = K, beta = beta, Riskset = Riskset,
                                           rho = rho, theta = theta, nu = nu)
      }
      beta_new <- beta + eps * gradient
    } else{

      if(is.null(K)){
        gradient <- cox_ll_lagr_gradient(X, d = d, beta = beta, Riskset = Riskset,
                                         rho = rho, theta = theta, nu = nu)
        H <- cox_ll_fisher(X, d = d, beta = beta, Riskset = Riskset) + rho * diag(ncol(X))
      } else{
        gradient <- cox_ll_lagr_gradient_K(X, d = d, K = K, beta = beta, Riskset = Riskset,
                                           rho = rho, theta = theta, nu = nu)
        H <- cox_ll_fisher(X, d = d, beta = beta, Riskset = Riskset) + rho * t(K) %*% K
      }

      # Pseudo-inverse of Fisher matrix
      M <- svd(H)
      zero_indices <- which(M$d == 0) # Check for zero singular values in the diagonal matrix M$d
      M$d[zero_indices] <- 1          # Replace zero singular values with 1 (fulfills D^-1=0 for d[ii]=0)
      M <- M$v %*% diag(1/M$d) %*% t(M$u)

      # Update coefficients: Newton-Raphson step
      beta_new <- beta + M %*% gradient
    }

    # Stopping criterion: partial log-likelihood
    ll_old <- ll
    ll <- partial_ll(beta = beta_new, X, d, risksetlist = Riskset)
    tol <- abs(ll - ll_old)

    # Step-halving if needed
    if(it > 1 & ll_old > ll){
      beta_new <- beta/2 + beta_new/2
      ll <- ll_old
    }

    it <- it + 1
    beta <- beta_new
  }
  # If max_iter is reached, print warning and exit loop
  if(it == max_iter){
    warning(paste("Newton-Raphson did not converge after", max_iter, "iterations\n"))
  }

  return(list(beta = beta, partial_loglik = ll, fisher = H, iter = it))
}
