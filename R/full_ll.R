# Negative full log-likelihood function:

## Input: X             [matrix]: Regression matrix of dimension n_obs x p_vars
##        d         [data frame]: Data set with variables Tstart, Tstop, trans and status
##        beta          [vector]: Regression parameter
##        Riskset       [matrix]: Risk set matrix
##
## Output:      loglik [numeric]: Negative full log-likelihood at beta

full_ll <- function(X, d, beta, Riskset){

  # Data extraction
  status <- d$status

  # Linear predictor
  lp <- X %*% beta
  ws <- drop(exp(lp))

  # Breslow estimate of baseline hazard: lambda_0(y_i)
  breslows <- drop(1 / ws %*% Riskset)
  # Weighted sum of Breslow estimates: Lambda_0(y_i)
  breslow <- drop(Riskset %*% breslows)

  # Full log-likelihood
  loglik <- -sum(ws * breslow) + sum(log(breslows)) + sum(lp[status==1])

  return(loglik)
}
