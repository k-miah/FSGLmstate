# Partial log-likelihood function:

## Input: X             [matrix]: Regression matrix of dimension n_obs x p_vars
##        d         [data frame]: Data set with variables Tstart, Tstop, trans and status
##        beta          [vector]: Regression parameter
##        risksetlist     [list]: Risk set list
##
## Output:     logplik [numeric]: Partial log-likelihood at beta

partial_ll <- function(X, d, beta, risksetlist){

  X <- as.matrix(X)
  event <- d$status
  n <- length(event)
  risk <- numeric(n)
  f <- as.numeric(X %*% beta)
  ef <- exp(f)

  for(i in which(event == 1)){
    risk[i] <- sum(ef[risksetlist[[i]]])
  }
  logplik <- sum(event * (f - log(risk)), na.rm = TRUE)
  return(logplik)
}
