#' partial_ll
#'
#' Partial log-likelihood function
#'
#' @param X [matrix]: Regression matrix of dimension n_obs x p_vars
#' @param d [data frame]: Data set with variables Tstart, Tstop, trans and status
#' @param beta [vector]: Regression parameter
#' @param Riskset [list]: Risk set list
#'
#' @returns loglik [numeric]: Partial log-likelihood at beta
#' @export
#'
#' @examples

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
