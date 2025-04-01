#' calculate Fisher information matrix of partial log-likelihood
#'
#' R-function calculating Fisher information matrix of partial log-likelihood
#'
#' @param X [matrix]: Regression matrix of dimension n_obs x p_vars
#' @param d [data frame]: Data set with variables Tstart, Tstop, trans and status
#' @param beta [vector]: Regression parameter
#' @param Riskset [list]: Risk set list
#'
#' @returns info: Fisher information matrix at beta
#' @export
#'
#' @examples

cox_ll_fisher <- function(X, d, beta, Riskset){

  X <- as.matrix(X)
  event <- d$status
  n <- length(event)

  P <- length(beta)
  f <- as.numeric(X %*% beta)
  ef <- exp(f)

  info <- matrix(nrow = P, ncol = P, 0)
  index <- which(event == 1)

  for (p in 1:P) {
    for (q in 1:P) {
      part1 <- part2 <- rep(0, n)
      for (i in index) {
        j <- Riskset[[i]]
        ef.j <- ef[j]
        risk <- sum(ef.j)
        X.j.p <- X[j, p]
        X.j.q <- X[j, q]
        part1[i] <- sum(ef.j * X.j.p * X.j.q)/risk
        part2[i] <- sum(ef.j * X.j.p) * sum(ef.j * X.j.q)/(risk * risk)
      }
      info[p, q] <- sum(event * part1) - sum(event * part2)
    }
  }
  return(info)
}
