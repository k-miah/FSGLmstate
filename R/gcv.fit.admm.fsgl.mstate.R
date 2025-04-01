## gcv.fit.admm.fsgl.mstate - R-function utilizing ADMM to fit FSGL-penalized multi-state models
##                            for beta estimation for optimal lambda with
##                            minimal general cross-validation statistic (GCV) via grid search
##
## Input: lambda.grid     [vector]: Candidate vector for overall regularization parameter in [0,1]
##        X           [data frame]: Regression matrix of dimension n x p (=P*Q)
##                                  with transition-specific covariates
##        d           [data frame]: Data set with variables Tstart, Tstop, trans and status
##                                  (long format data)
##        penalized   [data frame]: Regression matrix of dimension n x p (=P*Q)
##                                  with covariates that should be penalized
##        unpenalized [data frame]: Regression matrix of dimension n x p (=P*Q)
##                                  with additional covariates that should remain unpenalized
##        K               [matrix]: Penalty matrix of dimension M x p (=P*Q)
##        standardize      [logic]: Standardization of design matrix X
##                                  (TRUE: columns divided by standard deviation)

##        nl             [numeric]: Number of rows of K that encode the lasso penalty
##                                  (If lasso penalty is applied to all coefficients: p)
##        nf             [numeric]: Number of rows of K that encode the fused penalty
##        ng             [numeric]: Number of groups for the group penalty
##        groupsizes      [vector]: Vector of length ngroups that gives the size of each group
##                                  in the order they appear in the K matrix (Sum should equal ng)

##        penalty.factor  [vector]: Individual penalty scaling factor (default: 1)
##        alpha.grid      [vector]: Tuning parameter in [0,1]; controls degree of
##                                  group (alpha = 0) vs lasso (alpha=1) penalty
##        gamma.grid      [vector]: Tuning parameter in [0,1]; controls degree of
##                                  lasso (gamma=1) vs fused (gamma=0) penalty

##        rho            [numeric]: Augmented Lagrangian parameter (ADMM step size; default: 1)
##        beta.init       [vector]: Initial value of beta (default: 0)

##        step_size      [numeric]: Gradient ascent step size in (0,1) (default: .01)
##        est_tol        [numeric]: Tolerance of stopping criterion (partial log-likelihood)
##                                  for beta estimation (default: 1e-6)

##        eps_rel        [numeric]: Relative tolerance for ADMM stopping criterion (default: .01)
##        eps_abs        [numeric]: Absolute tolerance for ADMM stopping criterion (default: .0001)
##        max_iter       [numeric]: Maximum number of iterations (default: 1000)
##        n.cores        [numeric]: Number of cores for parallel computing
##
## Output: res.min.gcv      [list]: Beta estimation for optimal lambda (i.e. minimal GCV)

gcv.fit.admm.fsgl.mstate <- function(lambda.grid, X, d,
                                     penalized = NULL, unpenalized = NULL,
                                     K, standardize = TRUE,
                                     nl, nf, ng, groupsizes, penalty.factor = 1,
                                     alpha.grid = seq(0, 1, by = 0.25),
                                     gamma.grid = seq(0, 1, by = 0.25),
                                     rho = 1, beta.init = NULL,
                                     step_size = 0.01, est_tol = 1e-6,
                                     eps_rel = 1e-2, eps_abs = 1e-4,
                                     max_iter = 100, n.cores = 1){

  # all combinations of tuning parameters along grids
  alpha_gamma_lambda <- expand.grid(alpha = alpha.grid, gamma = gamma.grid,
                                    lambda = lambda.grid)

  #  res.fsgl.mstate_all <- lapply(seq_len(nrow(alpha_gamma_lambda)), function(i, ...){
  res.fsgl.mstate_all <- mclapply(seq_len(nrow(alpha_gamma_lambda)), function(i, ...){

    alpha <- alpha_gamma_lambda$alpha[i]
    gamma <- alpha_gamma_lambda$gamma[i]
    lambda <- alpha_gamma_lambda$lambda[i]

    fit.admm.fsgl.mstate(X = X, penalized = penalized, unpenalized = unpenalized,
                         d = d, K = K, nl = nl, nf = nf, ng = ng,
                         groupsizes = groupsizes, penalty.factor = penalty.factor,
                         standardize = standardize,
                         lambda = lambda, beta.init = beta.init,
                         alpha = alpha, gamma = gamma, rho = rho,
                         step_size = step_size, est_tol = est_tol,
                         eps_rel = eps_rel, eps_abs = eps_abs,
                         max_iter = max_iter)
  }, mc.cores = n.cores)
  #   })

  gcv <- sapply(res.fsgl.mstate_all, function(x) x$gcv)

  # Result for optimal lambda:
  index.min.gcv <- which.min(gcv)
  res.min.gcv <- res.fsgl.mstate_all[[index.min.gcv]]
  #lambda.min <- lambda.grid[index.min.gcv]

  return(list(res.min.gcv = res.min.gcv, res.all = res.fsgl.mstate_all))
}
