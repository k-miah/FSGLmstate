#' fit.admm.fsgl.mstate
#'
#' R-function utilizing ADMM for FSGL-penalized multi-state models for estimation of beta for one set of tuning parameters
#'
#' @param X [data frame]: Regression matrix of dimension n x p (=P*Q) with transition-specific covariates
#' @param d [data frame]: Data set with variables Tstart, Tstop, trans and status (long format data)
#' @param penalized [data frame]: Regression matrix of dimension n x p (=P*Q) with covariates that should be penalized
#' @param unpenalized [data frame]: Regression matrix of dimension n x p (=P*Q) with additional covariates that should remain unpenalized
#' @param K [matrix]: Penalty matrix of dimension M x p (=P*Q)
#' @param standardize [logic]: Standardization of design matrix X (TRUE: columns divided by standard deviation)
#' @param trace [logic]: Storage of updates/history at iteration k
#' @param nl [numeric]: Number of rows of K that encode the lasso penalty (If lasso penalty is applied to all coefficients: p)
#' @param nf [numeric]: Number of rows of K that encode the fused penalty
#' @param ng [numeric]: Number of groups for the group penalty
#' @param groupsizes [vector]: Vector of length ngroups that gives the size of each group in the order they appear in the K matrix (Sum should equal ng)
#' @param penalty.factor [vector]: Individual penalty scaling factor (default: 1)
#' @param lambda [numeric]: Overall tuning parameter in \[0,1\]
#' @param alpha [numeric]: Tuning parameter in \[0,1\]; controls degree of group (alpha = 0) vs lasso (alpha=1) penalty
#' @param gamma [numeric]: Tuning parameter in \[0,1\]; controls degree of lasso (gamma=1) vs fused (gamma=0) penalty
#' @param rho [numeric]: Augmented Lagrangian parameter (ADMM step size; default: 1)
#' @param beta.init [vector]: Initial value of beta (default: 0)
#' @param est_algorithm [character]: Cox estimation algorithm (default:'gradient.ascent')
#' @param step_size [numeric]: Cox estimation step size in (0,1) (default: .01)
#' @param est_tol [numeric]: Tolerance of stopping criterion (partial log-likelihood) for beta estimation (default: 1e-6)
#' @param eps_rel [numeric]: Relative tolerance for ADMM stopping criterion (default: .01)
#' @param eps_abs [numeric]: Absolute tolerance for ADMM stopping criterion (default: .0001)
#' @param max_iter[numeric]: Maximum number of iterations (default: 1000)
#'
#' @returns res [list]: Beta estimation at stopping iteration 'num.iter' with history
#' @export
#'
#' @examples

fit.admm.fsgl.mstate <- function(X, d, penalized = NULL, unpenalized = NULL, K,
                                 standardize = FALSE, trace = TRUE,
                                 nl, nf, ng, groupsizes, penalty.factor = 1,
                                 lambda = 1, alpha, gamma,
                                 rho = 1, beta.init = NULL,
                                 est_algorithm = "gradient.ascent",
                                 step_size = 0.01, est_tol = 1e-6,
                                 eps_rel = 1e-2, eps_abs = 1e-4,
                                 max_iter = 1000){

  n <- nrow(X)
  pq <- ncol(K)
  M <- nrow(K)

  # Total number of samples
  Ns <- nl + nf + ng

  # Standardize only continuous variables
  if(standardize){
    continuous_cols <- apply(X, 2, function(col) length(unique(col)) > 10)
    # Assume matrix columns with <= 10 unique values are categorical

    tmp <- scale(X, center = FALSE, scale = TRUE)
    tmp2 <- attributes(tmp)$`scaled:scale`
    #tmp2[!continuous_cols] <- 0
    scales <- tmp2[continuous_cols]
    X[, continuous_cols] <- scale(X[, continuous_cols], center = FALSE, scale = TRUE)

    if(!is.null(unpenalized)){
      continuous_cols <- apply(unpenalized, 2, function(col) length(unique(col)) > 10)
      tmp <- scale(unpenalized, center = FALSE, scale = TRUE)
      tmp2 <- attributes(tmp)$`scaled:scale`
      scale_unpen <- tmp2[continuous_cols]
      unpenalized[, continuous_cols] <- scale(unpenalized[, continuous_cols],
                                              center = FALSE, scale = TRUE)
    }
  }

  # Risk set list
  r <- buildrisksets(d$Tstart, d$Tstop, d$trans, d$status)
  riskset <- r$Ri

  # Calculate the cumulative sum of samples for each part of the design matrix
  ni <- c(rep(1, nl), rep(1, nf), groupsizes)
  cumsum_ni <- cumsum(ni)

  # Calculate the indices for each group based on groupsizes
  group_indices <- lapply(1:Ns, function(i) seq(max(1, cumsum_ni[i - 1] + 1), cumsum_ni[i]))

  # Overall tuning parameter vector
  Lambda <- c(rep(alpha * gamma * lambda, nl) * penalty.factor,
              rep((1 - gamma) * lambda, nf),
              rep((1 - alpha) * gamma * lambda, ng))

  # (1) Initialization:
  if(is.null(beta.init)){
    beta <- matrix(0, nrow = pq, ncol = 1)

    if(!is.null(unpenalized)){

      # Start with fully penalized model keeping only unpenalized covariates
      unpenalized.names <- colnames(unpenalized)
      long.data <- cbind(d, unpenalized)

      beta <- coefficients(coxph(as.formula(paste("Surv(Tstart, Tstop, status) ~ ",
                                                  paste(unpenalized.names, collapse = "+"),
                                                  "+ strata(trans)")),
                                 data = long.data))
      if(!is.null(penalized)){
        beta <- c(beta, rep(0, ncol(penalized)))
      }
    }
  } else{
    beta <- beta.init
  }
  theta <- matrix(0, nrow = M, ncol = 1)
  nu <- matrix(0, nrow = M, ncol = 1)


  # Orthogonalize penalized with respect to unpenalized
  if(!is.null(unpenalized) & !(is.null(penalized))){
    orthogonalizer <- solve(crossprod(unpenalized), crossprod(unpenalized, penalized))
    penalized <- penalized - unpenalized %*% orthogonalizer

    # Join penalized and unpenalized together
    X <- cbind(unpenalized, penalized)
  }
  X <- as.matrix(X)

  if(trace){
    history <- vector(mode = "list")
    updates <- vector(mode = "list")
  }
  for(k in 1:max_iter){
    # (2) Alternatingly update beta, theta, nu:

    # Update beta
    if(est_algorithm == "gradient.ascent"){
      beta_new <- cox_fixed_gradient_ascent(X = X, d = d, K = K, beta.init = beta,
                                            Riskset = riskset,
                                            eps = step_size, tolerance = est_tol,
                                            max_iter = max_iter, theta = theta, nu = nu)$beta
    } else{
      beta_new <- cox_newton_raphson(X = X, d = d, K = K, beta.init = beta,
                                     Riskset = riskset, tolerance = est_tol,
                                     max_iter = max_iter, theta = theta, nu = nu)$beta
    }

    # Calculate eta (intermediate variable) for updating theta

    eta <- K %*% beta_new + nu/rho

    # Update theta using soft-thresholding for each group
    theta_old <- theta
    theta_new <- theta
    for(i in 1:Ns){
      theta_new[group_indices[[i]]] <- S_kappa(a = eta[group_indices[[i]]],
                                               kappa = (Lambda[i] * sqrt(length(eta[group_indices[[i]]]))) / rho)
    }
    # Update nu using dual ascent
    nu_new <- nu + rho * (K %*% beta_new - theta_new)

    # Updated estimates
    beta <- as.matrix(beta_new)
    theta <- as.matrix(theta_new)
    nu <- as.matrix(nu_new)
    # residuals
    r_norm <- norm(theta - K %*% beta, type="F") # primal residuals
    s_norm <- norm(rho * t(K) %*% (theta - theta_old), type="F") # dual residuals

    # sufficiently small epsilons (Boyd et al. (2011): eps_abs = 10^-4, eps_rel = 10^-2)
    eps_pri <- sqrt(pq) * eps_abs + eps_rel * max(norm(K %*% beta, type="F"), norm(theta, type="F"))
    eps_dual <- sqrt(M) * eps_abs + eps_rel * norm(t(K) %*% nu, type="F")

    # (3) Storage at iteration k: Store updates, residuals, epsilons
    if(trace){
      updates[[k]] <- list(beta = drop(beta), theta = drop(theta), nu = drop(nu))

      history$r_norm[k] <- norm(theta - K %*% beta, type="F") # primal residuals
      history$s_norm[k] <- norm(rho * t(K) %*% (theta - theta_old), type="F") # dual residuals
      history$eps_pri[k] <- sqrt(pq) * eps_abs + eps_rel * max(norm(K %*% beta, type="F"), norm(theta, type="F"))
      history$eps_dual[k] <- sqrt(M) * eps_abs + eps_rel * norm(t(K) %*% nu, type="F")
    }

    # Adapt rho
    rho <- Adapt_rho(rho=rho, r_norm = r_norm, s_norm = s_norm)

    # (4) Stopping criterion: Sufficiently small primal & dual residuals
    if(r_norm < eps_pri && s_norm < eps_dual){
      break
    }
  }

  # Model diagnostics: Generalized cross-validation estimate (Wahba, 1980; Tibshirani, 1997)

  myTheta <- beta
  myTheta[1:nrow(beta),] <- theta[1:nrow(beta),]
  nlpl <- -partial_ll(X, d, myTheta, risksetlist = riskset)

  fisher <- fisherinfo(beta = myTheta, X = X, risksetlist = riskset, event = d$status)
  Lambda_K <- c(rep(alpha * gamma * lambda, nl),
                rep((1 - gamma) * lambda, nf),
                rep((1 - alpha) * gamma * lambda, ng*p))
  A <- penaltymatrix(lambda = Lambda_K, PSM = K, beta = myTheta, w = rep(1, M),
                     constant = 1e-08)
  M <- svd(fisher + A)
  M <- M$v %*% diag(1/M$d) %*% t(M$u)
  df <- sum(diag(fisher %*% M))

  gcv <- (1/n) * nlpl/(n * ((1 - df/n)^2))

  if(standardize){
    # Scale back estimates to original covariate scales
    beta[continuous_cols,] <- beta[continuous_cols,]/scales

    theta <- theta[1:(pq), ]
    theta <- as.matrix(theta)
    rownames(theta) <- colnames(X)
    theta[continuous_cols,] <- theta[continuous_cols,]/scales
  }
  if(!is.null(unpenalized) & !is.null(penalized)){
    # Scale back unpenalized estimates after orthogonalization
    beta[1:ncol(unpenalized), ] <- beta[1:ncol(unpenalized), ] - drop(orthogonalizer %*% beta[(ncol(unpenalized)+1):length(beta), ])
    theta[1:ncol(unpenalized), ] <- theta[1:ncol(unpenalized), ] - drop(orthogonalizer %*% theta[(ncol(unpenalized)+1):length(theta), ])

  }
  rownames(beta) <- colnames(X)
  beta <- cbind(beta, exp(beta))
  colnames(beta) <- c("Beta estimates", "exp(beta)")

  theta <- cbind(theta, exp(theta))
  colnames(theta) <- c("Theta estimates", "exp(theta)")


  res <- list(beta = beta,
              theta = theta,
              lambda = lambda,
              alpha = alpha,
              gamma = gamma,
              gcv = gcv,
              df = df,
              num.iter = k)
  if(trace){
    res$updates <- updates
    res$history <- history
  }

  return(res)
}
