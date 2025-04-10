% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cox_newton_raphson.R
\name{cox_newton_raphson}
\alias{cox_newton_raphson}
\title{Performing Newton-Raphson optimization algorithm for estimation in the Cox PH model}
\usage{
cox_newton_raphson(
  X,
  d,
  K,
  beta.init = NULL,
  Riskset,
  max_iter = 1000,
  tolerance = 1e-06,
  eps = 0.01,
  rho = 1,
  theta,
  nu
)
}
\arguments{
\item{X}{\link{matrix}: Regression matrix of dimension n_obs x p_vars}

\item{d}{\link{data frame}: Data set with variables Tstart, Tstop, trans and status}

\item{K}{\link{matrix}: Penalty matrix of dimension M x p}

\item{beta.init}{\link{vector}: Initial value of beta (default: 0)}

\item{Riskset}{\link{list}: Risk set list}

\item{max_iter}{\link{numeric}: Maximum number of iterations (default: 1000)}

\item{tolerance}{\link{numeric}: Tolerance for stopping criterion (default: 1e-6)}

\item{eps}{\link{numeric}: Step size in \link{0,1} (default: .01)}

\item{rho}{\link{numeric}: Augmented Lagrangian parameter (step size; default: 1)}

\item{theta}{\link{vector}: ADMM parameter theta = beta}

\item{nu}{\link{vector}: ADMM parameter nu = theta - beta}
}
\value{
res \link{list}: Beta estimation for Cox model at stopping iteration 'iter'
}
\description{
R-function implementing Newton-Raphson optimization for estimation in the Cox PH model
}
