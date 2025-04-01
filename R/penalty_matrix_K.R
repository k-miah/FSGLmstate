#' Construct combined penalty structure matrix
#'
#' R-function constructing a combined penalty structure matrix for use in penalized regression
#'
#' @param P [numeric]: Number of regression parameters, i.e. covariates
#' @param Q [numeric]: Number of transitions
#' @param fused [character or matrix]: Character string/matrix indicating whether all pairwise differences ("all") or adjacent differences ("neighbors") or user-specific pairs (matrix) shall be penalized for the fused penalty
#' @param D [matrix]: Difference matrix
#' @param groups [vector]: Vector indicating group membership of regression parameters for the group penalty
#'
#' @returns K [matrix]: General combined penalty matrix of dimension (P*Q + s + P*Q) x (P*Q)
#' @export
#'
#' @examples

penalty_matrix_K <- function(P, Q, fused = "all", D = NULL, groups){

  n.col <- P*Q

  # Lasso penalty: unit matrix
  I <- diag(n.col)

  # Fused penalty: difference matrix
  if(fused == "all"){
    D <- pairwise_contrast_matrix(n.col)
  }
  if(fused == "neighbors"){
    D <- -(diff(diag(n.row), diff = 1))
  } else{
    D <- D
  }

  # Group penalty: unit vectors where 1 indicates a sample of a group
  group_matrices <- lapply(1:max(groups), function(i) diag(groups == i))
  G <- do.call(rbind, group_matrices)
  G <- G[rowSums(G) == 1, ]

  # General combined penalty matrix
  K <- rbind(I, D, G)

  return(K)
}
