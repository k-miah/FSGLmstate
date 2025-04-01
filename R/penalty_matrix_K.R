## penalty_matrix_K - R-function constructing a combined penalty structure matrix
##                    for use in penalized regression
##
## Input:  P                 [numeric]: Number of regression parameters, i.e. covariates
##         Q                 [numeric]: Number of transitions
##         fused [character or matrix]: Character string/matrix indicating whether
##                                      - all pairwise differences ("all") or
##                                      - adjacent differences ("neighbors") or
##                                      - user-specific pairs (matrix)
##                                      shall be penalized for the fused penalty
##         D                 [matrix]: Difference matrix
##       groups             [vector]: Vector indicating group membership of regression parameters
##                                    for the group penalty
##
## Output: K                [matrix]: General combined penalty matrix of
##                                    dimension (P*Q + s + P*Q) x (P*Q)

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
