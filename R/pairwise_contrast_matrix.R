#' Construct fused penalty structure matrix
#'
#' R-function constructing the fused penalty structure matrix of all pairwise differences for use in penalized regression
#'
#' @param n.cols [numeric]: Number of regression parameters, i.e. covariates
#'
#' @returns K [matrix]: Fused penalty structure matrix of dimension (2^(P*Q)) x (P*Q)
#' @export
#'
#' @examples

pairwise_contrast_matrix <- function(n.cols){
  # Initialize an empty matrix to store the combinations
  K <- matrix(0, nrow = n.cols * (n.cols - 1) / 2, ncol = n.cols)

  # Generate all pairwise combinations
  index <- 1
  for (i in 1:(n.cols - 1)) {
    for (j in (i + 1):n.cols) {
      K[index, i] <- 1
      K[index, j] <- -1
      index <- index + 1
    }
  }
  return(K)
}
