#' S_kappa
#'
#' R-function implementing the vector soft thresholding operator S_kappa(a)
#'
#' @param a [vector]: Numeric vector
#' @param kappa [numeric]: Scalar
#'
#' @returns s [numeric]: Shrinked value
#' @export
#'
#' @examples

S_kappa <- function(a, kappa){
  a <- as.matrix(a)
  if(all(a==0) | any(is.na(a))){
    s <- 0
  } else {
    s <- max(0, 1-kappa/norm(a, type="F"))*a
  }
  return(s)
}
