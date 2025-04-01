## S_kappa - R-function implementing the vector soft thresholding operator S_kappa(a)

## Input: a      [vector]: Numeric vector
##        kappa [numeric]: Scalar

## Output:    s [numeric]: Shrinked value

S_kappa <- function(a, kappa){
  a <- as.matrix(a)
  if(all(a==0) | any(is.na(a))){
    s <- 0
  } else {
    s <- max(0, 1-kappa/norm(a, type="F"))*a
  }
  return(s)
}
