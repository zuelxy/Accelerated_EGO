affine <- function(x, lower, upper){
  x_affine <- matrix(ncol = ncol(x),nrow = nrow(x))
  for (i in 1:ncol(x)) {
    x_affine[,i] <- (upper[i]-lower[i])*x[,i]+lower[i]
  }
  return(x_affine)
}

experiments <- function(U, lower, upper){
  U <- (U-0.5)/nrow(U)
  # source("K:/A_EGO/affine.R")
  # X <- affine(U, lower = lower, upper = upper)
  X <- matrix(ncol = ncol(U),nrow = nrow(U))
  for (i in 1:ncol(U)) {
    X[,i] <- (upper[i]-lower[i])*U[,i]+lower[i]
  }
  return(X)
}
