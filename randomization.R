randomization <- function(QMC, lower = rep(0, ncol(QMC)), upper = rep(1, ncol(QMC)))
{
  ######################################################
  # This function is to randomize the QMC points using analogy modular operation
  # 'QMC' is a quasi Monte Carlo point pool ,such as Sobol sequence
  # 'lower' and 'upper' are the boundary of the domain
  ######################################################
  
  n <- nrow(QMC)
  p <- ncol(QMC)
  uu <- rep(NA, p)
  for (i in 1:p) {
    uu[i] <- runif(1, lower[i], upper[i])
  }
  QMC <- sweep(QMC, 2, uu, "+")
  lowerm <- matrix(rep(lower,rep(n,p)),n,p)
  upperm <- matrix(rep(upper,rep(n,p)),n,p)
  right <- QMC - upperm
  QMC[right>0] <- lowerm[right>0] + right[right>0]
  left <- lowerm - QMC
  QMC[left>0] <- upperm[left>0] - left[left>0]

  return(QMC)
}


resampling <- function(k, RQMC, model, minimization = minimization){
  EI_value <- apply(RQMC, 1, EI, model, minimization = minimization)
  weight <- EI_value/sum(EI_value)
  index <- sample(1:nrow(RQMC), k, prob = weight)
  add_points <- RQMC[index,]
  return(add_points)
}






