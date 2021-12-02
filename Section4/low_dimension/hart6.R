library(DiceOptim)
library(randtoolbox)

source("K:/Accelerated_EGO_code/section_4/low_dimension/algorithms_code/AEGO_eps.R")
source("K:/Accelerated_EGO_code/section_4/low_dimension/algorithms_code/qEGO_eps.R")
source("K:/Accelerated_EGO_code/section_4/low_dimension/algorithms_code/EGO_eps.R")
source("K:/Accelerated_EGO_code/affine.R")

hart6 <- function(xx){
  alpha <- c(1.0, 1.2, 3.0, 3.2)
  A <- c(10, 3, 17, 3.5, 1.7, 8,
         0.05, 10, 17, 0.1, 8, 14,
         3, 3.5, 1.7, 10, 17, 8,
         17, 8, 0.05, 10, 0.1, 14)
  A <- matrix(A, 4, 6, byrow=TRUE)
  P <- 10^(-4) * c(1312, 1696, 5569, 124, 8283, 5886,
                   2329, 4135, 8307, 3736, 1004, 9991,
                   2348, 1451, 3522, 2883, 3047, 6650,
                   4047, 8828, 8732, 5743, 1091, 381)
  P <- matrix(P, 4, 6, byrow=TRUE)
  
  xxmat <- matrix(rep(xx,times=4), 4, 6, byrow=TRUE)
  inner <- rowSums(A[,1:6]*(xxmat-P[,1:6])^2)
  outer <- sum(alpha * exp(-inner))
  
  y <- -outer
  return(y)
}
M <- hart6(c(0.20169, 0.150011, 0.476874, 0.275332, 0.311652, 0.6573))  # global optimum


U <- read.table("K:/Accelerated_EGO_code/section_4/low_dimension/UD21.txt", header = T)

lower <- rep(0, 6)
upper <- rep(1, 6)

X <- experiments(U, lower, upper)
response_X <- data.frame(apply(X, 1, hart6)) 
    
fitted_model <- km(design = X, response = response_X, optim.method = "gen",
                       control=list(trace = FALSE))

q <- 4             # updating points each stage, set to be 4 ,8 and 12 respectively  

recept <- matrix(ncol = 2, nrow = 100)   # store the number of iterations satisfied the stop criterion
CPU <- matrix(ncol = 2, nrow = 100)      # store the CPU time 

for (k in 1:100) {
  
################################# A_EGO #############################
QMC <- sobol(2500, dim = 6, scrambling = 3)
QMC <- affine(QMC, lower, upper)
begin_AEGO <- proc.time()
result <- AEGO_eps(QMC = QMC, fun = hart6, model = fitted_model, npoints = q, 
                    optimum = M, eps = 10^(-1), lower = lower, upper = upper,
                    control = list(pop.size = 70, max.generations = 120, 
                                   wait.generations = 7, BFGSburnin = 7))
    
end_AEGO <- proc.time()-begin_AEGO
################################ constant liar ########################
begin_CL <- proc.time()
result_CL <- qEGO_eps(fun = hart6, model = fitted_model, npoints = q, optimum = M, 
                          eps = 10^(-1), lower = lower, upper = upper, crit = "CL",
                          optimcontrol = list(L = "min", pop.size = 70, max.generations = 120, 
                                              wait.generations = 7, BFGSburnin = 7))
end_CL <- proc.time()-begin_CL
    
recept[k, ] <- c(result$nsteps, result_CL$nsteps) 
CPU[k, ] <- c(end_AEGO[1], end_CL[1])
}

##################################### EGO ##############################
repet_EGO <- c()
CPU_EGO <- c()
for (k in 1:100) {
  begin_EGO <- proc.time()
  result_EGO <- EGO_eps(model = fitted_model, fun = hart6, optimum = M, eps = 10^(-1),
                        lower = lower, upper = upper,
                         control=list(pop.size = 70, max.generations = 120, 
                                      wait.generations = 7, BFGSburnin = 7), kmcontrol=NULL)
  end_EGO <- proc.time()-begin_EGO
  repet_EGO[k] <- result_EGO$nsteps
  CPU_EGO[k] <- end_EGO[1]
}

