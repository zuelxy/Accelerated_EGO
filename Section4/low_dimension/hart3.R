library(DiceOptim)
library(randtoolbox)

source("K:/Accelerated_EGO_code/section_4/low_dimension/algorithms_code/AEGO_eps.R")
source("K:/Accelerated_EGO_code/section_4/low_dimension/algorithms_code/qEGO_eps.R")
source("K:/Accelerated_EGO_code/section_4/low_dimension/algorithms_code/EGO_eps.R")
source("K:/Accelerated_EGO_code/affine.R")

hart3 <- function(xx)
{
  # xx be a 3-dimensional vector
  alpha <- c(1.0, 1.2, 3.0, 3.2)
  A <- c(3.0, 10, 30,
         0.1, 10, 35,
         3.0, 10, 30,
         0.1, 10, 35)
  A <- matrix(A, 4, 3, byrow = TRUE)
  P <- 10^(-4) * c(3689, 1170, 2673,
                   4699, 4387, 7470,
                   1091, 8732, 5547,
                   381, 5743, 8828)
  P <- matrix(P, 4, 3, byrow = TRUE)
  
  xxmat <- matrix(rep(xx, times = 4), 4, 3, byrow = TRUE)
  inner <- rowSums(A[,1:3]*(xxmat-P[,1:3])^2)
  outer <- sum(alpha * exp(-inner))
  
  y <- -outer
  return(y)
}
M <- hart3(c(0.114614, 0.555649, 0.852547))  # global optimum

U <- read.table("K:/Accelerated_EGO_code/section_4/low_dimension/UD35.txt", header = T)

lower <- rep(0, 3)
upper <- rep(1, 3)

X <- experiments(U, lower, upper)
response_X <- apply(X, 1, hart3)

fitted_model <- km(design = X, response = response_X, optim.method = "gen",
                       control=list(trace = FALSE))

q <- 4             # updating points each stage, set to be 4 ,8 and 12 respectively  

recept <- matrix(ncol = 2, nrow = 100)   # store the number of iterations satisfied the stop criterion
CPU <- matrix(ncol = 2, nrow = 100)      # store the CPU time 

for (k in 1:100) {
  
################################# A_EGO #############################
QMC <- sobol(150, dim = 3, scrambling = 3)
QMC <- affine(QMC, lower, upper)
begin_AEGO <- proc.time()
result <- AEGO_eps(QMC = QMC, fun = hart3, model = fitted_model, npoints = q, 
                    optimum = M, eps = 10^(-4), lower = lower, upper = upper,
                    control = list(pop.size = 60, max.generations = 60, 
                                   wait.generations = 6, BFGSburnin = 6))
    
end_AEGO <- proc.time()-begin_AEGO
################################ constant liar ########################
begin_CL <- proc.time()
result_CL <- qEGO_eps(fun = hart3, model = fitted_model, npoints = q, optimum = M, 
                          eps = 10^(-4), lower = lower, upper = upper, crit = "CL",
                          optimcontrol = list(L = "min", pop.size = 60, max.generations = 60, 
                                              wait.generations = 6, BFGSburnin = 6))
end_CL <- proc.time()-begin_CL
    
recept[k, ] <- c(result$nsteps, result_CL$nsteps) 
CPU[k, ] <- c(end_AEGO[1], end_CL[1])
}

##################################### EGO ##############################
repet_EGO <- c()
CPU_EGO <- c()
for (k in 1:100) {
  begin_EGO <- proc.time()
  result_EGO <- EGO_eps(model = fitted_model, fun = hart3, optimum = M, eps = 10^(-4),
                        lower = lower, upper = upper, 
                        control = list(pop.size = 60, max.generations = 60, 
                                       wait.generations = 6, BFGSburnin = 6))
  end_EGO <- proc.time()-begin_EGO
  repet_EGO[k] <- result_EGO$nsteps
  CPU_EGO[k] <- end_EGO[1]
}

