library(DiceOptim)
library(randtoolbox)

source("K:/Accelerated_EGO_code/section_4/low_dimension/algorithms_code/AEGO_eps.R")
source("K:/Accelerated_EGO_code/section_4/low_dimension/algorithms_code/qEGO_eps.R")
source("K:/Accelerated_EGO_code/section_4/low_dimension/algorithms_code/EGO_eps.R")
source("K:/Accelerated_EGO_code/affine.R")

sin2 <- function(xx){
  x1 <- xx[1]
  x2 <- xx[2]
  y <- 1+(sin(x1))^2+(sin(x2))^2-0.1*exp(-x1^2-x2^2)
  return(y)
}

M <- sin2(c(0,0))   # global optimum

U <- read.table("K:/Accelerated_EGO_code/section_4/low_dimension/UD21.txt", header = T)

lower <- c(-5, -5)
upper <- c(5, 5)

X <- experiments(U, lower, upper)
response_X <- apply(X, 1, sin2)

fitted_model <- km(design = X, response = response_X, optim.method = "gen",
                   control=list(trace = FALSE))

q <- 4             # updating points each stage, set to be 4 ,8 and 12 respectively  

recept <- matrix(ncol = 2, nrow = 100)   # store the number of iterations satisfied the stop criterion
CPU <- matrix(ncol = 2, nrow = 100)      # store the CPU time 

for (k in 1:100) {            # repeat 100 times 
  
  ################################# A_EGO #############################
  QMC <- sobol(100, dim = 2, scrambling = 3)
  QMC <- affine(QMC, lower, upper)
  begin_AEGO <- proc.time()
  result <- AEGO_eps(QMC = QMC, fun = sin2, model = fitted_model, npoints = q, 
                  optimum = M, eps = 10^(-2), lower = lower, upper = upper,
                  control = list(pop.size = 50, max.generations = 50, 
                                 wait.generations = 5, BFGSburnin = 5))
  
  end_AEGO <- proc.time()-begin_AEGO
  ################################ constant liar ########################
  begin_CL <- proc.time()
  result_CL <- qEGO_eps(fun = sin2, model = fitted_model, npoints = q, optimum = M, 
                        eps = 10^(-2), lower = lower, upper = upper, crit = "CL",
                        optimcontrol = list(L = "min", pop.size = 50, max.generations = 50, 
                                            wait.generations = 5, BFGSburnin = 5))
  end_CL <- proc.time()-begin_CL
  
  recept[k, ] <- c(result$nsteps, result_CL$nsteps) 
  CPU[k, ] <- c(end_AEGO[1], end_CL[1])
}

##################################### EGO ##############################
repet_EGO <- c()
CPU_EGO <- c()
for (k in 1:100) {
  begin_EGO <- proc.time()
  result_EGO <- EGO_eps(model = fitted_model, fun = sin2, optimum = M, eps = 10^(-2),
                        lower = lower, upper = upper, 
                        control = list(pop.size = 50, max.generations = 50, 
                                       wait.generations = 5, BFGSburnin = 5))
  end_EGO <- proc.time()-begin_EGO
  repet_EGO[k] <- result_EGO$nsteps
  CPU_EGO[k] <- end_EGO[1]
}

