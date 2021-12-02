library(DiceOptim)
library(randtoolbox)

source("K:/Accelerated_EGO_code/section_4/low_dimension/algorithms_code/AEGO_eps.R")
source("K:/Accelerated_EGO_code/section_4/low_dimension/algorithms_code/qEGO_eps.R")
source("K:/Accelerated_EGO_code/section_4/low_dimension/algorithms_code/EGO_eps.R")
source("K:/Accelerated_EGO_code/affine.R")

GoldPrice <- function(xx)
{
  x1 <- xx[1]
  x2 <- xx[2]
  
  fact1a <- (x1+x2+1)^2
  fact1b <- 19-14*x1+3*x1^2-14*x2+6*x1*x2+3*x2^2
  fact1 <- 1+fact1a*fact1b
  
  fact2a <- (2*x1-3*x2)^2
  fact2b <- 18-32*x1+12*x1^2+48*x2-36*x1*x2+27*x2^2
  fact2 <- 30+fact2a*fact2b
  
  scprod <- fact1*fact2
  
  y <- (log(scprod)-8.693)/2.427
  return(y)
}

M <- GoldPrice(c(0,-1))     # global optimum

U <- read.table("K:/Accelerated_EGO_code/section_4/low_dimension/UD21.txt", header = T)

p <- ncol(U)
lower <- rep(-2, p)
upper <- rep(2, p)

X <- experiments(U, lower, upper)
response_X <- apply(X, 1, GoldPrice)

fitted_model <- km(design = X, response = response_X, optim.method = "gen",
                   control=list(trace = FALSE))

q <- 4             # updating points each stage, set to be 4 ,8 and 12 respectively  

recept <- matrix(ncol = 2, nrow = 100)   # store the number of iterations satisfied the stop criterion
CPU <- matrix(ncol = 2, nrow = 100)      # store the CPU time 

for (k in 1:100) {            # repeat 100 times
  
  ################################# A_EGO #############################
  QMC <- sobol(120, dim = p, scrambling = 3)
  QMC <- affine(QMC, lower, upper)
  begin_AEGO <- proc.time()
  result <- AEGO_eps(QMC = QMC, fun = GoldPrice, model = fitted_model, npoints = q, 
                     optimum = M, eps = 10^(-2), lower = lower, upper = upper,
                     control = list(max.generations = 70, 
                                    wait.generations = 7, BFGSburnin = 7))
  
  end_AEGO <- proc.time()-begin_AEGO
  ################################ constant liar ########################
  begin_CL <- proc.time()
  result_CL <- qEGO_eps(fun = GoldPrice, model = fitted_model, npoints = q, optimum = M, 
                        eps = 10^(-2), lower = lower, upper = upper, crit = "CL",
                        optimcontrol = list(L = "min", max.generations = 70, 
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
  result_EGO <- EGO_eps(model = fitted_model, fun = GoldPrice, optimum = M, eps = 10^(-2),
                        lower = lower, upper = upper, 
                        control = list(max.generations = 70, 
                                       wait.generations = 7, BFGSburnin = 7))
  end_EGO <- proc.time()-begin_EGO
  repet_EGO[k] <- result_EGO$nsteps
  CPU_EGO[k] <- end_EGO[1]
}


