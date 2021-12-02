library(DiceOptim)
library(randtoolbox)

source("K:/Accelerated_EGO_code/section_4/high_dimension/algorithms_code/AEGO_nsteps.R")
source("K:/Accelerated_EGO_code/section_4/high_dimension/algorithms_code/qEGO_nteps.R")
source("K:/Accelerated_EGO_code/section_4/high_dimension/algorithms_code/EGO_nsteps.R")
source("K:/Accelerated_EGO_code/affine.R")

ackley <- function(xx, a=20, b=0.2, c=2*pi)
{
  
  d <- length(xx)
  
  sum1 <- sum(xx^2)
  sum2 <- sum(cos(c*xx))
  
  term1 <- -a * exp(-b*sqrt(sum1/d))
  term2 <- -exp(sum2/d)
  
  y <- term1 + term2 + a + exp(1)
  return(y)
}

M <- 0

U <- read.table("K:/Accelerated_EGO_code/section_4/high_dimension/UD10_100.txt", header = T)

p <- ncol(U)

lower <- rep(-2, p)
upper <- rep(2, p)

q <- 10              # other case is set to be 15
nsteps <- 15         # other case is set to be 10

X <- experiments(U, lower, upper)
response_X <- apply(X, 1, ackley)

fitted_model <- km(design = X, response = response_X, optim.method = "gen",
                   control=list(trace = FALSE))

recept1 <- matrix(ncol = 10, nrow = 15)
CPU1 <- c()                                    # store the vlaues and CPU time for AEGO 
recept_CL <- matrix(ncol = 10, nrow = 15)
CPU_CL <- c()                                  # store the vlaues and CPU time for CL(min)  
recept_EGO <- matrix(ncol = 10, nrow = 150)
CPU_EGO <- c()                                 # store the vlaues and CPU time for EGO

for(k in 1:10){             # repeat 10 times  
  
################################# A_EGO #############################
  
QMC <- sobol(750, dim = p, scrambling = 3)
QMC <- affine(QMC, lower, upper)

begin_AEGO <- proc.time()
result <- AEGO_nsteps(QMC = QMC, fun = ackley, model = fitted_model, npoints = q, 
                        nsteps = nsteps, lower = lower, upper = upper,
                        control = list(max.generations = 80, 
                                       wait.generations = 8, BFGSburnin = 8))
end_AEGO <- proc.time()-begin_AEGO
  
recept1[, k] <- result$history
CPU1[k] <- end_AEGO[1]
  
############################# Constant Liar ##########################
begin_CL <- proc.time()
result_CL <- qEGO_nsteps(fun = ackley, model = fitted_model, npoints = q, nsteps = nsteps,
                           lower = lower, upper = upper, crit = "CL",
                           optimcontrol = list(L = "min", max.generations = 80, 
                                               wait.generations = 8, BFGSburnin = 8))
end_CL <- proc.time()-begin_CL
  
recept_CL[, k] <- result_CL$history
CPU_CL[k] <- end_CL[1]
  
##################################### EGO ##############################
nsteps_e <- 150
  
begin_EGO <- proc.time()
result_EGO <- EGO_nsteps(model = fitted_model, fun = ackley, nsteps = nsteps_e,
                           lower = lower, upper = upper, 
                           control = list(max.generations = 80, 
                                          wait.generations = 8, BFGSburnin = 8))
end_EGO <- proc.time()-begin_EGO
  
recept_EGO[, k] <- result_EGO$history
CPU_EGO[k] <- end_EGO[1]
  
}
