library(DiceOptim)
library(randtoolbox)

source("K:/Accelerated_EGO_code/section_4/low_dimension/algorithms_code/AEGO_eps.R")
source("K:/Accelerated_EGO_code/section_4/low_dimension/algorithms_code/qEGO_eps.R")
source("K:/Accelerated_EGO_code/section_4/low_dimension/algorithms_code/EGO_eps.R")
source("K:/Accelerated_EGO_code/affine.R")


bra <- function(xx, a=1, b=5.1/(4*pi^2), c=5/pi, r=6, s=10, t=1/(8*pi))
{
  x1 <- xx[1]
  x2 <- xx[2]
  
  term1 <- a * (x2 - b*x1^2 + c*x1 - r)^2
  term2 <- s*(1-t)*cos(x1)
  
  y <- term1 + term2 + s
  return(y)
}


M <- bra(c(pi, 2.275)) # global optimum

U <- read.table("K:/Accelerated_EGO_code/section_4/low_dimension/UD21.txt", header = T)

lower <- c(-5, 0)
upper <- c(10, 15)

X <- experiments(U, lower, upper)
response_X <- apply(X, 1, bra)

fitted_model <- km(design = X, response = response_X, optim.method = "gen",
                   control=list(trace = FALSE))

q <- 4           # updating points each stage, set to be 4 ,8 and 12 respectively

recept <- matrix(ncol = 2, nrow = 100)   # store the number of iterations satisfied the stop criterion
CPU <- matrix(ncol = 2, nrow = 100)      # store the CPU time 

for (k in 1:100) {           # the simulations are repeated 100 times
  
################################# A_EGO #############################
QMC <- sobol(100, dim = 2, scrambling = 3)
QMC <- affine(QMC, lower, upper)
begin_AEGO <- proc.time()
result <- AEGO_eps(QMC = QMC, fun = bra, model = fitted_model, npoints = q, 
                optimum = M, eps = 10^(-2), lower = lower, upper = upper,
                control = list(pop.size = 50, max.generations = 50, 
                               wait.generations = 5, BFGSburnin = 5))

end_AEGO <- proc.time()-begin_AEGO

################################ constant liar ########################
begin_CL <- proc.time()
result_CL <- qEGO_eps(fun = bra, model = fitted_model, npoints = q, optimum = M, 
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
  result_EGO <- EGO_eps(model = fitted_model, fun = bra, optimum = M, eps = 10^(-2),
                      lower = lower, upper = upper, 
                      control = list(pop.size = 50, max.generations = 50, 
                                     wait.generations = 5, BFGSburnin = 5))
  end_EGO <- proc.time()-begin_EGO
  repet_EGO[k] <- result_EGO$nsteps
  CPU_EGO[k] <- end_EGO[1]
}

