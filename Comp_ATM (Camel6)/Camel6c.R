library(DiceOptim)
source("affine.R")

camel6c <- function(xx){
  runsum <- 0
  ntimes <- length(xx)/2
  for (ll in 1:ntimes) {
    x1 <- xx[2 * ll - 1] 
    x2 <- xx[2 * ll]
    term1 <- (4 - 2.1 * x1^2 + (x1^4)/3) * x1^2
    term2 <- x1 * x2
    term3 <- (-4 + 4 * x2^2) * x2^2
    runsum <- runsum + term1 + term2 + term3
  }
  return(runsum)
}

UD <- read.table("UD41_8.txt", header = T)

lower <- c(-2,-1,-2,-1,-2,-1,-2,-1)
upper <- c(2,1,2,1,2,1,2,1)

design <- experiments(UD, lower, upper)
output <- apply(design,1,camel6c)

set.seed(124)
fitted_model <- km(design =  design, response = output)

nsteps <- 30

oEGO <- EGO.nsteps(model = fitted_model,fun = camel6c, nsteps = nsteps, 
                   lower = lower, upper = upper, control=list(pop.size=120,
                                                              max.generations=120, wait.generations=12, BFGSburnin=12))
min(oEGO$value)

########### Accelerated EGO
source("randomization.R")
source("AEGO_nsteps.R")

set.seed(999)
QMC <- sobol(750, dim = 8, scrambling = 3)
QMC <- affine(QMC, lower, upper)


result <- AEGO_nsteps(QMC = QMC, fun = camel6c, model = fitted_model, npoints = 5, 
                      nsteps = 6, lower = lower, upper = upper,
                      control = list(pop.size=120,max.generations=120, 
                                     wait.generations=12, BFGSburnin=12))


