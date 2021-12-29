setwd("F:/xiaoyao/Paviani")
library(DiceOptim)
library(randtoolbox)
source("affine.R")

Paviani <- function(xx){
  a <- (log(xx))^2+(log(11-xx))^2
  b <- (prod(xx))^0.2
  return(sum(a)-b)
}

lower <- rep(1,5)
upper <- rep(10,5)

axes_grid <- 1:10

design.grid <- expand.grid(axes_grid, axes_grid, 
                           axes_grid, axes_grid, axes_grid)
response.grid <- apply(design.grid,1,Paviani)
max(response.grid)

################ 
UD <- read.table("UD50_5.txt", header = T)
design <- experiments(UD, lower, upper)
output <- apply(design,1,Paviani)

set.seed(42)
fitted_model <- km(design =  design, response = output, nugget = 1e-5)

################ EGO
set.seed(23)
EGO_result <- c()
for (s in 1:100) {
  
oEGO <- EGO_nsteps(model = fitted_model,fun = Paviani, nsteps = 100,  minimization = FALSE,
                   lower = lower, upper = upper, control=list(pop.size=120,
                                                              max.generations=120, wait.generations=12, BFGSburnin=12))
EGO_result <- c(EGO_result, max(oEGO$value))
}

########### Accelerated EGO
source("randomization.R")
source("AEGO_nsteps.R")

set.seed(699)

QMC <- sobol(500, dim = 5, scrambling = 3)
QMC <- affine(QMC, lower, upper)

AEGO_result <- c()
for (ss in 1:100) {
result <- AEGO_nsteps(QMC = QMC, fun = Paviani, model = fitted_model, npoints = 10, 
                      nsteps = 10, lower = lower, upper = upper, minimization = FALSE,
                      control = list(pop.size=120,max.generations=120, 
                                     wait.generations=12, BFGSburnin=12))
AEGO_result <- c(AEGO_result, max(result$history))
}

