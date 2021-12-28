library(LHD)
source("random_LHD.R")
source("affine.R")
source("randomization.R")

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

M <- 0   # global optimum

U <- read.table("UD21.txt", header = T) # initial uniform design

lower <- rep(-2, 2)
upper <- rep(2, 2)

X <- experiments(U, lower, upper)
response_X <- apply(X, 1, ackley)

set.seed(1)
fitted_model <- km(design = X, response = response_X, optim.method = "gen",
                   control = list(trace = FALSE))

#@#######$$$$$$$$$$$ m=50
m=50
set.seed(13)
pool <- c()
for (k in 1:100) {
  result_LHD <- rLHD_EGO(fun = ackley, model = fitted_model, npoints = 5, 
                         optimum = M, eps = 10^(-1), lower = lower, upper = upper, m=m, 
                         control = list(pop.size = 50, max.generations = 60, 
                                        wait.generations = 5, BFGSburnin = 5))
  pool <- c(pool, result_LHD$nsteps)
}

#@####### m=100
m=100
set.seed(51)
pool_100 <- c()
for (k in 1:100) {
  result_LHD <- rLHD_EGO(fun = ackley, model = fitted_model, npoints = 5, 
                         optimum = M, eps = 10^(-1), lower = lower, upper = upper, m=m, 
                         control = list(pop.size = 50, max.generations = 60, 
                                        wait.generations = 5, BFGSburnin = 5))
  pool_100 <- c(pool_100, result_LHD$nsteps)
}

#@####### m=150
m=150
set.seed(44)
pool_150 <- c()
for (k in 1:100) {
  result_LHD <- rLHD_EGO(fun = ackley, model = fitted_model, npoints = 5, 
                         optimum = M, eps = 10^(-1), lower = lower, upper = upper, m=m, 
                         control = list(pop.size = 50, max.generations = 60, 
                                        wait.generations = 5, BFGSburnin = 5))
  pool_150 <- c(pool_150, result_LHD$nsteps)
}
