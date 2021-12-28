library(DiceOptim)
library(randtoolbox)
library(ggplot2)
library(LHD)

# setwd("~/Desktop/R_script/Accelerated_EGO_code")
source("AEGO_eg.R")
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
############################# the influence of pool size ###############################
set.seed(98)
m <- 50 ## 100, 150
pool_q_50 <- c()
QMC <- sobol(m, dim = 2, scrambling = 3)
QMC <- affine(QMC, lower, upper)

for (k in 1:100) {
result <- A_EGO(QMC = QMC, fun = ackley, model = fitted_model, npoints = 5, 
                optimum = M, eps = 10^(-1), lower = lower, upper = upper,
                control = list(pop.size = 50, max.generations = 60, 
                               wait.generations = 5, BFGSburnin = 5))
pool_q_50 <- c(pool_q_50, result$nsteps)
}

###################
set.seed(124)
m <- 100  #150
pool_100 <- c()
QMC <- sobol(m, dim = 2, scrambling = 3)
QMC <- affine(QMC, lower, upper)

for (k in 1:100) {
  
  result <- A_EGO(QMC = QMC, fun = ackley, model = fitted_model, npoints = 5, 
                  optimum = M, eps = 10^(-1), lower = lower, upper = upper,
                  control = list(pop.size = 50, max.generations = 60, 
                                 wait.generations = 5, BFGSburnin = 5))
  pool_100 <- c(pool_100, result$nsteps)
}

###################
set.seed(135)
m <- 150  #150
pool_150 <- c()
QMC <- sobol(m, dim = 2, scrambling = 3)
QMC <- affine(QMC, lower, upper)

for (k in 1:100) {
  
  result <- A_EGO(QMC = QMC, fun = ackley, model = fitted_model, npoints = 5, 
                  optimum = M, eps = 10^(-1), lower = lower, upper = upper,
                  control = list(pop.size = 50, max.generations = 60, 
                                 wait.generations = 5, BFGSburnin = 5))
  pool_150 <- c(pool_150, result$nsteps)
}


################## the iterated process of Accelerated EGO #############################
# QMC <- sobol(100, dim = 2, scrambling = 3)
# QMC <- affine(QMC, lower, upper)

# result <- A_EGO(QMC = QMC, fun = ackley, model = fitted_model, npoints = 5, 
#                optimum = M, eps = 10^(-1), lower = lower, upper = upper,
#                control = list(pop.size = 50, max.generations = 60, 
#                               wait.generations = 5, BFGSburnin = 5))

##################### initial 
par(mar = c(2, 2, 1, 1))
contour(x.grid, y.grid, z.grid)
points(0, 0, pch = 17, col = "blue")
points(X, pch = 20)

tu <- data.frame(design.grid, z = response.grid)
exp_point <- data.frame(X, z = response_X)

v <- ggplot(tu, aes(x, y, z = z)) + stat_contour(aes(colour = ..level..)) + 
  geom_point(aes(x=0,y=0),  shape = 17, col = "blue", size = 3) +
  geom_point(data = exp_point, aes(x = x1, y = x2), shape = 20, size = 3)
direct.label(v, method="bottom.pieces")

##################### step 1 
par(mar = c(2, 2, 1, 1))
contour(x.grid, y.grid, z.grid)
points(0, 0, pch = 17, col = "blue")
points(X, pch = 20)
points(-0.4656459, -0.5913320, pch = 22, bg = "green", cex = 0.8)
points(result$add[[1]], pch = 23, bg = "red", cex = 0.8)


add <- matrix(c(
  -0.2909267,  0.3562228,
  -0.7308759,  0.0777741,
  -0.2684399, -0.2681920,
  -0.5142521, -0.5710053), ncol = 2, byrow = T)
names(add) <- c("x1","x2")
exp_add1 <- data.frame(add, z = apply(add, 1, ackley))

v <- ggplot(tu, aes(x, y, z = z)) + stat_contour(aes(colour = ..level..)) + 
  geom_point(aes(x=0,y=0),  shape = 17, col = "blue", size = 3) +
  geom_point(data = exp_point, aes(x = x1, y = x2), shape = 20, size = 3) +
  geom_point(data = exp_add1, aes(x = X1, y = X2), shape = 23, bg = 'red', size = 3) +
  geom_point(aes(x=-0.4656459, y=-0.5913320),  shape = 22, bg = "orange", size = 3)
direct.label(v, method="bottom.pieces")


#################### step 2 
par(mar = c(2, 2, 1, 1))
contour(x.grid, y.grid, z.grid)
points(0, 0, pch = 17, col = "blue")
points(X, pch = 20)
points(-0.4656459, -0.5913320, pch = 20)
points(result$add[[1]], pch = 20)

points(0.05842766, -0.955879391, pch = 22, bg = "green", cex = 0.8)
points(result$add[[2]], pch = 23, bg = "red", cex = 0.8)


add <- matrix(c(
  0.01660588,  0.05130075,
  -0.05012373, -0.86955303,
  -0.15198911,  0.87749053,
  -0.31085746, -0.03720377), ncol = 2, byrow = T)
names(add) <- c("x1","x2")
exp_add2 <- data.frame(add, z = apply(add, 1, ackley))

v <- ggplot(tu, aes(x, y, z = z)) + stat_contour(aes(colour = ..level..)) + 
  geom_point(aes(x=0,y=0),  shape = 17, col = "blue", size = 3) +
  geom_point(data = exp_point, aes(x = x1, y = x2), shape = 20, size = 3) +
  geom_point(data = exp_add1, aes(x = X1, y = X2), shape = 20, size = 3) +
  geom_point(aes(x=-0.4656459, y=-0.5913320),  shape = 20, size = 3) +
  geom_point(data = exp_add2, aes(x = X1, y = X2), shape = 23, bg = 'red', size = 3) +
  geom_point(aes(x=0.05842766, y=-0.955879391),  shape = 22, bg = "orange", size = 3)
direct.label(v, method="bottom.pieces")


#################### step 3 
par(mar = c(2, 2, 1, 1))
contour(x.grid, y.grid, z.grid)
points(0, 0, pch = 17, col = "blue")
points(X, pch = 20)
points(-0.4656459, -0.5913320, pch = 20)
points(result$add[[1]], pch = 20)
points(0.05842766, -0.955879391, pch = 20)
points(result$add[[2]], pch = 20)

points(0.09004325, 0.001838837, pch = 22, bg = "green", cex = 0.8)
points(result$add[[3]], pch = 23, bg = "red", cex = 0.8)

add <- matrix(c(
  0.3870723,  0.9740065,
  1.8334596,  0.5650027,
  1.4243594, -0.7730981,
  0.8948653,  0.2365567), ncol = 2, byrow = T)
names(add) <- c("x1","x2")
exp_add3 <- data.frame(add, z = apply(add, 1, ackley))

v <- ggplot(tu, aes(x, y, z = z)) + stat_contour(aes(colour = ..level..)) + 
  geom_point(aes(x=0,y=0),  shape = 17, col = "blue", size = 3) +
  geom_point(data = exp_point, aes(x = x1, y = x2), shape = 20, size = 3) +
  geom_point(data = exp_add1, aes(x = X1, y = X2), shape = 20, size = 3) +
  geom_point(aes(x=-0.4656459, y=-0.5913320),  shape = 20, size = 3) +
  geom_point(data = exp_add2, aes(x = X1, y = X2), shape = 20, size = 3) +
  geom_point(aes(x=0.05842766, y=-0.955879391),  shape = 20, size = 3) +
  geom_point(data = exp_add3, aes(x = X1, y = X2), shape = 23, bg = 'red', size = 3) +
  geom_point(aes(x=0.09004325, y=0.001838837),  shape = 22, bg = "orange", size = 3)
direct.label(v, method="bottom.pieces")


#################### step 4
par(mar = c(2, 2, 1, 1))
contour(x.grid, y.grid, z.grid)
points(0, 0, pch = 17, col = "blue")
points(X, pch = 20)
points(-0.4656459, -0.5913320, pch = 20)
points(result$add[[1]], pch = 20)
points(0.05842766, -0.955879391, pch = 20)
points(result$add[[2]], pch = 20)
points(0.09004325,  0.001838837, pch = 20)
points(result$add[[3]], pch = 20)

points(0.04384017,  0.099930616, pch = 22, bg = "green", cex = 0.8)
points(result$add[[4]], pch = 23, bg = "red", cex = 0.8)

add <- matrix(c(
  0.1225721,  0.150049,
  1.8990791,  1.958523,
  -0.4770173, -1.141717,
  -0.9120836, -1.425537), ncol = 2, byrow = T)
names(add) <- c("x1","x2")
exp_add4 <- data.frame(add, z = apply(add, 1, ackley))

v <- ggplot(tu, aes(x, y, z = z)) + stat_contour(aes(colour = ..level..)) + 
  geom_point(aes(x=0,y=0),  shape = 17, col = "blue", size = 3) +
  geom_point(data = exp_point, aes(x = x1, y = x2), shape = 20, size = 3) +
  geom_point(data = exp_add1, aes(x = X1, y = X2), shape = 20, size = 3) +
  geom_point(aes(x=-0.4656459, y=-0.5913320),  shape = 20, size = 3) +
  geom_point(data = exp_add2, aes(x = X1, y = X2), shape = 20, size = 3) +
  geom_point(aes(x=0.05842766, y=-0.955879391),  shape = 20, size = 3) +
  geom_point(data = exp_add3, aes(x = X1, y = X2), shape = 20, size = 3) +
  geom_point(aes(x=0.09004325, y=0.001838837),  shape = 20, size = 3) +
  geom_point(data = exp_add4, aes(x = X1, y = X2), shape = 23, bg = 'red', size = 3) +
  geom_point(aes(x=0.04384017, y=0.099930616),  shape = 22, bg = "orange", size = 3)
direct.label(v, method="bottom.pieces")


#################### step 5
par(mar = c(2, 2, 1, 1))
contour(x.grid, y.grid, z.grid)
points(0, 0, pch = 17, col = "blue")
points(X, pch = 20)
points(-0.4656459, -0.5913320, pch = 20)
points(result$add[[1]], pch = 20)
points(0.05842766, -0.955879391, pch = 20)
points(result$add[[2]], pch = 20)
points(0.09004325,  0.001838837, pch = 20)
points(result$add[[3]], pch = 20)
points(0.04384017,  0.099930616, pch = 20)
points(result$add[[4]], pch = 20)

points(-0.02220093, -0.008929291, pch = 22, bg = "green", cex = 0.8)
points(result$add[[5]], pch = 23, bg = "red", cex = 0.8)

add <- matrix(c(
  -0.0979113, -0.02737333,
  0.6225918, -1.01618680,
  0.9640904, -0.40747161,
  -1.8089063, -1.77446276), ncol = 2, byrow = T)
names(add) <- c("x1","x2")
exp_add5 <- data.frame(add, z = apply(add, 1, ackley))

v <- ggplot(tu, aes(x, y, z = z)) + stat_contour(aes(colour = ..level..)) + 
  geom_point(aes(x=0,y=0),  shape = 17, col = "blue", size = 3) +
  geom_point(data = exp_point, aes(x = x1, y = x2), shape = 20, size = 3) +
  geom_point(data = exp_add1, aes(x = X1, y = X2), shape = 20, size = 3) +
  geom_point(aes(x=-0.4656459, y=-0.5913320),  shape = 20, size = 3) +
  geom_point(data = exp_add2, aes(x = X1, y = X2), shape = 20, size = 3) +
  geom_point(aes(x=0.05842766, y=-0.955879391),  shape = 20, size = 3) +
  geom_point(data = exp_add3, aes(x = X1, y = X2), shape = 20, size = 3) +
  geom_point(aes(x=0.09004325, y=0.001838837),  shape = 20, size = 3) +
  geom_point(data = exp_add4, aes(x = X1, y = X2), shape = 20, size = 3) +
  geom_point(aes(x=0.04384017, y=0.099930616),  shape = 20, size = 3) +
  
  geom_point(data = exp_add5, aes(x = X1, y = X2), shape = 23, bg = 'red', size = 3) +
  geom_point(aes(x=-0.02220093, y=-0.008929291),  shape = 22, bg = "orange", size = 3)
direct.label(v, method="bottom.pieces")
