
library(LHD)
source("section_3/random_LHD.R")


pool <- c()
for (k in 1:100) {
result_LHD <- rLHD_EGO(fun = ackley, model = fitted_model, npoints = 5, 
                optimum = M, eps = 10^(-1), lower = lower, upper = upper, m=100, 
                control = list(pop.size = 50, max.generations = 60, 
                               wait.generations = 5, BFGSburnin = 5))
pool <- c(pool, result_LHD$nsteps)
}

################## store result
# result_LHD$value
#[1] 2.62447560 1.65404859 0.85092124 0.19864693 0.19864693 0.19864693 0.11374153 0.11374153 0.11374153
#[10] 0.11374153 0.04944763

# add_5 <- rbind(result_LHD$add[[1]],result_LHD$add[[2]],result_LHD$add[[3]],
#             result_LHD$add[[4]],result_LHD$add[[5]], result_LHD$par[1:5,])
# write.csv(add_5, "section_3/addpoint_LHD.csv")


x.grid <- seq(-2,2, length.out =20)
y.grid <- seq(-2,2, length.out =20)
design.grid <- expand.grid(x.grid, y.grid)
colnames(design.grid) <- c("x1","x2")
response.grid <- apply(design.grid,  1, ackley)

tu <- data.frame(design.grid, z = response.grid)
exp_point <- data.frame(X, z = response_X)

##########3 step1
add <- result_LHD$add[[1]]
names(add) <- c("x1","x2")
exp_add1 <- data.frame(add, z = apply(add, 1, ackley))

v <- ggplot(tu, aes(x=x1, y=x2, z = z)) + stat_contour(aes(colour = ..level..)) + 
  geom_point(aes(x=0,y=0),  shape = 17, col = "blue", size = 3) +
  geom_point(data = exp_point, aes(x = X1, y = X2), shape = 20, size = 2) +
  geom_point(data = exp_add1, aes(x = X1, y = X2), shape = 23, bg = 'red', size = 3) +
  geom_point(aes(x=result_LHD$par[1, 1], y=result_LHD$par[1, 2]),  shape = 22, bg = "orange", size = 3)
direct.label(v, method="bottom.pieces")

########### step2
add <- result_LHD$add[[2]]
names(add) <- c("x1","x2")
exp_add2 <- data.frame(add, z = apply(add, 1, ackley))

v <- ggplot(tu, aes(x=x1, y=x2, z = z)) + stat_contour(aes(colour = ..level..)) + 
  geom_point(aes(x=0,y=0),  shape = 17, col = "blue", size = 3) +
  geom_point(data = exp_point, aes(x = X1, y = X2), shape = 20, size = 2) +
  geom_point(data = exp_add1, aes(x = X1, y = X2), shape = 20, size = 2) +
  geom_point(aes(x=result_LHD$par[1, 1], y=result_LHD$par[1, 2]),  shape = 20, size = 2) +
  geom_point(data = exp_add2, aes(x = X1, y = X2), shape = 23, bg = 'red', size = 3) +
  geom_point(aes(x=result_LHD$par[2, 1], y=result_LHD$par[2, 2]),  shape = 22, bg = "orange", size = 3)
direct.label(v, method="bottom.pieces")

########### step3

add <- result_LHD$add[[3]]
names(add) <- c("x1","x2")
exp_add3 <- data.frame(add, z = apply(add, 1, ackley))

v <- ggplot(tu, aes(x=x1, y=x2, z = z)) + stat_contour(aes(colour = ..level..)) + 
  geom_point(aes(x=0,y=0),  shape = 17, col = "blue", size = 3) +
  geom_point(data = exp_point, aes(x = X1, y = X2), shape = 20, size = 2) +
  geom_point(data = exp_add1, aes(x = X1, y = X2), shape = 20, size = 2) +
  geom_point(aes(x=result_LHD$par[1, 1], y=result_LHD$par[1, 2]),  shape = 20, size = 2) +
  geom_point(data = exp_add2, aes(x = X1, y = X2), shape = 20, size = 2) +
  geom_point(aes(x=result_LHD$par[2, 1], y=result_LHD$par[2, 2]),  shape = 20, size = 2) +
  geom_point(data = exp_add3, aes(x = X1, y = X2), shape = 23, bg = 'red', size = 3) +
  geom_point(aes(x=result_LHD$par[3, 1], y=result_LHD$par[3, 2]),  shape = 22, bg = "orange", size = 3)
direct.label(v, method="bottom.pieces")


############### step 4
add <- result_LHD$add[[4]]
names(add) <- c("x1","x2")
exp_add4 <- data.frame(add, z = apply(add, 1, ackley))

v <- ggplot(tu, aes(x=x1, y=x2, z = z)) + stat_contour(aes(colour = ..level..)) + 
  geom_point(aes(x=0,y=0),  shape = 17, col = "blue", size = 3) +
  geom_point(data = exp_point, aes(x = X1, y = X2), shape = 20, size = 2) +
  geom_point(data = exp_add1, aes(x = X1, y = X2), shape = 20, size = 2) +
  geom_point(aes(x=result_LHD$par[1, 1], y=result_LHD$par[1, 2]),  shape = 20, size = 2) +
  geom_point(data = exp_add2, aes(x = X1, y = X2), shape = 20, size = 2) +
  geom_point(aes(x=result_LHD$par[2, 1], y=result_LHD$par[2, 2]),  shape = 20, size = 2) +
  geom_point(data = exp_add3, aes(x = X1, y = X2), shape = 20, size = 2) +
  geom_point(aes(x=result_LHD$par[3, 1], y=result_LHD$par[3, 2]),  shape = 20, size = 2) +
  geom_point(data = exp_add4, aes(x = X1, y = X2), shape = 23, bg = 'red', size = 3) +
  geom_point(aes(x=result_LHD$par[4, 1], y=result_LHD$par[4, 2]),  shape = 22, bg = "orange", size = 3)
direct.label(v, method="bottom.pieces")


############## step 5
add <- result_LHD$add[[5]]
names(add) <- c("x1","x2")
exp_add5 <- data.frame(add, z = apply(add, 1, ackley))

v <- ggplot(tu, aes(x=x1, y=x2, z = z)) + stat_contour(aes(colour = ..level..)) + 
  geom_point(aes(x=0,y=0),  shape = 17, col = "blue", size = 3) +
  geom_point(data = exp_point, aes(x = X1, y = X2), shape = 20, size = 2) +
  geom_point(data = exp_add1, aes(x = X1, y = X2), shape = 20, size = 2) +
  geom_point(aes(x=result_LHD$par[1, 1], y=result_LHD$par[1, 2]),  shape = 20, size = 2) +
  geom_point(data = exp_add2, aes(x = X1, y = X2), shape = 20, size = 2) +
  geom_point(aes(x=result_LHD$par[2, 1], y=result_LHD$par[2, 2]),  shape = 20, size = 2) +
  geom_point(data = exp_add3, aes(x = X1, y = X2), shape = 20, size = 2) +
  geom_point(aes(x=result_LHD$par[3, 1], y=result_LHD$par[3, 2]),  shape = 20, size = 2) +
  geom_point(data = exp_add4, aes(x = X1, y = X2), shape = 20, size = 2) +
  geom_point(aes(x=result_LHD$par[4, 1], y=result_LHD$par[4, 2]),  shape = 20, size = 2) +
  geom_point(data = exp_add5, aes(x = X1, y = X2), shape = 23, bg = 'red', size = 3) +
  geom_point(aes(x=result_LHD$par[5, 1], y=result_LHD$par[5, 2]),  shape = 22, bg = "orange", size = 3)
direct.label(v, method="bottom.pieces")


