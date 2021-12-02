library(DiceKriging)
library(DiceOptim)
library(rgenoud)
library(randtoolbox)
library(readxl)     # load excel file
library(e1071)      #  package for SVM
library(parallel)
library(doParallel)

source("K:/Accelerated_EGO_code/section_5/algorithms_code/AEGO_hpo.R")
source("K:/Accelerated_EGO_code/section_5/algorithms_code/EGO_hpo.R")
source("K:/Accelerated_EGO_code/section_5/algorithms_code/qEGO_hpo.R")
source("K:/Accelerated_EGO_code/affine.R")

######### data prepare #######
glass <- read_excel("K:/Accelerated_EGO_code/section_5/svm_data/glass.xlsx")
glass$class <- factor(glass$class)    # transform 'class' to factor type
str(glass) # check the type of data

set.seed(305)
index <- sample(1:nrow(glass), size = 0.75*nrow(glass))
train <- as.data.frame(glass[index,])    
test <- as.data.frame(glass[-index,])    # divide the train and test data 

x <- subset(train, select = -class)
y <- train$class

################################### grid ##################################
rep_train <- c()
rep_test <- c()
rep_time <- c()
for (k in 1:10) {
  begin_g <- proc.time()
  obj <- tune.svm(x, y, gamma = 2^(-15:5), cost = 2^(0:20), 
                  tune.control(sampling = "cross", cross = 5))
  end_g <- proc.time()-begin_g

  p <- obj$best.parameters
###### predict 
  svm_model_g <- svm(x, y, kernel = "radia", gamma = p[1] , cost = p[2])

  pred_g <- predict(svm_model_g, test[,-1])
  ta_g <- table(pred_g, test[,1])
  result_g <- sum(diag(ta_g))/sum(ta_g)
  
  rep_train[k] <- 1-obj$best.performance
  rep_test[k] <- result_g
  rep_time[k] <- end_g[1]
}


###############################  response function #############################
svm_accuracy <- function(xx, svmx, svmy){
  x1 <- xx[1]
  x2 <- xx[2]
  svm_model <- svm(svmx, svmy, kernel = "radia", gamma = 2^(x1), cost = 2^(x2), cross = 5)
  return(svm_model$tot.accuracy)
}

# a <- c(-20, 0)      # the range of gamma 
# b <- c(-5, 15)        # the range of cost
lower <- c(-15, 0)
upper <- c(5, 20)

U <- read.table("K:/Accelerated_EGO_code/section_5/UD21.txt", header = T) 
X <- experiments(U, lower, upper) 
accuracy_X <- as.matrix(apply(X, 1, svm_accuracy, x, y))

fitted_model <- km(design = X, response = accuracy_X, optim.method = "gen",
                   nugget=1e-4*var(accuracy_X), 
                   control=list(trace = FALSE))

q <- 5
train_accuracy <- matrix(nrow = 100, ncol = 3)
test_accuracy <- matrix(nrow = 100, ncol = 3)
CPU_time <- matrix(nrow = 100, ncol = 3)

for(k in 1:100){
  
################################### A_EGO ################################
  QMC <- sobol(100, dim = 2, scrambling = 3)
  QMC <- affine(QMC, lower, upper)
  
  begin_AEGO <- proc.time()
  result <- AEGO_nsteps(QMC = QMC, fun = svm_accuracy, svmx = x, svmy = y,
                        model = fitted_model, npoints = q, 
                        nsteps = 4, lower = lower, upper = upper, minimization = FALSE,
                        control = list(pop.size = 50, max.generations = 50, 
                                       wait.generations = 5, BFGSburnin = 5))
  end_AEGO <- proc.time()-begin_AEGO
  
  best_point <- matrix(result$par[which(result$value==max(result$value)),], ncol = 2)
  
  svm_model_A <- svm(x, y, kernel = "radia", 
                     gamma = 2^(best_point[1,1]), cost = 2^(best_point[1,2]))
  pred_A <- predict(svm_model_A, test[,-1])
  ta_A <- table(pred_A, test[,1])
  accuracy_A <- sum(diag(ta_A))/sum(ta_A)
  
  ################################# constant liar ################################
  begin_CL <- proc.time()
  
  result_CL <- qEGO_nsteps(fun = svm_accuracy, svmx = x, svmy = y, 
                           model = fitted_model, npoints = q,
                           nsteps = 4, lower = lower, upper = upper, crit = "CL", 
                           minimization = FALSE, 
                           optimcontrol = list(L = "min",pop.size = 50, max.generations = 50, 
                                               wait.generations = 5, BFGSburnin = 5))
  end_CL <- proc.time()-begin_CL
  
  best_point1 <- matrix(result_CL$par[which(result_CL$value==max(result_CL$value)),],ncol = 2)
  
  svm_model_C <- svm(x, y, kernel = "radia", 
                     gamma = 2^(best_point1[1,1]), cost = 2^(best_point1[1,2]))
  pred_C <- predict(svm_model_C, test[,-1])
  ta_C <- table(pred_C, test[,1])
  accuracy_C <- sum(diag(ta_C))/sum(ta_C)
  
  ############################# EGO #############################
  nsteps <- 21
  begin_EGO <- proc.time()
  result_EGO <- EGO_nsteps(model = fitted_model, fun = svm_accuracy, svmx = x, svmy =y,
                           lower = lower, upper = upper, nsteps = nsteps, minimization = FALSE,
                           control = list(pop.size = 50, max.generations = 50, 
                                          wait.generations = 5, BFGSburnin = 5))
  end_EGO <- proc.time()-begin_EGO
  
  
  best_point_E <- matrix(result_EGO$par[which(result_EGO$value==max(result_EGO$value)),],ncol = 2) 
  svm_model_E <- svm(x, y, kernel = "radia", 
                     gamma = 2^(best_point_E[1,1]), cost = 2^(best_point_E[1,2]))
  pred_E <- predict(svm_model_E, test[,-1])
  ta_E <- table(pred_E, test[,1])
  accuracy_E <- sum(diag(ta_E))/sum(ta_E)
  
  
  train_accuracy[k, ] <- c(max(result$value), max(result_CL$value), max(result_EGO$value))
  test_accuracy[k, ] <- c(accuracy_A, accuracy_C, accuracy_E)
  CPU_time[k, ] <- c(end_AEGO[1], end_CL[1], end_EGO[1])
}
