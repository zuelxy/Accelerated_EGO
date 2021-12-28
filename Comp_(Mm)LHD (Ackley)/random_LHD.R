rLHD_EGO <- function(fun, model, npoints, optimum, eps, lower, upper, m,
                  minimization = TRUE, control = NULL)
{
  value <- c()             # store the  global optimum each step
  EI_point <- c()          # the point with the lagest EI value each step
  sampling_point <- list() # the points obatained through sampling process each step
  
  i <- 1
  while(TRUE){ 
    x_new <- max_EI(model = model, lower = lower, upper = upper, control = control, 
                    minimization = minimization) # maximing the EI function
    
    LHD <- rLHD(m, length(lower))
    LHD_domain <- experiments(LHD, lower, upper)
    x_sampling <- resampling(k = npoints-1, RQMC = LHD_domain, model = model, minimization = minimization) 
    add <- rbind(x_new$par, x_sampling)
    model <- update(object = model, newX = add, 
                    newy = as.matrix(apply(add, 1, fun)), 
                    newX.alreadyExist = FALSE, cov.reestim = TRUE, 
                    kmcontrol = list(control = list(trace = FALSE)))  # refit the model 
    
    value[i] <- min(model@y)*minimization + max(model@y)*(!minimization)
    EI_point <- rbind(EI_point, x_new$par)
    sampling_point[[i]] <- x_sampling
    
    if(abs(fun(x_new$par)-optimum)<eps)
      break                               # stop criterion
    i <- i+1
  }
  return(list(nsteps = i, par = EI_point, value = value, add = sampling_point))
}