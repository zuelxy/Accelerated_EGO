AEGO_eps <- function(QMC, fun, model, npoints, optimum, eps, lower, upper, 
                  minimization = TRUE, control = NULL)
{
  value <- c()
  i <- 1
  
  while(TRUE){
    x_new <- max_EI(model = model, lower = lower, upper = upper, control = control, 
                    minimization = minimization)

    if(i==1)
      x_sampling <- resampling(k = npoints-1, RQMC = QMC, model = model, minimization = minimization)
    else{
      RQMC <- randomization(QMC, lower = lower, upper = upper)
      x_sampling <- resampling(k = npoints-1, RQMC = RQMC, model = model, minimization = minimization)
    }
    add <- rbind(x_new$par, x_sampling)
    model <- update(object = model, newX = add, 
                   newy = as.matrix(apply(add, 1, fun)), 
                   newX.alreadyExist = FALSE, cov.reestim = TRUE, 
                   kmcontrol = list(control = list(trace = FALSE)))
    
    value[i] <- min(model@y)*minimization+max(model@y)*(!minimization)
    
    if(abs(value[i]-optimum)<eps)
      break
    i <- i+1
  }
  return(list(nsteps = i, value = value))
}
