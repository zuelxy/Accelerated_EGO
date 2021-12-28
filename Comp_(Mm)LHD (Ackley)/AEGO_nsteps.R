AEGO_nsteps <- function(QMC, fun, model, npoints, nsteps, lower, upper, 
                  minimization = TRUE, control = NULL)
{
  history <- c()
  n <- nrow(model@X)
  for (i in 1:nsteps) {
    print(paste("step :", i, "/", nsteps))
    
    x_new <- max_EI(model = model, lower = lower, upper = upper, control = control, minimization = minimization)
    
    if(i==1)
      x_sampling <- resampling(k = npoints-1, RQMC = QMC, model = model, minimization = minimization)
    else{
      RQMC <- randomization(QMC, lower = lower, upper = upper)
      x_sampling <- resampling(k = npoints-1, RQMC = RQMC, model = model, minimization = minimization)
    }
    
    add <- rbind(x_new$par, x_sampling)
    
    model <- update(object = model, newX = add, 
                    newy = as.matrix(apply(add,1,fun)), 
                    newX.alreadyExist = FALSE, cov.reestim = TRUE, 
                    kmcontrol = list(control = list(trace = FALSE)))
    history[i] <- min(model@y)*minimization + max(model@y)*(!minimization)
    
  }
  return(list(par = model@X[(n+1):(n+nsteps*npoints),,drop = FALSE],
              value = model@y[(n+1):(n+nsteps*npoints),,drop = FALSE], history = history))
}