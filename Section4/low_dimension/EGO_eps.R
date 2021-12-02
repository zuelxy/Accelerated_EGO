EGO_eps <- function (model, fun, optimum, eps, lower, upper, parinit = NULL, control = NULL, 
          kmcontrol = NULL) 
{
  history <- c()
  nsteps <- 1
  while(TRUE){
    oEGO <- max_EI(model = model, lower = lower, upper = upper, 
                   parinit = parinit, control = control)
    history[nsteps] <- fun(t(oEGO$par))
    if(abs(history[nsteps]-optimum)<eps)
      break
    model <- update(object = model, newX = oEGO$par, 
                    newy = as.matrix(fun(t(oEGO$par))), 
                    newX.alreadyExist = FALSE, cov.reestim = TRUE, 
                    kmcontrol = list(control = list(trace = FALSE)))
    nsteps <- nsteps + 1
  }
  return(list(nsteps = nsteps, history = history))
}