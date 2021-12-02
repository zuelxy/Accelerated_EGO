qEGO_eps <- function (fun, model, npoints, optimum, eps, lower = rep(0, model@d), 
          upper = rep(1, model@d), crit = "CL", minimization = TRUE, 
          optimcontrol = NULL, cov.reestim = TRUE, ...) 
{
  history <- c()
  d <- model@d
  parinit <- optimcontrol$parinit
  i <- 1
  while(TRUE){
    res <- max_qEI(model = model, npoints = npoints, lower = lower, 
                   upper = upper, crit = crit, minimization = minimization, 
                   optimcontrol = optimcontrol)
    model <- update(object = model, newX = res$par, newy = as.matrix(apply(res$par, 
                                                                           1, fun, ...)), newX.alreadyExist = FALSE, 
                    cov.reestim = cov.reestim, kmcontrol = list(control = list(trace = FALSE)))
    history[i] <- min(model@y) * minimization + max(model@y) * 
      (!minimization)
    if(abs(history[i]-optimum)<eps)
      break
    i <- i+1  
  }
  return(list(nsteps = i, history = history))
}