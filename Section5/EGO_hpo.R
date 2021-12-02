EGO_nsteps <- function (model, fun, svmx, svmy, nsteps, lower, upper, 
                        minimization = TRUE, parinit = NULL, control = NULL, kmcontrol = NULL) 
{
  value <- c()
  par <- c()
  for (i in 1:nsteps) {
    oEGO <- max_EI(model = model, lower = lower, upper = upper, minimization = minimization,
                   parinit = parinit, control = control)
    # model@X <- rbind(model@X, oEGO$par)
    # model@y <- rbind(model@y, fun(t(oEGO$par)))
    value[i] <- fun(oEGO$par, svmx, svmy)
    par <- rbind(par, oEGO$par)
    model <- update(object = model, newX = oEGO$par, 
                    newy = as.matrix(value[i]), 
                    newX.alreadyExist = FALSE, cov.reestim = TRUE, 
                    kmcontrol = list(control = list(trace = FALSE)))
  }
  return(list(nsteps = i, par = par, value = value))
}
