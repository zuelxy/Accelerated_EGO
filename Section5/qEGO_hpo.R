qEGO_nsteps <- function (fun,svmx,svmy, model, npoints, nsteps, lower = rep(0, model@d), 
          upper = rep(1, model@d), crit = "exact", minimization = TRUE, 
          optimcontrol = NULL, cov.reestim = TRUE, ...) 
{
  history <- rep(NaN, nsteps)
  d <- model@d
  n <- nrow(model@X)
  parinit <- optimcontrol$parinit
  for (i in 1:nsteps) {
    print(paste("step :", i, "/", nsteps))

    res <- max_qEI(model = model, npoints = npoints, lower = lower, 
                   upper = upper, crit = crit, minimization = minimization, 
                   optimcontrol = optimcontrol)
    
    cl <- makeCluster(5)
    registerDoParallel(cl)
    o <- foreach(l=1:5, .combine = 'rbind', .packages = 'e1071') %dopar% {
      fun(res$par[l,], svmx, svmy)
    }
    stopCluster(cl)
    
    model <- update(object = model, newX = res$par, 
                    newy = o, newX.alreadyExist = FALSE, 
                    cov.reestim = cov.reestim, kmcontrol = list(control = list(trace = FALSE)))
    
    history[i] <- min(model@y) * minimization + max(model@y) * (!minimization)
    }
  return(list(par = model@X[(n + 1):(n + nsteps * npoints), , drop = FALSE], 
              value = model@y[(n + 1):(n + nsteps * npoints), , drop = FALSE], 
              npoints = npoints, nsteps = nsteps, lastmodel = model, history = history))
}