####### simulate data with common regression coefficients #####
##### x dimension is nx and are iid normal, meanx mean of x, sdx sd of x
#### etavec is the common coefficents 

#etavec = runif(n = 2,1,2)


#### time is randomly sampled from uniform distribution Unif(0,1)


simdatv2 = function(xlist = list(nx = 2, meanx = 0, sdx = 1, eta),
                   sig2, lamj, mvec = c(5,20), ncl = 50, 
                   funlist, eigenlist, seed = 2228)
{
  set.seed(seed)
  ngroup = length(funlist)
  ## number of observations for each subject

  nsub = sample(mvec[1]:mvec[2], ncl*ngroup, replace = TRUE)
  ntotal = sum(nsub)
  group = rep(1:ngroup, each = ncl)
  
  xmat = matrix(rnorm(ntotal*xlist$nx,xlist$meanx, xlist$sdx), nrow = ntotal)
  
  # dat = data.frame(group = rep(group, nsub),
  #                  x = xmat,
  #                  ind = rep(1:(ngroup*ncl),nsub),
  #                  time = sample(grids,ntotal,replace = TRUE),
  #                  obs = rep(0, ntotal),
  #                  mean = rep(0, ntotal))
  
  dat = data.frame(group = rep(group, nsub),
                   x = xmat,
                   ind = rep(1:(ngroup*ncl),nsub),
                   time = runif(ntotal),
                   obs = rep(0, ntotal),
                   mean = rep(0, ntotal))
  
  colnames(dat)[2:(xlist$nx+1)] = paste0("x",1:xlist$nx)
  
  rande = sqrt(sig2) * rnorm(ntotal)
  
  for(j in 1:ngroup)
  {
    dat$mean[dat$group==j] = funlist[[j]](dat$time[dat$group==j])
  }
  
  vi = rep(0, ntotal)
  for(i in 1:(ncl*ngroup))
  {
    timei = dat$time[dat$ind==i]
    vi[dat$ind ==i] = sqrt(lamj[1]) * rnorm(1) * eigenlist[[1]](timei) +
      sqrt(lamj[2]) * rnorm(1) * eigenlist[[2]](timei)
  }
  
  dat$obs = xmat %*% xlist$eta + dat$mean + vi + rande
  return(dat)
}