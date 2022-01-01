######## simulation for x grouping ######

### time is simulated from a uniform distribution Unif(0,1) 

### yih = xietai + mui(tih) + eih
## nx is the number of x 
## etag group parameter for x ## nx \times ng matrix, each column is for one group 

##### parametric part is the same for the same individual ####
simdatx = function(xlist = list(nx = 2, meanx = 0, sdx = 1, etag),
                              sig2, lamj, mvec = c(10,20), ncl = 50, 
                              funlist, eigenlist, seed = 2228)
{
  set.seed(seed)
  ngroup = length(funlist)
  
  nobs = ncl*ngroup
  if(mvec[1]== mvec[2])
  {
    nsub = rep(mvec[1], nobs)
  }else{
    nsub = sample(mvec[1]:mvec[2], nobs, replace = TRUE)
  }
  ntotal = sum(nsub)
  group = rep(1:ngroup, each = ncl)
  
  xmat = matrix(rnorm(nobs*xlist$nx,xlist$meanx, xlist$sdx), nrow = nobs)
  
  dat = data.frame(group = rep(group, nsub),
                   ind = rep(1:(ngroup*ncl),nsub),
                   time = runif(ntotal),
                   x = xmat[rep(1:nobs, nsub),],
                   obs = rep(0, ntotal),
                   mean = rep(0, ntotal),
                   meanx = rep(0, ntotal))
  
  
  ### mean is the mean curve function #####
  
  colnames(dat)[4:(xlist$nx+3)] = paste0("x",1:xlist$nx)

  rande = sqrt(sig2) * rnorm(ntotal)
  
  for(j in 1:ngroup)
  {
    dat$mean[dat$group==j] =  funlist[[j]](dat$time[dat$group==j])
    dat$meanx[dat$group==j] = as.matrix(dat[dat$group==j,4:(xlist$nx+3)]) %*% xlist$etag[,j]
  }
  
  vi = rep(0, ntotal)
  for(i in 1:(ncl*ngroup))
  {
    timei = dat$time[dat$ind==i]
    vi[dat$ind ==i] = sqrt(lamj[1]) * rnorm(1) * eigenlist[[1]](timei) +
      sqrt(lamj[2]) * rnorm(1) * eigenlist[[2]](timei)
  }
  dat$obs = dat$meanx + dat$mean + vi + rande
  
  return(dat)
}


####### parametric part is not the same for the same individual ####

### yih = xih etai + mui(tih) + eih


simdatx2 = function(xlist = list(nx = 2, meanx = 0, sdx = 1, etag),
                   sig2, lamj, mvec = c(10,20), ncl = 50, 
                   funlist, eigenlist, seed = 2228)
{
  set.seed(seed)
  ngroup = length(funlist)
  
  nobs = ncl*ngroup
  if(mvec[1]== mvec[2])
  {
    nsub = rep(mvec[1], nobs)
  }else{
    nsub = sample(mvec[1]:mvec[2], nobs, replace = TRUE)
  }
  ntotal = sum(nsub)
  group = rep(1:ngroup, each = ncl)
  
  xmat = matrix(rnorm(ntotal*xlist$nx,xlist$meanx, xlist$sdx), nrow = ntotal)
  
  dat = data.frame(group = rep(group, nsub),
                   ind = rep(1:(ngroup*ncl),nsub),
                   time = runif(ntotal),
                   x = xmat,
                   obs = rep(0, ntotal),
                   mean = rep(0, ntotal),
                   meanx = rep(0, ntotal))
  
  
  ### mean is the mean curve function #####
  
  colnames(dat)[4:(xlist$nx+3)] = paste0("x",1:xlist$nx)
  
  rande = sqrt(sig2) * rnorm(ntotal)
  
  for(j in 1:ngroup)
  {
    dat$mean[dat$group==j] =  funlist[[j]](dat$time[dat$group==j])
    dat$meanx[dat$group==j] = as.matrix(dat[dat$group==j,4:(xlist$nx+3)]) %*% xlist$etag[,j]
  }
  
  vi = rep(0, ntotal)
  for(i in 1:(ncl*ngroup))
  {
    timei = dat$time[dat$ind==i]
    vi[dat$ind ==i] = sqrt(lamj[1]) * rnorm(1) * eigenlist[[1]](timei) +
      sqrt(lamj[2]) * rnorm(1) * eigenlist[[2]](timei)
  }
  dat$obs = dat$meanx + dat$mean + vi + rande
  
  return(dat)
}