##### simulate data based on two clusters ####
mu1 = function(x)
{
  2*(x/50)^(0.2)
}

mu1 = Vectorize(mu1)

mu2 = function(x)
{
  -2*(x/50)^(0.2)
}

mu2 = Vectorize(mu2)


# eigen functions
psi1 = function(x)
{
  sqrt(2)*sin(2*pi*x)
}

psi1 = Vectorize(psi1)


psi2 = function(x)
{
  sqrt(2)*cos(2*pi*x)
}

psi2 = Vectorize(psi2)




### ncl number of observations for each group
### mvec is the range of number of observations for each subject
## time is simulated from uniform(0,1)

# sig2 = 0.1
# lamj = c(0.1,0.2)
# mvec = c(5, 20)
# ncl = 50
simdat = function(sig2, lamj, mvec, ncl = 50, seed = 2228 )
{
  set.seed(seed)
  
  ## number of observations for each subject
  
  nsub = sample(mvec[1]:mvec[2], ncl*2, replace = TRUE)
  ntotal = sum(nsub)
  group = rep(1:2, each = ncl)
  
  dat = data.frame(group = rep(group, nsub),
                    ind = rep(1:(2*ncl),nsub),
                    time = runif(ntotal),
                    obs = rep(0, ntotal))
  
  rande = sqrt(sig2) * rnorm(ntotal)
  
  mean1 = mu1(dat$time[dat$group==1])
  mean2 = mu2(dat$time[dat$group==2])
  
  vi = rep(0, ntotal)
  for(i in 1:(ncl*2))
  {
    timei = dat$time[dat$ind==i]
    vi[dat$ind ==i] = sqrt(lamj[1]) * rnorm(1) * psi1(timei) +
      sqrt(lamj[2]) * rnorm(1) * psi2(timei)
  }
  
  
  dat$obs = c(mean1, mean2) + vi + rande
  
  return(dat)

}



##### simulate data based on spline and coefficients not curves####
library(orthogonalsplinebasis)
simdat1 = function(knots, betag, sig2, lamj, mvec, ncl = 50, seed = 2228, boundary = c(0,1) )
{
  set.seed(seed)
  
  ## number of observations for each subject
  
  nsub = sample(mvec[1]:mvec[2], ncl*2, replace = TRUE)
  ntotal = sum(nsub)
  group = rep(1:2, each = ncl)
  
  dat = data.frame(group = rep(group, nsub),
                   ind = rep(1:(2*ncl),nsub),
                   time = runif(ntotal),
                   obs = rep(0, ntotal))

  knotsall = c(rep(boundary[1],4),knots, rep(boundary[2],4))
  obasisobj = OBasis(knotsall)
  Bm = evaluate(obasisobj,dat$time)  ## orthogonal
  
  rande = sqrt(sig2) * rnorm(ntotal)
  
  meanvec = rep(0, ntotal)
  meanvec[dat$group==1] = Bm[dat$group==1,] %*% betag[1,]
  meanvec[dat$group==2] = Bm[dat$group==2,] %*% betag[2,]
  
  vi = rep(0, ntotal)
  for(i in 1:(ncl*2))
  {
    timei = dat$time[dat$ind==i]
    vi[dat$ind ==i] = sqrt(lamj[1]) * rnorm(1) * psi1(timei) +
      sqrt(lamj[2]) * rnorm(1) * psi2(timei)
  }
  
  
  dat$obs = meanvec + vi + rande
  
  return(dat)
  
}









