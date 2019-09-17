
###### two groups setting #####
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



#### ISE of two functions (without mean now) ####
### lower is the lower bound and upper is the upper bound 
## obj is estimate obj of FDAsubgroup 
ISEFDA = function(obj, lower = 0, upper = 1)
{
  obasisobj = obj$obasisobj
  betaest = obj$betaest
  groupest = obj$groupest
  fundiff = function(x, ind)
  {
    bmx = evaluate(res$obasisobj,x)
    meanest = bmx %*% betaest[ind,]
    groupi = groupest[ind]
    if(groupi == 1)
    {
      mean0 = mu1(x)
    }
    if(groupi == 2)
    {
      mean0 = mu2(x)
    }
    
    meandiff = (meanest - mean0)^2
    return(meandiff)
  }
  

  ISEvec = unlist(lapply(1:nrow(betaest),function(ind){stats::integrate(function(x){fundiff(x,ind)},lower = 0,upper = 1)$value}))
 
  return(ISEvec)
}



#### betaest is the estimate of Zhu and qu, ncol of columns is 100
ISEind = function(betaest, timerange, order, nknots)
{
  fundiff = function(x,ind)
  {
    bmx = bsplineS(x, knots_eq3(timerange, k = order, m = nknots), norder = order)
    meanest = bmx %*% (betaest[,ind])
    groupi = groupest[ind]
    if(groupi == 1)
    {
      mean0 = mu1(x)
    }
    if(groupi == 2)
    {
      mean0 = mu2(x)
    }
    
    meandiff = (meanest - mean0)^2
    return(meandiff)
  }
  
  ISEvec = unlist(lapply(1:ncol(betaest),function(ind){stats::integrate(function(x){fundiff(x,ind)},lower = 0,upper = 1)$value}))
  return(ISEvec)
  
}


#### ISE for James and 
