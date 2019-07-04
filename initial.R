####### initial value based on GCV for beta when assuming no correlation ####


library(orthogonalsplinebasis)

# ind is ID
# tm is the time column
# y is the observation
initialcoef = function(ind, tm, y, knots, boundary = c(0,1), lamv)
{
  ntotal = length(y)
  knotsall = c(rep(boundary[1],4),knots, rep(boundary[2],4))
  obasisobj = OBasis(knotsall)
  Bm = evaluate(obasisobj,tm)  ## orthogonal
  nlam = length(lamv)
  
  ### ni for each subject ###
  uind = unique(ind)
  n = length(uind)
  nvec = rep(0,n)
  for(i in 1:n)
  {
    nvec[i] = sum(ind==uind[i])
  }
  
  p = ncol(Bm)
  Imp = diag(1, p, p)
  
  gcvmat = matrix(0, n, nlam)
  for(i in 1:n)
  {
    indi = ind == uind[i]
    yi = y[indi]
    Bmi = Bm[indi,]
    for(j in 1:nlam)
    {
      gcvmat[i,j] = gcvi(yi, Bmi, Imp, nvec[i], lamv[j])
    }
  }
  
  gcvvec = colSums(gcvmat)
  
  lamselect = lamv[which.min(gcvvec)]
  
  betam = matrix(0, n, p)
  for(i in 1:n)
  {
    indi = ind == uind[i]
    yi = y[indi]
    Bmi = Bm[indi,]
    betam[i,] = solve(t(Bmi) %*% Bmi + lamselect * Imp) %*% t(Bmi) %*%yi
  }
  
  
  return(betam)
  
}

# Imp is the identity matrix
gcvi = function(yi, Bmi, Imp, ni, lam1)
{
  Dmlam = lam1 * Imp
  Hmat = Bmi %*% solve(t(Bmi)%*% Bmi + Dmlam) %*% t(Bmi)
  Ini = diag(1, ni, ni)
  Ih = Ini - Hmat
  trh = sum(diag(Ih))
  yhi = Ih %*% yi
  gcvvalue = t(yhi) %*% yhi* ni/(trh^2)
  return(gcvvalue)
}


