####### initial value based on GCV for beta when assuming no correlation ####


library(orthogonalsplinebasis)

# ind is ID
# tm is the time column
# y is the observation
# based on GCV
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


#### each ind has its own coefficients ####

initial_indiv = function(ind, y, xm, lam = 0.01)
{
  ntotal = length(y)
  
  uind = unique(ind)
  n = length(uind)
  
  p = ncol(xm)
  Imp = diag(1, p, p)

  betam = matrix(0, n, p)
  
  for(i in 1:n)
  {
    indi = ind == uind[i]
    yi = y[indi]
    xmi = xm[indi,]
    betam[i,] = solve(t(xmi) %*% xmi + lam * Imp) %*% t(xmi) %*%yi
  }
  
  return(betam)
}


###### initial value based on the difference (Laplacian) matrix #####
# 
# ind = dat$ind
# tm = dat$time
# y = dat$obs
# knots = seq(0,1,length.out = 5)[2:4]

## initialcoef2 --> initiallap
initiallap = function(ind, tm, y, knots, boundary = c(0,1), lam = 0.001)
{
  ntotal = length(y)
  knotsall = c(rep(boundary[1],4),knots, rep(boundary[2],4))
  obasisobj = OBasis(knotsall)
  Bm = evaluate(obasisobj,tm)  ## orthogonal

  
  ### ni for each subject ###
  uind = unique(ind)
  n = length(uind)
  nvec = rep(0,n)
  for(i in 1:n)
  {
    nvec[i] = sum(ind==uind[i])
  }
  
  p = ncol(Bm)

  np = n*p
  Ip = 1/lam*diag(p)
  nIp = lam*n*diag(p)
  
  DB = matrix(0, p, p)
  AB = matrix(0, np, p)
  matinv = matrix(0, np, np)
  Bty = rep(0,np)
  
  idp1 = 1
  idp2 = p
  for(i in 1:n)
  {
    Bmi = Bm[ind == uind[i],]
    ni = nrow(Bmi)
    
    mati = solve(t(Bmi)%*% Bmi + nIp)
    DB = DB + mati
    AB[idp1:idp2,] = mati
    matinv[idp1:idp2,idp1:idp2] = mati
    Bty[idp1:idp2] = t(Bmi) %*% y[ind == uind[i]]
    
    idp1 = idp1 + p
    idp2 = idp2 + p
  }
  
  IB = solve(Ip - DB)
  
  matinv = matinv + AB %*% IB %*% t(AB)
  
  
  betam0 = matrix(matinv%*%Bty,n,p,byrow= TRUE)
  
  return(betam0)
  
#### check the result with the original definition
  ### correct 
  # library(Matrix)
  # Xm <- matrix(0, ntotal, np)
  # for(i in 1:n)
  # {
  #   Xm[ind == uind[i],(p*(i-1) + 1) : (p*i)] = Bm[ind == uind[i],]
  # }
  # 
  # D <- matrix(0,n*(n-1)/2,n)
  # for(j in 1:(n-1))
  # {
  #   indexj <- (n-1 + n-j+1)*(j-1)/2
  #   indexvj <- indexj + (1:(n-j))
  #   D[indexvj,j] <- 1
  #   D[cbind(indexvj,(j+1):n)] <- -1
  # }
  # 
  # 
  # AtA <- t(D)%*%D %x% diag(p)
  # 
  # betam01 = solve(t(Xm)%*%Xm + lam*AtA)%*%t(Xm)%*%y
  
  
}


####### initial value based on Laplacian matrix for model with x grouping #####
## assume that the number of rows of x is the same as the length of y 

initiallapx = function(ind, tm, x, y, knots, boundary = c(0,1), lam = 0.01)
{
  ntotal = length(y)
  knotsall = c(rep(boundary[1],4),knots, rep(boundary[2],4))
  obasisobj = OBasis(knotsall)
  Bm = evaluate(obasisobj,tm)  ## orthogonal
  Bm = cbind(x, Bm)
  
  
  ### ni for each subject ###
  uind = unique(ind)
  n = length(uind)
  nvec = rep(0,n)
  for(i in 1:n)
  {
    nvec[i] = sum(ind==uind[i])
  }
  
  p = ncol(Bm)
  
  np = n*p
  Ip = 1/lam*diag(p)
  nIp = lam*n*diag(p)
  
  DB = matrix(0, p, p)
  AB = matrix(0, np, p)
  matinv = matrix(0, np, np)
  Bty = rep(0,np)
  
  idp1 = 1
  idp2 = p
  for(i in 1:n)
  {
    Bmi = Bm[ind == uind[i],,drop = FALSE]
    ni = nrow(Bmi)
    
    mati = solve(t(Bmi)%*% Bmi + nIp)
    DB = DB + mati
    AB[idp1:idp2,] = mati
    matinv[idp1:idp2,idp1:idp2] = mati
    Bty[idp1:idp2] = t(Bmi) %*% y[ind == uind[i]]
    
    idp1 = idp1 + p
    idp2 = idp2 + p
  }
  
  IB = solve(Ip - DB)
  
  matinv = matinv + AB %*% IB %*% t(AB)
  
  betam0 = matrix(matinv%*%Bty,n,p,byrow= TRUE)
  return(betam0)
}



##### input not knots and boundary, should be the design matrix ####

initiallap_mat = function(ind, y, xm, lam = 0.001)
{
  ntotal = length(y)
 
  
  ### ni for each subject ###
  uind = unique(ind)
  n = length(uind)
  nvec = rep(0,n)
  for(i in 1:n)
  {
    nvec[i] = sum(ind==uind[i])
  }
  
  p = ncol(xm)
  
  np = n*p
  Ip = 1/lam*diag(p)
  nIp = lam*n*diag(p)
  
  DB = matrix(0, p, p)
  AB = matrix(0, np, p)
  matinv = matrix(0, np, np)
  Bty = rep(0,np)
  
  idp1 = 1
  idp2 = p
  for(i in 1:n)
  {
    xmi = xm[ind == uind[i],]
    ni = nrow(xmi)
    
    mati = solve(t(xmi)%*% xmi + nIp)
    DB = DB + mati
    AB[idp1:idp2,] = mati
    matinv[idp1:idp2,idp1:idp2] = mati
    Bty[idp1:idp2] = t(xmi) %*% y[ind == uind[i]]
    
    idp1 = idp1 + p
    idp2 = idp2 + p
  }
  
  IB = solve(Ip - DB)
  
  matinv = matinv + AB %*% IB %*% t(AB)
  
  
  betam0 = matrix(matinv%*%Bty,n,p,byrow= TRUE)
  
  return(betam0)
}


### common intercept: doesn't work #####
### common eta?


library(Matrix)
est_initial = function(ind, tm, x, y, knots, boundary = c(0,1), lam = 0.001)
{
  ntotal = length(y)
  knotsall = c(rep(boundary[1],4),knots, rep(boundary[2],4))
  obasisobj = OBasis(knotsall)
  Bm = evaluate(obasisobj,tm)  ## orthogonal

  uind = unique(ind)
  n = length(uind)
  p = ncol(Bm)
  np = n * p

  Xm = matrix(0, ntotal, np)
  for(i in 1:n)
  {
    Xm[ind == uind[i],(p*(i-1) + 1) : (p*i)] = Bm[ind == uind[i],]
  }

  D = matrix(0,n*(n-1)/2,n)
  for(j in 1:(n-1))
  {
    indexj = (n-1 + n-j+1)*(j-1)/2
    indexvj = indexj + (1:(n-j))
    D[indexvj,j] = 1
    D[cbind(indexvj,(j+1):n)] = -1
  }

  AtA = t(D)%*%D %x% diag(p)
  
  xtxinv = solve(t(x)%*%x)

  Qx = diag(ntotal) - x%*%xtxinv%*%t(x)

  betam0 = solve(t(Xm)%*%Qx%*%Xm + lam*AtA)%*%t(Xm)%*%Qx%*%y
  eta0 = xtxinv%*%t(x)%*%(y - Xm %*%betam0)

  bmx = Xm%*%betam0

  outlist = list(betam = matrix(betam0,n,p,byrow = TRUE), eta0 = eta0, bmx = bmx)
  return(outlist)

}
