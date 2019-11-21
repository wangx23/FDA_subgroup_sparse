##### ADMM based on known theta, sig2 and lamj 
### update beta and eta(common coefficients)
library(expm)

#### calculate (B^TOmegaB + nuA^TA)^-1 and B^TOmega^-1y
## Omgea = Sigma^-1 - Sigma^-1X(X^TOmegaX)^-1 X^TOmega^-1
## q is the dimension of xm
inverseB2 = function(Bmlist, xmlist, ylist, theta, lamj, sig2, n, p, q,nu = 1)
{

  n0 = n*p;
  Ip = 1/lam * diag(p)
  nIp = lam * n * diag(p)
  matinv = matrix(0, n0,n0)
  Bty = rep(0,n0)
  DB = matrix(0,p,p)
  xAx = matrix(0,q, q);
  Ax = matrix(0,p,q);
  AB = matrix(0,n0, p);
  AX = matrix(0,n0,q);
  
  xSx = matrix(0, q, q)
  BSx = matrix(0, n0, q)
  xSy = rep(0,q)
  
  mati = matrix(0,p,p);

  idp1 = 1
  idp2 = p
  for(i in 1:n)
  {
    Bmi = Bmlist[[i]]
    xmi = xmlist[[i]]
    ni = nrow(Bmi)
    Bthetai = Bmi %*% theta
    #Sigmainvi  = solve(sig2 * diag(ni) + Bthetai %*% diag(lamj) %*% t(Bthetai))
    Sigmainvi  = solve( diag(ni) + 1/sig2*Bthetai %*% diag(lamj) %*% t(Bthetai))
    Sigmainvisr = sqrtm(Sigmainvisr)
    Bmi =  Sigmainvisr%*% Bmi
    xmi = Sigmainvisr %*% xmi
    yi = Sigmainvisr %*% ylist[[i]]
    
    mati = solve(t(Bmi)%*% Bmi + nIp)
    xAx = xAx + t(xmi) %*%xi -  t(xmi) %*% Bmi %*% mati%*% t(Bmi) %*% xmi
    
    temp = mati %*% t(Bmi) %*% xmi
    Ax = Ax + temp
    AB[idp1:idp2,] = mati
    AX[idp1:idp2,] = temp
    matinv[idp1:idp2,idp1:idp2] = mati
    Bty[idp1:idp2] = t(Bmi) %*% yi
    xSx = xSx + t(xmi) %*%xmi
    BSx[ind1:ind2,] = t(Bmi) %*% xmi
    xSy = xSy + t(xmi) %*% yi
    
    idp1 = idp1 + p
    idp2 = idp2 + p
  }
  
  IB = solve(Ip - DB)
  
  matinv = matinv + AB %*% IB %*% t(AB)
  
  res = list(matinv, Bty)
  return(res)
  
}


####SCAD ####
sfun <- function(x, th)
{
  xn <- sqrt(sum(x^2))
  thval <- 1 - th/xn
  thval*((thval) >0)*x
}
scad <- function(x,lam,nu,gam)
{
  temp1 <- lam/nu
  temp2 <- gam * lam
  xn <- sqrt(sum(x^2))
  
  if(xn <= lam + temp1)
  {
    z <- sfun(x, temp1)
  }else if(xn <= temp2 & xn >= lam + temp1)
  {
    z <- sfun(x, temp2/((gam-1)*nu))/(1- 1/((gam - 1 )*nu))
  }else{
    z <- x
  }
  
  return(z)
}





##### ADMM algorithm ####
revisedadmm = function(Bmlist, ylist, n, p, npair,
                       Dmat,deltam, vm, lam, nu = 1, gam = 3,
                       theta, lamj, sig2)
{
  
  ### update beta
  resinv = inverseB(Bmlist, ylist, theta, lamj, sig2, n, p, nu = nu)
  matinv = resinv[[1]]
  Bty = resinv[[2]]
  
  betanew = matinv%*%(Bty + c((deltam -  vm/nu) %*% Dmat))
  betam = matrix(betanew, n, p, byrow = TRUE)
  betadiff = t(Dmat %*% betam)
  deltam = betadiff + (1/nu) * vm
  
  for(i in 1:npair)
  {
    deltam[,i] = scad(deltam[,i],lam,nu,gam)
  }
  
  vm =  vm + nu * (betadiff - deltam)
  
  res = list(betam = betam, deltam = deltam, vm = vm, betadiff = betadiff)
  return(res)
}