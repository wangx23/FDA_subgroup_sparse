##### ADMM based on known theta, sig2 and lamj 
### update beta and eta(common coefficients)
library(expm)

#### calculate (B^TOmegaB + nuA^TA)^-1 and B^TOmega^-1y
## Omgea = Sigma^-1 - Sigma^-1X(X^TOmegaX)^-1 X^TOmega^-1
## q is the dimension of xm
inverseB2 = function(Bmlist, xlist, ylist, theta, lamj, sig2, n, p, q,nu = 1)
{

  n0 = n*p;
  Ip = 1/nu * diag(p)
  nIp = nu * n * diag(p)
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
    xmi = xlist[[i]]
    ni = nrow(Bmi)
    Bthetai = Bmi %*% theta
    #Sigmainvi  = solve(sig2 * diag(ni) + Bthetai %*% diag(lamj) %*% t(Bthetai))
    Sigmainvi  = solve( diag(ni) + 1/sig2*Bthetai %*% diag(lamj) %*% t(Bthetai))
    Sigmainvisr = sqrtm(Sigmainvi)
    Bmi =  Sigmainvisr%*% Bmi
    xmi = Sigmainvisr %*% xmi
    yi = Sigmainvisr %*% ylist[[i]]
    
    mati = solve(t(Bmi)%*% Bmi + nIp)
    xAx = xAx + t(xmi) %*%xmi -  t(xmi) %*% Bmi %*% mati%*% t(Bmi) %*% xmi
    
    DB = DB + mati
    temp = mati %*% t(Bmi) %*% xmi
    Ax = Ax + temp
    AB[idp1:idp2,] = mati
    AX[idp1:idp2,] = temp
    matinv[idp1:idp2,idp1:idp2] = mati
    
    
    Bty[idp1:idp2] = t(Bmi) %*% yi
    xSx = xSx + t(xmi) %*%xmi
    BSx[idp1:idp2,] = t(Bmi) %*% xmi
    xSy = xSy + t(xmi) %*% yi
    
    idp1 = idp1 + p
    idp2 = idp2 + p
  }
  
  xAxinv = solve(xAx)
  Axt = t(Ax)
  temp1 = solve(Ip - DB - Ax %*% xAxinv %*% Axt);
  
  A1 = AB %*% temp1%*% t(AB)
  A2 = AX %*% xAxinv %*% Axt
  A3 = A2 %*% temp1
  A22 = A3 %*% t(A2)
  A12 = A3 %*% t(AB)
  
  matinv = matinv + AX %*% xAxinv %*% t(AX) + A1 + A22 + A12 + t(A12)

  # solve(xSx) %*% (xSy - t(BSx)%*%betam02)
  # Bty - BSx %*% solve(xSx)%*% xSy
  res = list(matinv = matinv, Bty = Bty,
             BSx = BSx, xSxinv = solve(xSx),
             xSy = xSy)
  return(res)
  
}

#### check the result with the original definition
### correct 
# library(Matrix)
# Xm <- matrix(0, ntotal, n0)
# SigmaInv = matrix(0, ntotal, ntotal)
# for(i in 1:n)
# {
#   Xm[ind == uind[i],(p*(i-1) + 1) : (p*i)] = Bm[ind == uind[i],]
#   Bmi = Bmlist[[i]]
#   xmi = xlist[[i]]
#   ni = nrow(Bmi)
#   Bthetai = Bmi %*% theta
#   SigmaInv[ind == uind[i],ind == uind[i]] = solve( diag(ni) + 1/sig2*Bthetai %*% diag(lamj) %*% t(Bthetai))
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
# AtA <- t(D)%*%D %x% diag(p)
# 
# Omega = SigmaInv - SigmaInv %*% x %*% solve(t(x) %*% SigmaInv %*% x) %*% t(x) %*% SigmaInv
# betam01 = solve(t(Xm) %*% Omega %*% Xm + nu * AtA)%*% t(Xm) %*% Omega %*%y
# eta01 = solve(t(x) %*% SigmaInv %*% x) %*% t(x) %*%SigmaInv %*% (y - Xm %*%betam01)
# 
# 
# betam02 = matinv%*%(Bty - BSx %*% xSxinv%*% xSy)
# max(abs(betam02 - betam01))

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
revisedadmmv2 = function(Bmlist, xlist, ylist, n, p, q, npair,
                       Dmat,deltam, vm, lam, nu = 1, gam = 3,
                       theta, lamj, sig2)
{
  
  ### update beta
  resinv = inverseB2(Bmlist,xlist, ylist, theta, lamj, sig2, n, p,q, nu = nu)
  matinv = resinv$matinv
  Bty = resinv$Bty
  BSx = resinv$BSx
  xSxinv = resinv$xSxinv
  xSy = resinv$xSy
  
  
  betanew = matinv%*%(Bty - BSx %*% xSxinv%*% xSy  + c((deltam -  vm/nu) %*% Dmat))
  etanew = xSxinv %*% (xSy - t(BSx)%*%betanew)
  
  betam = matrix(betanew, n, p, byrow = TRUE)
  betadiff = t(Dmat %*% betam)
  deltam = betadiff + (1/nu) * vm
  
  for(i in 1:npair)
  {
    deltam[,i] = scad(deltam[,i],lam,nu,gam)
  }
  
  vm =  vm + nu * (betadiff - deltam)
  
  res = list(betam = betam, eta = etanew, deltam = deltam, vm = vm, betadiff = betadiff)
  return(res)
}