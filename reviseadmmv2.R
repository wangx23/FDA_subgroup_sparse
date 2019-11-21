##### ADMM based on known theta, sig2 and lamj

#### calculate (B^TSigma^-1B + nuA^TA)^-1 and B^TSigma^-1y
inverseB = function(Bmlist, ylist, theta, lamj, sig2, n, p, nu = 1)
{
  np = n*p
  Ip = 1/nu*diag(p)
  nIp = nu*n*diag(p)
  
  DB = matrix(0, p, p)
  AB = matrix(0, np, p)
  matinv = matrix(0, np, np)
  Bty = rep(0,np)
  
  idp1 = 1
  idp2 = p
  for(i in 1:n)
  {
    Bmi = Bmlist[[i]]
    ni = nrow(Bmi)
    Bthetai = Bmi %*% theta
    #Sigmainvi  = solve(sig2 * diag(ni) + Bthetai %*% diag(lamj) %*% t(Bthetai))
    Sigmainvi  = solve( diag(ni) + 1/sig2*Bthetai %*% diag(lamj) %*% t(Bthetai))
    mati = solve(t(Bmi) %*% Sigmainvi %*% Bmi + nIp)
    DB = DB + mati
    AB[idp1:idp2,] = mati
    matinv[idp1:idp2,idp1:idp2] = mati
    Bty[idp1:idp2] = t(Bmi) %*% Sigmainvi %*% ylist[[i]]
    
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