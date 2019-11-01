######## EM based on known group information ###
## using this as initial values 

library(orthogonalsplinebasis)

## ind: is the subject ID
## tm is time
## y is measurement or observation
## group0 a vector is the known group information
## maxiter is the number of interations


## a function to keep the first element with the largest magnitude to be positive
# input is a vector 
max2pos = function(x)
{
  maxx = which.max(abs(x))
  x = x * sign(x[maxx])
  return(x)
  
}

EMgroup = function(ind, tm, y, knots, group0, P, betam0, boundary = c(0,1),
                   maxiter = 50, tol = 1e-3)
{
  
  ntotal = length(y)
  knotsall = c(rep(boundary[1],4),knots, rep(boundary[2],4))
  obasisobj = OBasis(knotsall)
  Bm = evaluate(obasisobj,tm)  ## orthogonal
  
  p = ncol(Bm)
  uind = unique(ind)
  n = length(uind)
  
  ng = length(unique(group0))
  Ux = matrix(0, ntotal, ng*p)
  
  BtB = array(0, dim = c(p,p,n))
  ylist = list()
  Blist = list()
  
  for(i in 1:n)
  {
    groupi = group0[i]
    indi = ind == uind[i]
    Bmi = Bm[indi,]
    Ux[indi,((groupi - 1)*p+1):(groupi*p)] = Bmi
    BtB[,,i] = t(Bmi) %*% Bmi
    ylist[[i]] = y[indi]
    Blist[[i]] = Bmi
  }
  
  est = solve(t(Ux)%*%Ux) %*% t(Ux) %*%y
  alpm = matrix(est, ng, p, byrow = TRUE)
  
  Cm = matrix(0, p, p)
  for(i in 1:n)
  {
    cv = betam0[i,] - alpm[group0[i],]
    Cm = Cm + cv %*% t(cv)/n
  }
  
  decomp = eigen(Cm)
  lamj = decomp$values[1:P]
  theta = decomp$vectors[,1:P,drop = FALSE]
  theta = apply(theta,2,max2pos)


  sig2 = mean((y - Ux %*%est)^2)
  
  
  mhat = matrix(0, n,P)
  Vhat = array(0, dim = c(P,P,n))
  Bty = matrix(0, n, p )
  thetacross = array(0, dim = c(n,p,P))
  residv = rep(0,n)
  Btheta = rep(0,ntotal)

  niteration = 0
  
  for(m in 1:maxiter)
  {
    lamjold = lamj 
    thetaold = theta
    sig2old = sig2
    alpmold = alpm
    
    Laminv = diag(1/lamj,P,P)
    ### E step 
    Sigma = matrix(0, P,P)
    E2mat = array(0, dim = c(p,p, P))
    
    for(i in 1:n)
    {
      Bmi = Blist[[i]]
      Vi = solve(t(theta) %*% BtB[,,i] %*% theta/sig2 + Laminv)
      Btyi = t(Bmi) %*% (ylist[[i]] - Bmi%*%alpm[group0[i],])
      mi = 1/sig2 * Vi %*% t(theta) %*% Btyi
      Sigma = Sigma + mi %*% t(mi) + Vi
      residv[i] =  sum((ylist[[i]] - Bmi%*%alpm[group0[i],] - Bmi %*%theta %*% mi)^2) + sum(diag(Bmi %*% theta %*% Vi %*% t(theta) %*% t(Bmi)))
      
      for(j in 1:P)
      {
        E2mat[,,j] = E2mat[,,j] + BtB[,,i] * (mi[j]^2 + Vi[j,j])
        thetacross[i,,j] = BtB[,,i] %*% theta[,j]
      }
      mhat[i,] = mi
      Vhat[,,i] = Vi
      Bty[i,] = Btyi
    }
    
    ## update theta ###
    
    for(j in 1:P)
    {
      temp = rep(0,p)
      for(j1 in (1:P)[(1:P)!=j])
      {
        temp = temp + t(thetacross[,,j1]) %*% (mhat[,j1] * mhat[,j] + Vhat[j1,j,])
      }
      invEj = solve(E2mat[,,j]) 
      BtyEj = t(Bty)%*% mhat[,j] - temp
      theta[,j] = invEj %*% BtyEj
    }
    
    sig2 = sum(residv)/ntotal
    
    # update theta and lamj
    Sigma = Sigma/n
    M0 = theta %*% Sigma %*%t(theta)
    decompm = eigen(M0, symmetric = TRUE)
    lamj = decompm$values[1:P]
    theta = decompm$vectors[,1:P,drop=FALSE]
    theta = apply(theta,2,max2pos)
    
    
    
    
    # update est 
    
    for(i in 1:n)
    {
      indi = ind == uind[i]
      Btheta[indi] = Blist[[i]] %*% theta %*% mhat[i,]
    }

    est = solve(t(Ux)%*%Ux) %*% t(Ux) %*%(y-Btheta)
    alpm = matrix(est, ng, p, byrow = TRUE)
    
    
    diffvalue = sqrt(sum((lamj - lamjold)^2) + sum((thetaold - theta)^2) + (sig2old - sig2)^2)
    niteration = niteration + 1
    if(diffvalue <= tol){break}
    
  }
  
  residf = rep(0,n)
  for(i in 1:n)
  {
    residf[i] = sum((ylist[[i]] - Blist[[i]] %*% alpm[group0[i],] -
                       Blist[[i]] %*% theta %*% mhat[i,])^2)
  }
  
  residsum = sum(residf)
  
  
  # res = list(sig2 = sig2, theta = theta, alpm = alpm, diffvalue = diffvalue,
  #        lamj = lamj, residsum = residsum, thetaold = thetaold, sig2old = sig2old,lamjold = lamjold, alpmold = alpmold,
  #        niteration = niteration)
  
  res = list(sig2 = sig2, theta = theta, alpm = alpm, diffvalue = diffvalue,
         lamj = lamj, residsum = residsum, niteration = niteration)
  
  return(res)

}
  





