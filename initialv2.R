##### initial values including the common covariates #####


initialv2 = function(ind, x,  tm, y, knots, boundary = c(0,1), lam = 0.001)
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
  q = ncol(x)
  
  
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
  
  mati = matrix(0,p,p);
  
  

  idp1 = 1
  idp2 = p
  
  Hy = y - x %*% solve(t(x) %*%x) %*%t(x)%*%y
  
  for(i in 1:n)
  {
    indexi = ind == uind[i]
    Bmi = Bm[indexi,]
    xi = x[indexi,]
    ni = nrow(Bmi)
    mati = solve(t(Bmi)%*% Bmi + nIp)
    xAx = xAx + t(xi) %*%xi -  t(xi) %*% Bmi %*% mati%*% t(Bmi) %*% xi
    
    DB = DB + mati
    
    temp = mati %*% t(Bmi) %*% xi
    Ax = Ax + temp
    AB[idp1:idp2,] = mati
    AX[idp1:idp2,] = temp
    matinv[idp1:idp2,idp1:idp2] = mati
    Bty[idp1:idp2] = t(Bmi) %*% Hy[indexi]
    
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
  

  betam0 = matrix(matinv%*%Bty,n,p,byrow= TRUE)
  
  
  Bbeta = rep(0, ntotal)
  for(i in 1:n)
  {
    Bbeta[ind == uind[i]] = Bm[ind == uind[i],] %*% betam0[i,]
  }
  
  eta0 = solve(t(x) %*% x) %*% t(x) %*%(y - Bbeta)
  
  
  initlist = list(betam = betam0, eta = eta0)
  return(initlist)
  
  #### check the result with the original definition
  ### correct 
  # library(Matrix)
  # Xm <- matrix(0, ntotal, n0)
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
  # # 
  # # 
  # AtA <- t(D)%*%D %x% diag(p)
  # # 
  # betam01 = solve(t(Xm)%*% Xm - t(Xm) %*% x %*% solve(t(x) %*%x) %*%t(x) %*%Xm + lam*AtA)%*%t(Xm)%*%Hy
  
  #eta01 = solve(t(x) %*% x) %*% t(x) %*%(y - Xm %*%betam01)
  
  
}



