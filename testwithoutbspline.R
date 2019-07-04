beta1 = c(-2,-2)
beta2 = c(2,2)
betam12 = rbind(beta1,beta2)

xm = cbind(1,rnorm(length(group)*20))
ind1 = rep(1:length(group), each = 20)
y1 = rep(0, length(group)*20)
for(i in 1:length(group))
{
  indi = ind1 == i
  y1[indi] = xm[indi,] %*% betam12[group[i],] + rnorm(20)*0.5
}


betam = matrix(0,length(group),2)
for(i in 1:length(group))
{
  indi = ind1 ==i
  betam[i,] = coef(lm(y1[indi]~0+xm[indi,]))
}



Bmlist = list()
ylist = list()
for(i in 1:n)
{
  indi = ind1 == i
  Bmlist[[i]] = xm[indi,]
  ylist[[i]] = y1[indi]
}


mat1 = matrix(0,2*n, 2*n)
Bty1 = matrix(0, 2*n)
for(i in 1:n)
{
  mat1[((i-1)*2+1):(i*2),((i-1)*2+1):(i*2)] = t(Bmlist[[i]]) %*% Bmlist[[i]]/sig2
  Bty1[((i-1)*2+1):(i*2)]  = 
}


Am = kronecker(Dmat, diag(2))
matinv1 = solve(mat1 + nu * t(Am) %*% Am)

max(abs(matinv1 - matinv))


res2 = concavefusionSpatial_v4_x(ind1, y1, xm, rep(1,npair),
                                      sol.fun = scad, nu = 1, gam = 3, lam = 0.5,
                                      betam,
                                      maxiter = 500, tol = 1e-5)
