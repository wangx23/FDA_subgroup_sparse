###### parallel test ##### 
library(flexclust)
library(funcy)
library(doParallel)

source("simdat.R")
source("initial.R")
source("EMgroup.R")
source("FDAsubgroup.R")
source("BICvalue.R")

library(plyr);
library(fda);#bspline basis
library(Matrix);
library(matrixcalc);
library(igraph);
library(Rcpp);
library(RcppArmadillo);
library(Lclust)


#sourceCpp("../LongitudinalClustering-master/admmmcp_code.cpp")
#source("../LongitudinalClustering-master/clustering_functions.R")
#source("../LongitudinalClustering-master/BICfunctions.R")

subfun = function(mm)
{
  dat = simdat(sig2 = 0.1,lamj = c(0.2,0.1),mvec = c(5,20),ncl = 50,seed = mm + 4452)
  group0 = rep(1:2,each = 50)

  knots1 = seq(0,1,length.out = 6)[2:5]
  betam0 = initialcoef(ind = dat$ind,tm = dat$time,y = dat$obs,knots = knots1,
                       lamv = seq(0,20,by = 0.5)[-1])
  
  resem = EMgroup(ind = dat$ind,tm = dat$time,y = dat$obs,
                  knots =knots1 ,group0 = group0,P = 2,betam0 = betam0)
  
  alpor = resem$alpm
  betaor = alpor[group0,]
  
  #### proposed algorithm ###
  lamvec = seq(0.2,1,by = 0.05)
  BICm = matrix(0,length(lamvec), 3)
  
  t1 = Sys.time()
  for(Pv in 1:3)
  {
    for(j in 1:length(lamvec))
    {
      resi = try(FDAsubgroup(ind = dat$ind,tm = dat$time,y = dat$obs,P = Pv,
                         betam0 = betam0, knots = knots1,
                         lam = lamvec[j],maxiter = 50,tolabs = 1e-4,tolrel = 1e-2))
      errori = inherits(resi,"try-error")
      if(errori)
      {
        BICm[j, Pv] = 99999
      }else{
        BICm[j, Pv] = BICvalue(resi)
      }
    }
  }
  t2 = Sys.time()
  
  
  inds = which(BICm == min(BICm), arr.ind = TRUE)
  
  res = FDAsubgroup(ind = dat$ind,tm = dat$time,y = dat$obs,P = inds[2],
                    betam0 = betam0, knots = seq(0,1,length.out = 6)[2:5],
                    lam = lamvec[inds[1]],maxiter = 50,tolabs = 1e-4,tolrel = 1e-2)
  
  ari1 = randIndex(res$groupest, group0)
  meanest1 = res$meanfunest
  rmse1 = sqrt(mean((meanest1 - dat$mean)^2))
  ng1 = length(unique(res$groupest))

  
  
  ####### without covariance structure #####
  nknots = 4;
  order = 4;
  p =  order + nknots;
  #time points 
  TT = max(table(dat$ind))
  timerange = seq(0, 1, length.out = TT);
  
  n = 100
  
  basis = dlply(dat, .(ind), function(xx) bsplineS(xx$time, knots_eq3(timerange, k = order, m = nknots), norder = order))
  X = bdiag(basis)
  X = as.matrix(X)
  
  C <- matrix(0, nrow=nknots+order-2, ncol=nknots+order)
  for (j in 1:(nknots+order-2)){
    d_j <- c(rep(0,j-1),1,-2,1,rep(0,(nknots+order)-3-(j-1)))
    e_j <- c(rep(0,j-1), 1 ,rep(0,(nknots+order)-3-(j-1)))
    C <- C + e_j%*%t(d_j)
  }
  
  D = t(C)%*%C;
  diagD <- kronecker(diag(1, n), D);
  
  
  
  gamma1 = 0.005
  index = t(combn(n,2));
  B_ini0 = update_B_ini(X, diagD, as.vector(dat$obs), n, gamma1, index, lambda0 = gamma1)
  
  
  t3 = Sys.time()
  lamvec2 = seq(0.5,2,by=0.1)
  
  BIC2vec = rep(0, length(lamvec2))
  
  for(j in 1:length(lamvec2))
  {
    sol_finalj = prclust_admm(X, y=as.vector(dat$obs), diagD, B_ini0, index,
                              gamma1 = 0.005, gamma2 = lamvec2[j], 
                              theta=1, tau = 2, 
                              n, p,  max_iter= 200,
                              eps_abs=1e-4, eps_rel=1e-2)
    BIC2vec[j] = BIC2(obj = sol_finalj,dat$obs, dat$ind, X, basis, n, 0.005, D)
    
  }
  t4 = Sys.time()
  
  
  inds2 = which.min(BIC2vec)
  
  sol_final = prclust_admm(X, y=as.vector(dat$obs), diagD, B_ini0, index,
                           gamma1 = 0.005, gamma2 = lamvec2[inds2], 
                           theta=1, tau = 2, n, p,  max_iter=200,
                           eps_abs=1e-4, eps_rel=1e-2)
  #result:
  Ad_final <- create_adjacency(sol_final$V, n);
  G_final <- graph.adjacency(Ad_final, mode = 'upper')
  #clustering membership
  cls_final <- components(G_final);
  #number of clusters
  k_final <- cls_final$no
  
  meanest2 = X %*% c(sol_final$B)
  
  
  ari2 = randIndex(cls_final$membership, group0)
  rmse2 = sqrt(mean((meanest2 - dat$mean)^2))
  
  #### James and Sugar(2003)
  
  resjs = funcit(data = as.matrix(dat[,c("ind","obs","time")]), 
                 methods = "fitfclust",k = 2)
  groupjs = resjs@allClusters
  
  
  time1 = difftime(t2, t1, units = "mins")
  time2 = difftime(t4, t3, units = "mins")
  
  result = rbind(c(ari1, rmse1, time1,lamvec[inds[1]], ng1, inds[2]), 
                 c(ari2, rmse2, time2, lamvec2[inds2], k_final,0))
  colnames(result) = c("ARI","RMSE","time(mins)","lam","numgroup","numP")
  return(result)
  
}


cl <- makeCluster(28)  
registerDoParallel(cl)  
resultcomp1 <- foreach(mm=1:100,
                       .packages=c("flexclust","orthogonalsplinebasis","plyr","fda","Lclust","igraph")) %dopar%  subfun(mm)
stopCluster(cl) 
save(resultcomp1,file = "result/resultcomp1.RData")

# res1 = subfun(1)
# res10 = subfun(10)
# 
# res20 = subfun(20)


