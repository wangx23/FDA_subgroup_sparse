##### 100 simulations ####
library(flexclust)
library(doParallel)

subfun = function(mm)
{
  source("simdat.R")
  source("initial.R")
  source("EMgroup.R")
  source("FDAsubgroup.R")
  source("BICvalue.R")
  
  dat = simdat(sig2 = 0.1,lamj = c(0.2,0.1),mvec = c(5,20),ncl = 50, seed = mm + 4452)
  group = unique(dat[,1:2])[,1]
  
  lamvec = seq(0.2,1,by = 0.05)
  
  BICvec = rep(0,length(lamvec))
  for(j in 1:length(lamvec))
  {
    resi = FDAsubgroup(ind = dat$ind,tm = dat$time,y = dat$obs,P = 2,
                       knots = seq(0,1,length.out = 6)[2:5],
                       lam = lamvec[j],maxiter = 50)
    BICvec[j] = BICvalue(resi)
    
  }
  
  res = FDAsubgroup(ind = dat$ind,tm = dat$time,y = dat$obs,P = 2,
                    knots = seq(0,1,length.out = 6)[2:5],
                    lam = lamvec[which.min(BICvec)],maxiter = 100)
  
  sig2est = res$sig2
  ari = randIndex(table(group,res$groupest))
  
  return(c(ari, sig2est))
  
}


cl <- makeCluster(28)  
registerDoParallel(cl)  
resulttest1 <- foreach(mm=1:100,.packages=c("flexclust")) %dopar% subfun(mm)
stopCluster(cl) 