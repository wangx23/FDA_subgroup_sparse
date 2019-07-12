##### 100 simulations ####
#### in this simulation, the number of components is also selected based on BIC
####
library(flexclust)
library(doParallel)

subfunP = function(mm)
{
  source("simdat.R")
  source("initial.R")
  source("EMgroup.R")
  source("FDAsubgroup.R")
  source("BICvalue.R")
  
  dat = simdat(sig2 = 0.1,lamj = c(0.2,0.1),mvec = c(5,20),ncl = 50, seed = mm + 4452)
  #dat = read.csv("dat.csv")
  group = unique(dat[,1:2])[,1]
  
  lamvec = seq(0.1,1,by = 0.025)
  
  BICm = matrix(0,length(lamvec),3)
  for(j in 1:length(lamvec))
  {
    for(Pv in 1:3)
    {
      resi = try(FDAsubgroup(ind = dat$ind,tm = dat$time,y = dat$obs,P = Pv,
                             knots = seq(0,1,length.out = 6)[2:5],
                             lam = lamvec[j],maxiter = 50))
      errori = inherits(resi,"try-error")
      if(errori)
      {
        BICm[j, Pv] = 99999
      }else{
        BICm[j, Pv] = BICvalue(resi)
      }
    }
  }
  
  
  inds = which(BICm == min(BICm), arr.ind = TRUE)
  
  res = FDAsubgroup(ind = dat$ind,tm = dat$time,y = dat$obs,P = inds[2],
                    knots = seq(0,1,length.out = 6)[2:5],
                    lam = lamvec[inds[1]],maxiter = 50)
  
  sig2est = res$sig2
  ari = randIndex(table(group,res$groupest))
  
  return(c(ari, sig2est, inds[2]))
  
}


#temp = subfun(10)
cl <- makeCluster(28)  
registerDoParallel(cl)  
resulttestP1 <- foreach(mm=1:100,.packages=c("flexclust")) %dopar% subfunP(mm)
stopCluster(cl) 


save(resulttestP1, file = "resulttestP1.RData")




