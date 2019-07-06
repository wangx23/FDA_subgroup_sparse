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
  #dat = read.csv("dat.csv")
  group = unique(dat[,1:2])[,1]
  
  lamvec = seq(0.1,1,by = 0.025)
  
  BICvec = rep(0,length(lamvec))
  for(j in 1:length(lamvec))
  {
    resi = try(FDAsubgroup(ind = dat$ind,tm = dat$time,y = dat$obs,P = 2,
                       knots = seq(0,1,length.out = 6)[2:5],
                       lam = lamvec[j],maxiter = 50))
    errori = inherits(resi,"try-error")
    if(errori)
    {
      BICvec[j] = 99999
    }else{
      BICvec[j] = BICvalue(resi)
    }
    
  }
  
  res = FDAsubgroup(ind = dat$ind,tm = dat$time,y = dat$obs,P = 2,
                    knots = seq(0,1,length.out = 6)[2:5],
                    lam = lamvec[which.min(BICvec)],maxiter = 50)
  
  sig2est = res$sig2
  ari = randIndex(table(group,res$groupest))
  
  return(c(ari, sig2est))
  
}


#temp = subfun(10)
cl <- makeCluster(28)  
registerDoParallel(cl)  
resulttest1 <- foreach(mm=1:100,.packages=c("flexclust")) %dopar% subfun(mm)
stopCluster(cl) 


save(resulttest1, file = "resulttest1.RData")




