source("BICvalue.R")
source("simdatx.R")
source("initial.R")
source("refitFDA.R")
source("FDAsubgroup.R")
library(flexclust)
library(doParallel)

funlist = list(Vectorize(function(x){4*(x/50)^(0.2)}),
               Vectorize(function(x){-4*(x/50)^(0.2)}))

eigenlist = list(Vectorize(function(x){sqrt(2)*sin(pi*x)}),
                 Vectorize(function(x){ sqrt(2)*cos(pi*x)}))

xlist = list(nx = 2, meanx = 0, sdx = 1, etag = matrix(c(-1,-1,1,1),nrow = 2))
sig200 = 0.04
lamj00 = c(0.2,0.1)
mvec00 = c(10,20)
ncl00 = 100
lamvec = seq(0.2,2,by = 0.1)

subfunx = function(mm, sig200, lam00, mvec00, ncl00,
                   lamvec, funlist, eigenlist, xlist)
{
  
  datx = simdatx(xlist = xlist,
                 sig2 = sig200,lamj = lamj00,mvec = mvec00,ncl = ncl00,
                 funlist = funlist, eigenlist = eigenlist, seed = mm+ 4452)
  
  group0 = rep(1:length(funlist),each = ncl00)
  ng0 = length(unique(group0))
  
  
  x = as.matrix(datx[,c("x1","x2")])
  
  ntotal = length(datx$obs)
  knots = seq(0,1,length.out = 5)[2:4]
  knotsall = c(rep(0,4),knots, rep(1,4))
  obasisobj = OBasis(knotsall)
  Bm = evaluate(obasisobj,datx$time)  ## orthogonal and include intercept 
  
  betam002 = initiallap_mat(ind = datx$ind,y = datx$obs,xm = cbind(1, Bm[,-1]), lam = 0.001) 
  
  betam021median = apply(betam002[,-1], 1, median)
  groupb2 = as.numeric(cut(betam021median,quantile(betam021median,seq(0,1,by=0.1)), include.lowest = TRUE))
  res022 = refitFDAX(ind = datx$ind,tm = datx$time,x = x,y = datx$obs,group0 = groupb2, knots = knots)
  betam022 = res022$alpha[groupb2,]
  
  
  BICm = BICm2 = BICm3 = matrix(0,length(lamvec), 3)
  arim = matrix(0, length(lamvec),3)
  
  t1 = Sys.time()
  for(Pv in 1:3)
  {
    for(j in 1:length(lamvec))
    {
      resi = try(FDAXsubgroup(ind = datx$ind,x = x, tm = datx$time,y = datx$obs,P = Pv,
                              betam0 = betam022,group0 = groupb2,knots = knots, K0 = 10,
                              max.step = 5,
                              lam = lamvec[j],maxiter = 50,tolabs = 1e-4,tolrel = 1e-2))
      errori = inherits(resi,"try-error")
      if(errori)
      {
        BICm[j, Pv] = 99999
        BICm2[j, Pv] = 99999
        BICm3[j, Pv] = 99999
      }else{
        BICm[j, Pv] = BICvalue(resi)
        BICm2[j, Pv] = BICvalue2(resi)
        BICm3[j, Pv] = BICvalue3(resi)
        arim[j, Pv] = randIndex(resi$groupest,group0)
      }
    }
  }
  t2 = Sys.time()
  
  inds = which(BICm == min(BICm), arr.ind = TRUE)
  
  # res1 = FDAXsubgroup(ind = datx$ind,x = x, tm = datx$time,y = datx$obs,P = 3,
  #              betam0 = betam022,group0 = groupb2,knots = knots, K0 = 10,max.step = 4,
  #              lam = lamvec[9],maxiter = 50,tolabs = 1e-4,tolrel = 1e-2)
  # 
  # res2 = FDAXsubgroup(ind = datx$ind,x = x, tm = datx$time,y = datx$obs,P = 2,
  #                    betam0 = betam022,group0 = groupb2,knots = knots, K0 = 10,max.step = 4,
  #                    lam = lamvec[9],maxiter = 50,tolabs = 1e-4,tolrel = 1e-2)
  
  # BICvalue2(res2)
  
  res = FDAXsubgroup(ind = datx$ind,x = x, tm = datx$time,y = datx$obs,P = 2,
                     betam0 = betam022,group0 = groupb2,knots = knots, K0 = 10,max.step = 4,
                     lam = lamvec[which.min(BICm[,2])],maxiter = 50,tolabs = 1e-4,tolrel = 1e-2)
  
  
  ari = randIndex(res$groupest,group0)
  lamest = res$lamj
  Pest = c(which(BICm == min(BICm), arr.ind = TRUE)[2],
           which(BICm2 == min(BICm2), arr.ind = TRUE)[2],
           which(BICm3 == min(BICm3), arr.ind = TRUE)[2])
  
  aribest = max(arim)
  indaribest = which(arim == max(ari), arr.ind = TRUE)
  
  output = list(ari = ari, lamest = lamest, Pest = Pest,
                aribest = aribest, indaribest = indaribest)
  return(output)
}


# t1 = Sys.time()
# restest = subfunx(1, sig200 = sig200, lam00 = lamj00, mvec00 = mvec00,
#                   ncl00 = ncl00, lamvec = lamvec, 
#                   funlist = funlist, eigenlist = eigenlist,xlist= xlist)
# t2 = Sys.time()

subfunx2cluster = function(mm){
  subfunx(mm, sig200 = sig200, lam00 = lamj00, mvec00 = mvec00,
          ncl00 = ncl00, lamvec = lamvec, 
          funlist = funlist, eigenlist = eigenlist,xlist= xlist)
}

cl <- makeCluster(2)  
registerDoParallel(cl)  
result_x_2cluster <- foreach(mm=1:2,
                       .packages=c("flexclust","orthogonalsplinebasis","plyr"),
                       .errorhandling = "remove") %dopar%subfunx2cluster(mm)
stopCluster(cl) 
save(result_x_2cluster,file = "../result/result_x_2cluster.RData")