source("BICvalue.R")
source("simdatx.R")
source("initial.R")
source("refitFDA.R")
source("FDAsubgroup.R")
library(flexclust)
library(doParallel)
library(cluster)

subfunx = function(mm, sig200, lam00, mvec00, ncl00,
                   lamvec, funlist, eigenlist, xlist,
                   method = "median", K0 = 10)
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
  
  if(method == "median")
  {
    betam021median = apply(betam002[,-1], 1, median)
    groupb2 = as.numeric(cut(betam021median,quantile(betam021median,seq(0,1,by=0.1)), include.lowest = TRUE))
  }
  
  if(method == "kmeans")
  {
    repeat{
      groupb2 = kmeans(betam002[,-1],centers = K0, iter.max = 20)$cluster
      if(min(table(groupb2))>1){break}
    }
  }
  
  if(method == "ydist")
  {
    ny = length(unique(datx$ind))
    distmaty = matrix(0, ny, ny)
    xm = unique(as.matrix(datx[,c("x1","x2")]))
    for(i in 1:(ny-1))
    {
      tsi = sort(datx$time[datx$ind==i],index.return = TRUE)
      tmi = tsi$x
      yi = datx$obs[datx$ind==i][tsi$ix]
      
      sfuni = stepfun(tmi, c(yi,yi[length(yi)]),f=1, right = TRUE)
      for(j in (i+1):ny)
      {
        tsj = sort(datx$time[datx$ind==j],index.return = TRUE)
        tmj = tsj$x
        yj = datx$obs[datx$ind==j][tsj$ix]
        
        sfunj = stepfun(tmj, c(yj,yj[length(yj)]),f=1, right = TRUE)
        
        funij = function(x)
        {
          (sfuni(x) - sfunj(x))^2
        }
        funij = Vectorize(funij)
        
        value1 = stats::integrate(funij,0,1, subdivisions = 500)$value
        valueij = sqrt(value1) + sqrt(sum((xm[i,] - xm[j,])^2))
        distmaty[i,j] = valueij
      }
    }
    
    distmaty = distmaty + t(distmaty)
    groupb2 = pam(distmaty, K0, diss = TRUE)$clustering
  }
  
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

