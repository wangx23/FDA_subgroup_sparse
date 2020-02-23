source("BICvalue.R")
source("simdatx.R")
source("initial.R")
source("refitFDA.R")
source("FDAsubgroup.R")
source("initialgroup.R")
library(flexclust)
library(cluster)


# funlist = list(Vectorize(function(x){cos(2*pi*x)}),
#                 Vectorize(function(x){1-2*exp(-6*x)}),
#                 Vectorize(function(x){-1.5*x})
# )
# 
# eigenlist = list(Vectorize(function(x){sqrt(2)*sin(pi*x)}),
#                  Vectorize(function(x){ sqrt(2)*cos(pi*x)}))
# 
# xlist = list(nx = 2, meanx = 0, sdx = 1, etag = matrix(c(-0.5,-0.5,0.5,0.5,1,1),nrow = 2))
# sig200 = 0.04
# lamj00 = c(0.2,0.1)
# lamj00 = c(0.15, 0.1)
# mvec00 = c(10,20)
# lamvec = seq(0.2,0.8,by = 0.025)
# 

subfunx_all = function(mm, sig200, lam00, mvec00, ncl00,
                       lamvec, funlist, eigenlist, xlist, 
                       K0 = 15, ming = 1, max.step)
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
  
  betam002 = initiallap_mat(ind = datx$ind,y = datx$obs,xm = cbind(1, Bm[,-1]), lam = 0.0001) 
  
  # if(method == "median")
  #{
  betam021median = apply(betam002[,-1], 1, median)
  groupb21 = as.numeric(cut(betam021median,quantile(betam021median,seq(0,1,length.out = K0+1)), include.lowest = TRUE))
  #}
  
  # if(method == "kmeans")
  #{
  groupb22 = initialgroup_ks(betam002[,-1],K0 = K0, ming = ming)
  
  #}
  
  #if(method == "pam")
  #{
  groupb23 = initialgroup_pam(betam002[,-1],K0 = K0, diss = FALSE, ming = ming)
  #}
  
  #if(method == "ydist")
  #{
  ny = length(unique(datx$ind))
  distmaty = matrix(0, ny, ny)
  xm = unique(as.matrix(datx[,c("x1","x2")]))
  t1 = Sys.time()
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
  t2 = Sys.time()
  
  distmaty = distmaty + t(distmaty)
  #groupb24 = pam(distmaty, K0, diss = TRUE)$clustering
  groupb24 = initialgroup_pam(distmaty,K0 = K0, diss = TRUE)
  # }
  
  
  subestfun = function(groupbb)
  {
    
    nllmat = matrix(0, length(lamvec),3)
    arim = matrix(0, length(lamvec),3)
    ngest = matrix(0, length(lamvec),3)
    groupest_mat = array(0, dim = c(length(group0),length(lamvec),3))
    lamest1 = lamest2 = lamest3 = matrix(0, length(lamvec),3)
    
    #betam022 = refitINDX(ind = datx$ind, tm = datx$time, x=x, y = datx$obs, group0 = groupbb, knots = knots)
    
    betam022 = refitFDAX(ind = datx$ind,tm = datx$time,x = x,y = datx$obs,group0 = groupbb, knots = knots)$alpha[groupbb,]
    
    t1 = Sys.time()
    for(Pv in 1:3)
    {
      for(j in 1:length(lamvec))
      {
        resi = try(FDAXsubgroup(ind = datx$ind,x = x, tm = datx$time,y = datx$obs,P = Pv,
                                betam0 = betam022,group0 = groupbb,knots = knots, K0 = K0,
                                max.step = max.step,
                                lam = lamvec[j],maxiter = 50,tolabs = 1e-4,tolrel = 1e-2))
        errori = inherits(resi,"try-error")
        if(errori)
        {
          nllmat[j,Pv] = 9999
        }else{
          nllmat[j,Pv] = resi$likevalue
          arim[j, Pv] = randIndex(resi$groupest,group0)
          lamest1[j,Pv]= resi$lamj[1]
          lamest2[j,Pv] = resi$lamj[2]
          lamest3[j,Pv] = resi$lamj[3]
          ngest[j,Pv] = length(unique(resi$groupest))
          groupest_mat[,j,Pv] =  resi$groupest
        }
      }
    }
    t2 = Sys.time()
    
    
    outlist = list(nllmat = nllmat, arim = arim, lamest1= lamest1, lamest2 = lamest2,
                   lamest3 = lamest3, ngest = ngest, groupest = groupest_mat)
    return(outlist)
  }
  
  
  res1 = subestfun(groupbb = groupb21)
  res2 = subestfun(groupbb = groupb22)
  res3 = subestfun(groupbb = groupb23)
  res4 = subestfun(groupbb = groupb24)
  
  
  output = list(res_median = res1, res_kmeans = res2, res_pam = res3,
                res_ydist = res4, p = ncol(Bm) + ncol(x), ntotal = length(datx$obs))
  return(output)
}