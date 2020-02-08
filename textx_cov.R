#### test for examples with x involved
source("simdatx.R")
source("initial.R")
source("refitFDA.R")
source("FDAsubgroup.R")
library(flexclust)


# funlist21 = list(Vectorize(function(x){4*(x-0.5)^2 + 1}),
#                Vectorize(function(x){2.5*exp(-25*(x-0.25)^2) + 2*exp(-50*(x-0.75)^2)}))

funlist21 = list(Vectorize(function(x){4*(x/50)^(0.2)}),
                 Vectorize(function(x){-4*(x/50)^(0.2)}))

eigenlist21 = list(Vectorize(function(x){sqrt(2)*sin(pi*x)}),
                   Vectorize(function(x){ sqrt(2)*cos(pi*x)}))

xlist = list(nx = 2, meanx = 0, sdx = 1, etag = matrix(c(-1,-1,1,1),nrow = 2))
sig200 = 0.04
lamj00 = c(0.2,0.1)

mm = 1
datx = simdatx(xlist = xlist,
               sig2 = sig200,lamj = lamj00,mvec = c(10,20),ncl = 100,
               funlist = funlist21, eigenlist = eigenlist21, seed = 1 + 4452)

ind = datx$ind
tm = datx$time
x = as.matrix(datx[,c("x1","x2")])
y = datx$obs
knots = seq(0,1,length.out = 5)[2:4]

group00 = unique(datx[,c("group","ind")])[,1] ## true group 


res00 = refitFDAX(ind = datx$ind,tm = datx$time,x = x, y = datx$obs,P = 2,group0 = group00, 
                  knots = knots)
betatd00 = res00$alpha[group00,] ### true eta and beta 



### try 1 
betatd01 = initiallapx(ind = datx$ind,x = x,tm = datx$time,
                       y = datx$obs,knots = knots,lam = 0.01)

res1= FDAXsubgroup(ind = datx$ind,x = x, tm = datx$time,y = datx$obs,P = 2,
                   betam0 = betatd00,knots = knots, K0 = 10,
                   lam = 0.25,maxiter = 50,tolabs = 1e-4,tolrel = 1e-2)
plot(res1$betam)

### if put true plus randomness as initial ?

betam02 = betatd00 + 1*rnorm(nrow(betatd00)*ncol(betatd00))


res2= FDAXsubgroup(ind = datx$ind,x = x, tm = datx$time,y = datx$obs,P = 3,
                   betam0 = betam02,knots = knots, K0 = 10, max.step = 5,
                   lam = 1,maxiter = 50,tolabs = 1e-4,tolrel = 1e-2)
plot(res2$betam)
res2$lamj






### note: traditional method does not work

########### no covariance#######
### cal_initialrx
library(Spgr)
library(orthogonalsplinebasis)

ntotal = length(datx$obs)
knots = seq(0,1,length.out = 5)[2:4]
knotsall = c(rep(0,4),knots, rep(1,4))
obasisobj = OBasis(knotsall)
Bm = evaluate(obasisobj,datx$time)  ## orthogonal and include intercept 

betam002 = initiallap_mat(ind = datx$ind,y = datx$obs,xm = cbind(1, Bm[,-1]), lam = 0.001) 

## intercept is the covariate part, which is not identifiable
plot(betam002[group00==1,1],unique(datx$meanx)[group00==1] )
abline(0,1)
plot(betam002[group00==2,1],unique(datx$meanx)[group00==2] )
abline(0,1)

betam002[group00==1,1] - unique(datx$meanx)[group00==1]
betam002[group00==2,1] - unique(datx$meanx)[group00==2]


##### based on betam002[,-1] !!!!! works! #####

x1 = x + 0.05*matrix(rnorm(nrow(x)*2),nrow(x),2)

repeat{
  groupb = kmeans(betam002[,-1],centers = 10, iter.max = 20)$cluster
  if(min(table(groupb))>1){break}
}

betam021 = refitFDAX(ind = datx$ind,tm = datx$time,x = x,y = datx$obs,group0 = groupb, knots = knots)

res021 = Spgrrx(indexy = datx0$ind,y = datx0$obs,x = cbind(x,Bm),weights = wts,
                betam0 = betam021,lam = 1.1)

plot(res021$beta)
plot(res021$group)
randIndex(group00, as.numeric(res021$group))


### use median? ####
betam021median = apply(betam002[,-1], 1, median)
groupb2 = as.numeric(cut(betam021median,quantile(betam021median,seq(0,1,by=0.1)), include.lowest = TRUE))
res022 = refitFDAX(ind = datx$ind,tm = datx$time,x = x,y = datx$obs,group0 = groupb2, knots = knots)
betam022 = res022$alpha[groupb2,]

res022= FDAXsubgroup(ind = datx$ind,x = x, tm = datx$time,y = datx$obs,P = 2,
                   betam0 = betam022,group0 = groupb2,knots = knots, K0 = 10,max.step = 5,
                   lam = 1.1,maxiter = 50,tolabs = 1e-4,tolrel = 1e-2)


plot(res022$betam)
plot(res022$groupest)
randIndex(group00, as.numeric(res022$group), correct = FALSE)


betam023 = betam022 + 1*matrix(rnorm(nrow(betam022)*ncol(betam022)),nrow(betam022),ncol(betam022))

res023= FDAXsubgroup(ind = datx$ind,x = x, tm = datx$time,y = datx$obs,P = 2,
                     betam0 = betam023,group0 = groupb2,knots = knots, K0 = 10,max.step = 10,
                     lam = 1,maxiter = 50,tolabs = 1e-4,tolrel = 1e-2)

plot(res023$betam)
plot(res023$group)
randIndex(group00, as.numeric(res023$group), correct = FALSE)


########### BIC ##########
source("BICvalue.R")
source("simdatx.R")
source("initial.R")
source("refitFDA.R")
source("FDAsubgroup.R")
library(flexclust)

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


t1 = Sys.time()
restest = subfunx(1, sig200 = sig200, lam00 = lamj00, mvec00 = mvec00,
                  ncl00 = ncl00, lamvec = lamvec, 
                  funlist = funlist, eigenlist = eigenlist,xlist= xlist)
t2 = Sys.time()

xlist2 = list(nx = 2, meanx = 0, sdx = 1, etag = matrix(c(-0.5,-0.5,0.5,0.5),nrow = 2))
restest2 = subfunx(10, sig200 = sig200, lam00 = lamj00, mvec00 = mvec00,
                  ncl00 = ncl00, lamvec = lamvec, 
                  funlist = funlist, eigenlist = eigenlist,xlist= xlist2)
