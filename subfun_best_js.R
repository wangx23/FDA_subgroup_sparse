###### best soultion three groups with other values can be changed ######
###### parallel test, time is randomly selected from grids ##### 
#### three groups #####
#### compare the best results one algorithm can have 
#### only compare with JS 
library(flexclust)
library(doParallel)
library(mcclust)

source("simdat.R")
source("initial.R")
source("EMgroup.R")
source("FDAsubgroup.R")
source("BICvalue.R")
source("fclust2.R")
source("IMSE.R")

library(plyr);
library(fda);#bspline basis





#### mean functions ###
funlist31 = list(Vectorize(function(x){sqrt(2)*sin(4*pi*x)}),
                Vectorize(function(x){exp(-10*(x-0.25)^2)}),
                Vectorize(function(x){1.5*x - 1})
)


eigenlist = list(Vectorize(function(x){sqrt(2)*sin(2*pi*x)}),
                 Vectorize(function(x){ sqrt(2)*cos(2*pi*x)}))

grids = seq(0,1, by = 0.001)

lamvec = seq(0.23,0.32,length.out = 10)

lam00 = c(0.15,0.1)
sig200 = 0.05
mvec00 = c(30, 40)
ncl00 = 50
# lamvec2 = exp(seq(-0.5,0.3,length.out = 50))

subfun_best_js = function(mm, sig200, lam00 = c(0.1, 0.05), mvec00 = c(5,20), ncl00= 50,
                       lamvec, funlist, eigenlist, Kjs, hjs)
{
  dat = simdat2(sig2 = sig200,lamj = lam00,mvec = mvec00,ncl = ncl00,
                funlist = funlist, eigenlist = eigenlist, 
                grids = grids, seed = mm + 4452)
  group0 = rep(1:length(funlist),each = ncl00)
  ng0 = length(unique(group0))
  
  ##### output matrix for indicators
  outputmat = matrix(0, 4, 2)
  rownames(outputmat) = c("ng","ari","vi","lam") ### ng = Kjs
  colnames(outputmat) = c("fda","js")
  
  
  isemeanmat = matrix(0, length(group0), 2)
  colnames(isemeanmat) = c("fda","js")
  
  
  
  #### J&S method  #####
  set.seed(mm + 1058)
  datlist = list(x = dat$obs,
                 curve = dat$ind,
                 timeindex = match(dat$time,grids))
  
  fit_js = fitfclust(data=datlist,q=7,h=1,p=8,K=2,maxit=30,grid=grids,plot=F,trace=F)
  
  group_js = fclust.pred(fit_js)$class.pred
  ng_js = length(unique(group_js))
  ari_js= randIndex(group_js, group0)
  #ari_js
  vi_js = vi.dist(group_js, group0)
  ise_js = ISEJSmean(obj = fit_js,group0,group_js,grids,funlist = funlist)
  outputmat[,2] = c(ng_js, ari_js, vi_js,0)
  isemeanmat[,2] = ise_js
  
  
  
  ##### proposed algorithm ###
  knots1 = seq(0,1,length.out = 5)[2:4]
  betam0 = initialcoef(ind = dat$ind,tm = dat$time,y = dat$obs,knots = knots1,
                       lamv =seq(0,20,by = 0.5)[-1])
  #betam0 = initialcoef2(ind = dat$ind,tm = dat$time,y = dat$obs,knots = knots1)

  BICfda = rep(0,length(lamvec))
  ng_fda_vec = ari_fda_vec = vi_fda_vec = rep(-99, length(lamvec))
  groupmat = matrix(0,length(group0),length(lamvec))
  

  ### fixed P = 2
  #betamj = betam0
  for(j in 1:length(lamvec))
  {
    resj = try(FDAsubgroup(ind = dat$ind,tm = dat$time,y = dat$obs,P = 2,
                           betam0 = betam0, knots = knots1,
                           lam = lamvec[j],maxiter = 50,tolabs = 1e-4,tolrel = 1e-2))
    errorj = inherits(resj,"try-error")
    if(errorj)
    {
      BICfda[j] = 99999
    }else{
      BICfda[j] = BICvalue(resj)
      group_fdaj = resj$groupest
      ng_fdaj = length(unique(group_fdaj))
      ari_fdaj = randIndex(group_fdaj, group0)
      vi_fdaj = vi.dist(group_fdaj, group0)

      ng_fda_vec[j] = ng_fdaj
      ari_fda_vec[j] = ari_fdaj
      vi_fda_vec[j] = vi_fdaj
      groupmat[,j] = resj$groupest
      betaarray[,,j] = betamj
    }
  }

  ### maximum ari_fda_vec with certain groups

  inds1 = which(ari_fda_vec == max(ari_fda_vec))[1]
  # 
  res = FDAsubgroup(ind = dat$ind,tm = dat$time,y = dat$obs,P = 2,
                    betam0 = betam0, knots = knots1,
                    lam = lamvec[inds1],maxiter = 50,tolabs = 1e-4,tolrel = 1e-2,K0 = 10)
  #res$groupest
  # # 
  #randIndex(res$groupest, group0)

  ise_fda = ISEFDAmean(obj = res,group0, grids = grids,funlist = funlist)

  outputmat[,1] = c(ng_fda_vec[inds1], ari_fda_vec[inds1], vi_fda_vec[inds1], lamvec[inds1])
  isemeanmat[,1] = ise_fda
  # 
  # 
  ise_eig = ISEFDAeig(obj = res,grids = grids,eigenlist = eigenlist) ### ise of eigenfunction estimates
  mse_lamj = (res$lamj -  lam00)^2
  # # 
  # # 
  # # 
  # 
  output = list(outputmat = outputmat, isemeanmat = isemeanmat,
                ise_eig = ise_eig, mse_lamj = mse_lamj)
  
  # output = list(outputmat = outputmat, isemeanmat = isemeanmat)
  return(output)
}


res1 = subfun_best_js(4,sig200, lam00, mvec, ncl00,lamvec, funlist31, eigenlist, Kjs = 3, hjs = 2)
res1$outputmat

#### another set of mean functions #####

funlist21 = list(Vectorize(function(x){4*(x-0.5)^2 + 1}),
                Vectorize(function(x){2.5*exp(-25*(x-0.25)^2) + 2*exp(-50*(x-0.75)^2)}))


eigenlist21 = list(Vectorize(function(x){sqrt(2)*sin(pi*x)}),
                 Vectorize(function(x){ sqrt(2)*cos(pi*x)}))

res2 = subfun_best_js(4,sig200 = 0.02, lam00 = c(0.15,0.5), mvec00 = c(20,50),ncl00 = 50,lamvec, 
                      funlist = funlist21, eigenlist = eigenlist21, Kjs = 2, hjs = 1)
res1$outputmat

testari = matrix(0,50,2)
testng = matrix(0,50,2)
for(mm  in 1:50)
{
  res2 = subfun_best_js(mm,sig200 = 0.05, lam00 = c(0.15,0.1), mvec00 = c(30,40),ncl00 = 50,lamvec, 
                        funlist = funlist21, eigenlist = eigenlist21, Kjs = 2, hjs = 1)
  testng[mm,] = res2$outputmat[1,]
  testari[mm,]= res2$outputmat[2,]
  print(mm)
}

