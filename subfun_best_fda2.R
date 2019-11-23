########### some tests #####
library(flexclust)
library(doParallel)
library(mcclust)

source("simdatv2.R")
source("initialv2.R")
source("FDAsubgroup2.R")


library(plyr);
library(fda);#bspline basis


grids = seq(0,1, by = 0.001)
lamvec = seq(0.2,0.30, by= 0.01)


#### Kjs is a vector #####
subfun_best_fda2 = function(mm, eta00, sig200, lam00, mvec00, ncl00,
                          lamvec, funlist, eigenlist, Kjs= 2)
{
  datv2 = simdatv2(xlist = list(nx = 2, meanx =0, sdx = 1, eta = eta00),
                   sig2 = sig200,lamj = lam00,mvec = mvec00,ncl = ncl00,
                   funlist = funlist, eigenlist = eigenlist, 
                   grids = grids, seed = 20 + 4452)
  
  
  group00 = rep(1:length(funlist),each = ncl00)
  ng0 = length(unique(group00))
  
  ##### output matrix for indicators
  # outputmat = matrix(0, 4, 2)
  # rownames(outputmat) = c("ng","ari","vi","lam") ### ng = Kjs
  # colnames(outputmat) = c("fda","js")
  # 
  # 
  # isemeanmat = matrix(0, length(group0), 2)
  # colnames(isemeanmat) = c("fda","js")
  # 
  
  
  #### J&S method  #####
  # set.seed(mm + 1058)
  # datlist = list(x = dat$obs,
  #                curve = dat$ind,
  #                timeindex = match(dat$time,grids))
  
  
  # nk = length(Kjs)
  # ar_jsvec = rep(0,nk)
  # for(j in 1:nk)
  # {
  #   fit_js = fitfclust(data=datlist,q=7,h=Kjs[j]-1,p=8,K=Kjs[j],maxit=30,grid=grids,plot=F,trace=F)
  #   group_js = fclust.pred(fit_js)$class.pred
  #   ar_jsvec[j] = randIndex(group_js, group0)
  # }
  # 
  # indexjs = which.max(ar_jsvec)
  
  #fit_js = fitfclust(data=datlist,q=7,h=Kjs[indexjs]-1,p=8,K=Kjs[indexjs],maxit=30,grid=grids,plot=F,trace=F)
  
  # fit_js = fitfclust(data=datlist,q=7,h=Kjs-1,p=8,K=Kjs,maxit=30,grid=grids,plot=F,trace=F)
  # group_js = fclust.pred(fit_js)$class.pred
  # ng_js = length(unique(group_js))
  # ari_js= randIndex(group_js, group0)
  # vi_js = vi.dist(group_js, group0)
  # ise_js = ISEJSmean(obj = fit_js,group0,group_js,grids,funlist = funlist)
  # outputmat[,2] = c(ng_js, ari_js, vi_js,0)
  # isemeanmat[,2] = ise_js
  # 
  
  ##### proposed algorithm ###
  knots1 = seq(0,1,length.out = 5)[2:4]
  xm = as.matrix(datv2[,c("x1","x2")])
  inilist = initialv2(ind = datv2$ind,x = xm,tm = datv2$time,
                      y = datv2$obs,knots = knots1,lam = 0.001 )
  
  betam0 = inilist$betam
  eta0 = inilist$eta
  
  
 # BICfda = rep(0,length(lamvec))
  ng_fda_vec = ari_fda_vec = rep(-99, length(lamvec))
  vi_fda_vec = rep(99, length(lamvec))
  groupmat = matrix(0,length(group00),length(lamvec))
  
  ng_fda_vec2 = ari_fda_vec2 = rep(-99, length(lamvec))
  vi_fda_vec2 = rep(99, length(lamvec))
  
  groupmat2 = matrix(0,length(group00),length(lamvec))
  
  # ### fixed P = 2
  #betamj = betam0
  for(j in 1:length(lamvec))
  {
    resj1 = try(FDAsubgroupv2(ind = datv2$ind,x = xm, tm = datv2$time,y = datv2$obs,P = 2,
                             betam0 = betam0, eta0 = eta0, knots = knots1, initial = "EMgroup", K0 = 10,
                             lam = lamvec[j],maxiter = 50,tolabs = 1e-4,tolrel = 1e-2))
    
    errorj = inherits(resj1,"try-error")
    if(!errorj)
    {
      group_fdaj1 = resj1$groupest
      ng_fdaj1 = length(unique(group_fdaj1))
      ari_fdaj1 = randIndex(group_fdaj1, group00)
      vi_fdaj1 = vi.dist(group_fdaj1, group00)
      
      ng_fda_vec[j] = ng_fdaj1
      ari_fda_vec[j] = ari_fdaj1
      vi_fda_vec[j] = vi_fdaj1
      groupmat[,j] = resj1$groupest
    }
    
    resj2 = try(FDAsubgroupv2(ind = datv2$ind,x = xm, tm = datv2$time,y = datv2$obs,P = 2,
                              betam0 = betam0, eta0 = eta0, knots = knots1, initial = "EMcov", K0 = 10,
                              lam = lamvec[j],maxiter = 50,tolabs = 1e-4,tolrel = 1e-2))
    
    errorj2 = inherits(resj2,"try-error")
    if(!errorj2)
    {
      group_fdaj2 = resj2$groupest
      ng_fdaj2 = length(unique(group_fdaj2))
      ari_fdaj2 = randIndex(group_fdaj2, group00)
      vi_fdaj2 = vi.dist(group_fdaj2, group00)
      
      ng_fda_vec2[j] = ng_fdaj2
      ari_fda_vec2[j] = ari_fdaj2
      vi_fda_vec2[j] = vi_fdaj2
      groupmat2[,j] = resj2$groupest
    }
    
  }
  
  inds11 = which(ari_fda_vec == max(ari_fda_vec))[1]
  inds12 = which(vi_fda_vec == min(vi_fda_vec))[1]
  ng11 = ng_fda_vec[inds11]
  ng12 = ng_fda_vec[inds12]
  
  
  inds21 = which(ari_fda_vec2 == max(ari_fda_vec2))[1]
  inds22 = which(vi_fda_vec2 == min(vi_fda_vec2))[1]
  
  ng21 = ng_fda_vec2[inds21]
  ng22 = ng_fda_vec2[inds22]
  
  
  output = cbind(c(ari_fda_vec[inds11], vi_fda_vec[inds12]),c(ng11, ng12),
        c(ari_fda_vec2[inds21], vi_fda_vec2[inds22]),c(ng21, ng22))

  return(output)
}


funlist21 = list(Vectorize(function(x){4*(x-0.5)^2 + 1}),
                 Vectorize(function(x){2.5*exp(-25*(x-0.25)^2) + 2*exp(-50*(x-0.75)^2)}))


eigenlist21 = list(Vectorize(function(x){sqrt(2)*sin(pi*x)}),
                   Vectorize(function(x){ sqrt(2)*cos(pi*x)}))
# eta00 = c(1,2)
# lam00 = c(0.15,0.1)
# sig200 = 0.04
# mvec00 = c(20, 30)
# ncl00 = 50
# 
# subfun_best_fda2(1, eta00 = c(1,2),sig200 = 0.04,lam00 = c(0.1,0.05),mvec00 = c(10,30),ncl00 = 50,lamvec = lamvec,
#                  funlist = funlist21,eigenlist = eigenlist21)


#### parallel ####
library(doParallel)
cl <- makeCluster(18)  
registerDoParallel(cl)  
outputfda2_test = foreach(mm=1:50,
                       .packages=c("flexclust","orthogonalsplinebasis","plyr","fda","Lclust","igraph","mcclust"),.errorhandling = "remove") %dopar%  subfun_best_fda2(mm, eta00 = c(1,2),sig200 = 0.04,lam00 = c(0.1,0.05),mvec00 = c(10,30),ncl00 = 50,lamvec = lamvec,funlist = funlist21,eigenlist = eigenlist21)
stopCluster(cl) 
save(outputfda2_test, file = "../result/outputfda2_test.RData")