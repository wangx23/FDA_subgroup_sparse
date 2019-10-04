###### parallel test, time is randomly selected from grids ##### 
#### three groups #####
#### compare the best results one algorithm can have 
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
library(Matrix);
library(matrixcalc);
library(igraph);
library(Rcpp);
library(RcppArmadillo);
library(Lclust)


#sourceCpp("../LongitudinalClustering-master/admmmcp_code.cpp")
#source("../LongitudinalClustering-master/clustering_functions.R")
#source("../LongitudinalClustering-master/BICfunctions.R")


#### mean functions ###
funlist3 = list(Vectorize(function(x){sqrt(2)*sin(4*pi*x)}),
                Vectorize(function(x){exp(-10*(x-0.25)^2)}),
                Vectorize(function(x){1.5*x - 1})
)

# curve(funlist3[[1]](x),0,1)
# curve(funlist3[[2]](x),0,1,add = TRUE, col="red")
# curve(funlist3[[3]](x),0,1,add = TRUE, col="blue")


eigenlist = list(Vectorize(function(x){sqrt(2)*sin(2*pi*x)}),
                 Vectorize(function(x){ sqrt(2)*cos(2*pi*x)}))

grids = seq(0,1, by = 0.001)



subfun = function(mm)
{
  dat = simdat2(sig2 = 0.01,lamj = c(0.1,0.06),mvec = c(5,20),ncl = 50,
                funlist = funlist3, grids = grids, seed = mm + 4452)
  group0 = rep(1:length(funlist3),each = 50)
  ng0 = length(unique(group0))
  
  ##### output matrix for indicators
  outputmat = matrix(0, 4, 3)
  rownames(outputmat) = c("ng","ari","vi","lam")
  colnames(outputmat) = c("fda","ind","js")
  
  isemeanmat = matrix(0, length(group0), 3)
  colnames(isemeanmat) = c("fda","ind","js")
  
  confmat = matrix(0, 3, ng0*2) ## results based on the confusion matrix
  colnames(confmat) = paste0(rep(c("true","true_rate"),each=2),rep(1:ng0,2))
  rownames(confmat) = c("fda","ind","js")
  
  
  
  
  ##### proposed algorithm ###
  knots1 = seq(0,1,length.out = 5)[2:4]
  betam0 = initialcoef(ind = dat$ind,tm = dat$time,y = dat$obs,knots = knots1,
                       lamv = seq(0,20,by = 0.5)[-1])
  
  lamvec = seq(0.15,0.4,length.out = 50)
  BICfda = rep(0,length(lamvec))
  ng_fda_vec = ari_fda_vec = vi_fda_vec = rep(-99, length(lamvec))
  
  ### fixed P = 2 
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
    }
    

    
    
  }
  
  ### maximum ari_fda_vec
  
  inds1 = which(ari_fda_vec == max(ari_fda_vec[ng_fda_vec==3]))[1]
  
  res = FDAsubgroup(ind = dat$ind,tm = dat$time,y = dat$obs,P = 2,
                   betam0 = betam0, knots = knots1,
                   lam = lamvec[inds1],maxiter = 50,tolabs = 1e-4,tolrel = 1e-2)
  
  ise_fda = ISEFDAmean(obj = res,group0, grids = grids,funlist = funlist3)
  
  outputmat[,1] = c(ng_fda_vec[inds1], ari_fda_vec[inds1], vi_fda_vec[inds1], lamvec[inds1])
  isemeanmat[,1] = ise_fda
  
  
  ise_eig = ISEFDAeig(obj = res,grids = grids,eigenlist = eigenlist) ### ise of eigenfunction estimates 
  mse_lamj = (res$lamj -  c(0.1,0.05))^2
  
  
  #### If the true number of groups is selected, then calculate the following ####
  # true_fda = true_rate_fda = rep(0, ng0)
  # if(ng_fda ==ng0 )
  # {
  #   confusion_fda = table(group_fda,group0)
  #   true_fda = diag(confusion_fda) ##  number of correct pairs
  #   true_rate_fda = true_fda/colSums(confusion_fda)
  # }
  # 
  # 
  # confmat[1,] = c(true_fda, true_rate_fda)
  
  
  ####### without covariance structure ######
  nknots = 3;
  order = 4;
  p =  order + nknots;
  #time points 
  TT = max(table(dat$ind))
  timerange = seq(0, 1, length.out = TT);
  n = 150
  
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
  
  lamvec2 = seq(0.5,1.5,by=0.025)
  
  BIC2vec = rep(0, length(lamvec2))
  ng_ind_vec = ari_ind_vec = vi_ind_vec = rep(0, length(lamvec2))
  
  for(j in 1:length(lamvec2))
  {
    sol_finalj = prclust_admm(X, y=as.vector(dat$obs), diagD, B_ini0, index,
                              gamma1 = 0.005, gamma2 = lamvec2[j], 
                              theta=1, tau = 2, 
                              n, p,  max_iter= 200,
                              eps_abs=1e-4, eps_rel=1e-2)
    BIC2vec[j] = BIC2(obj = sol_finalj,dat$obs, dat$ind, X, basis, n, 0.005, D)
    
    Ad_final = create_adjacency(sol_finalj$V, n);
    G_final = graph.adjacency(Ad_final, mode = 'upper')
    #clustering membership
    cls_final = components(G_final);
    #number of clusters
    ng_ind = cls_final$no
    group_ind = cls_final$me
    ari_ind= randIndex(group_ind, group0)
    vi_ind = vi.dist(group_ind, group0)
    
    ng_ind_vec[j] = ng_ind
    ari_ind_vec[j] = ari_ind
    vi_ind_vec[j] = vi_ind
    
  }

  
  inds2 = which(ari_ind_vec == max(ari_ind_vec[ng_ind_vec==3]))[1]
  
  sol_final = prclust_admm(X, y=as.vector(dat$obs), diagD, B_ini0, index,
                           gamma1 = 0.005, gamma2 = lamvec2[inds2], 
                           theta=1, tau = 2, n, p,  max_iter=200,
                           eps_abs=1e-4, eps_rel=1e-2)
  #result:
  Ad_final = create_adjacency(sol_final$V, n);
  G_final = graph.adjacency(Ad_final, mode = 'upper')
  #clustering membership
  cls_final = components(G_final);
  #number of clusters
  ng_ind = cls_final$no
  group_ind = cls_final$me
  
  ari_ind= randIndex(group_ind, group0)
  vi_ind = vi.dist(group_ind, group0)
  ise_ind = ISEINDmean(t(sol_final$B), group0, grids,timerange, funlist3, 
                       nknots = nknots, order = 4)
  
  outputmat[,2] = c(ng_ind, ari_ind, vi_ind, lamvec2[inds2],0)
  isemeanmat[,2] = ise_ind
  
  
  # true_ind = true_rate_ind = rep(0, ng0)
  # if(ng_ind ==ng0)
  # {
  #   confusion_ind = table(group_ind,group0)
  #   true_ind = diag(confusion_ind) ##  number of correct pairs
  #   true_rate_ind = true_ind/colSums(confusion_ind)
  # }
  # confmat[2,]= c(true_ind, true_rate_ind)
  
  
  #### J&S method  #####
  
  datlist = list(x = dat$obs,
                 curve = dat$ind,
                 timeindex = match(dat$time,grids))
  
  fit.js = fitfclust(data=datlist,q=7,h=2,p=8,K=3,maxit=30,grid=grids,plot=F,trace=F)
  
  group_js = fclust.pred(fit.js)$class.pred
  ng_js = length(unique(group_js))
  ari_js= randIndex(group_js, group0)
  vi_js = vi.dist(group_js, group0)
  ise_js = ISEJSmean(obj = fit.js,group0,group_js,grids,funlist = funlist3)
  outputmat[,3] = c(ng_js, ari_js, vi_js,0,0)
  isemeanmat[,3] = ise_js
  
  
  # true_js = true_rate_js = rep(0, ng0)
  # if(ng_js ==ng0)
  # {
  #   confusion_js = table(group_js,group0)
  #   true_js = diag(confusion_js) ##  number of correct pairs
  #   true_rate_js = true_js/colSums(confusion_js)
  # }
  # 
  # confmat[3,] = c(true_js, true_rate_js)
  
  output = list(outputmat = outputmat, isemeanmat = isemeanmat, 
                ise_eig = ise_eig, mse_lamj = mse_lamj)
  return(output)
}


# t1 = Sys.time()
# res1 = subfun(10)
# t2 = Sys.time()

cl <- makeCluster(24)  
registerDoParallel(cl)  
resultcomp3 <- foreach(mm=1:120,
                       .packages=c("flexclust","orthogonalsplinebasis","plyr","fda","Lclust","igraph","mcclust"),.errorhandling = "remove") %dopar%  subfun(mm)
stopCluster(cl) 
save(resultcomp3,file = "../result/resultcompgrids3.RData")

