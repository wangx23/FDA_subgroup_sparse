###### parallel test, time is randomly selected from grids ##### 
library(flexclust)
library(doParallel)
library(mcclust)

source("simdat.R")
source("initial.R")
source("EMgroup.R")
source("FDAsubgroup.R")
source("BICvalue.R")
source("fclust.R")

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
funlist2 = list(Vectorize(function(x){2*(x/50)^(0.2)}),
                Vectorize(function(x){-2*(x/50)^(0.2)})
)


eigenlist = list(Vectorize(function(x){sqrt(2)*sin(2*pi*x)}),
                 Vectorize(function(x){ sqrt(2)*cos(2*pi*x)}))


grids = seq(0,1, by = 0.001)

mm = 1

dat = simdat2(sig2 = 0.1,lamj = c(0.2,0.1),mvec = c(5,20),ncl = 50,
              funlist = funlist2, grids = grids, seed = mm + 4452)
group0 = rep(1:length(funlist2),each = 50)



##### proposed algorithm ###
knots1 = seq(0,1,length.out = 6)[2:5]
betam0 = initialcoef(ind = dat$ind,tm = dat$time,y = dat$obs,knots = knots1,
                     lamv = seq(0,20,by = 0.5)[-1])

lamvec = seq(0.2,1,by = 0.05)
BICm = matrix(0,length(lamvec), 3)

t1 = Sys.time()
for(Pv in 1:3)
{
  for(j in 1:length(lamvec))
  {
    resi = try(FDAsubgroup(ind = dat$ind,tm = dat$time,y = dat$obs,P = Pv,
                           betam0 = betam0, knots = knots1,
                           lam = lamvec[j],maxiter = 50,tolabs = 1e-4,tolrel = 1e-2))
    errori = inherits(resi,"try-error")
    if(errori)
    {
      BICm[j, Pv] = 99999
    }else{
      BICm[j, Pv] = BICvalue(resi)
    }
  }
}
t2 = Sys.time()


inds = which(BICm == min(BICm), arr.ind = TRUE)

res = FDAsubgroup(ind = dat$ind,tm = dat$time,y = dat$obs,P = inds[2],
                  betam0 = betam0, knots = seq(0,1,length.out = 6)[2:5],
                  lam = lamvec[inds[1]],maxiter = 50,tolabs = 1e-4,tolrel = 1e-2)

group_fda = res$groupest

ari_fda = randIndex(group_fda, group0)
vi_fda = vi.dist(group_fda, group0)
ise_fda = ISEFDAmean(obj = res,group0, grids = grids,funlist = funlist2)

ise_eig  ### ise of eigenfunction estimates 

####### without covariance structure #####
nknots = 4;
order = 4;
p =  order + nknots;
#time points 
TT = max(table(dat$ind))
timerange = seq(0, 1, length.out = TT);
n = 100

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

t3 = Sys.time()
lamvec2 = seq(0.5,2,by=0.1)

BIC2vec = rep(0, length(lamvec2))

for(j in 1:length(lamvec2))
{
  sol_finalj = prclust_admm(X, y=as.vector(dat$obs), diagD, B_ini0, index,
                            gamma1 = 0.005, gamma2 = lamvec2[j], 
                            theta=1, tau = 2, 
                            n, p,  max_iter= 200,
                            eps_abs=1e-4, eps_rel=1e-2)
  BIC2vec[j] = BIC2(obj = sol_finalj,dat$obs, dat$ind, X, basis, n, 0.005, D)
  
}
t4 = Sys.time()

inds2 = which.min(BIC2vec)

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
ise_ind = ISEINDmean(t(sol_final$B), group0, grids, funlist2, 
                     nknots =4, order = 4)



#### J&S method  #####

datlist = list(x = dat$obs,
                curve = dat$ind,
                timeindex = match(dat$time,grids))

fit.js = fitfclust(data=datlist,q=8,h=1,p=8,K=2,maxit=30,grid=grids,plot=F,trace=F)
group_js = fclust.pred(fit.js)$class
ari_js= randIndex(group_js, group0)
vi_js = vi.dist(group_js, group0)
ise_js = ISEJSmean(obj = fit.js,group0,grids,funlist = funlist2)
