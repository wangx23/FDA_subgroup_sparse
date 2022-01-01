source("BICvalue.R")
source("simdatx.R")
source("initial.R")
source("refitFDA.R")
source("FDAsubgroup.R")
source("FDAXsubgroup.R")
source("initialgroup.R")
library(flexclust)
library(cluster)


funlist1 = list(Vectorize(function(x){cos(2*pi*x)}),
                Vectorize(function(x){1-2*exp(-6*x)}),
                Vectorize(function(x){-1.5*x})
)

eigenlist = list(Vectorize(function(x){sqrt(2)*sin(pi*x)}),
                 Vectorize(function(x){ sqrt(2)*cos(pi*x)}))

xlist1 = list(nx = 2, meanx = 0, sdx = 1, 
              etag = matrix(c(-0.5,-0.5,0.5,0.5,1.5,1.5),nrow = 2))
sig200 = 0.04
# lamj00 = c(0.2,0.1)
lamj00 = c(0.3, 0.2)
mvec00 = c(10,20)
#lamvec = seq(0.2,0.8,by = 0.025)

lamvec1 = seq(0.2,0.6, by = 0.01)
lamvec2 = seq(0.3,0.6,  by = 0.01)
lamvec3 = seq(0.555,0.6, by = 0.001)

lamvec1 = 0.2
lamvec2 = 0.3
#lamvec3 = 0.2
lamvec4 = 0.2

ncl00 = 50

datx = simdatx2(xlist = xlist1,
                sig2 = sig200,lamj = lamj00,mvec = mvec00,ncl = ncl00,
                funlist = funlist1, eigenlist = eigenlist, seed = 20+ 4452)




group0 = rep(1:length(funlist1),each = ncl00)
ng0 = length(unique(group0))


x = as.matrix(datx[,c("x1","x2")])

ntotal = length(datx$obs)
knots = seq(0,1,length.out = 5)[2:4]
knotsall = c(rep(0,4),knots, rep(1,4))
obasisobj = OBasis(knotsall)
Bm = evaluate(obasisobj,datx$time)  ## orthogonal and include intercept 
K0 = 15
ming = 3

betam002 = initiallap_mat(ind = datx$ind,y = datx$obs,
                          xm = cbind(x, Bm), lam = 0.001) 

groupb22 = initialgroup_ks(betam002,K0 = K0, ming = ming)

groupbb = groupb22

resi = FDAXsubgroup(ind = datx$ind,x = x, tm = datx$time,y = datx$obs,P = 2,
                    betam0 = betam002,group0 = groupbb,knots = knots, K0 = K0,
                    max.step = 5,
                    lam = 0.43,maxiter = 500,tolabs = 1e-4,tolrel = 1e-2)
plot(resi$betam)
resi$groupest
randIndex(resi$groupest,group0)
