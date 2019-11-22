######## test for FDA_subgroup simdat2 #######
###### test with common coefficients ######


source("simdat.R")
source("initial.R")
source("FDAsubgroup.R")

funlist21 = list(Vectorize(function(x){4*(x-0.5)^2 + 1}),
                 Vectorize(function(x){2.5*exp(-25*(x-0.25)^2) + 2*exp(-50*(x-0.75)^2)}))


eigenlist21 = list(Vectorize(function(x){sqrt(2)*sin(pi*x)}),
                   Vectorize(function(x){ sqrt(2)*cos(pi*x)}))

grids = seq(0,1, by = 0.001)

dat1 = simdat2(sig2 = 0.04,lamj = c(0.1,0.05),mvec = c(20,30),ncl = 50,
                 funlist = funlist21, eigenlist = eigenlist21, 
                 grids = grids, seed = 1+ 4452)

ind = dat1$ind
tm = dat1$time
y = dat1$obs

knots = seq(0,1,length.out = 5)[2:4]

betam0 = initialcoef(ind = dat1$ind,tm = dat1$time,y = dat1$obs,knots = knots,
                     lamv =seq(0,20,by = 0.5)[-1])

betam0 = initialcoef2(ind = dat1$ind,tm = dat1$time,
                    y = dat1$obs,knots = knots,lam = 0.001)


res1 = FDAsubgroup(ind = dat1$ind,tm = dat1$time,y = dat1$obs,P = 2,
                       betam0 = betam0, knots = knots, initial = "EMgroup",
                       lam = 0.25,maxiter = 50,tolabs = 1e-4,tolrel = 1e-2)

res12 = FDAsubgroup(ind = dat1$ind,tm = dat1$time,y = dat1$obs,P = 2,
                   betam0 = betam0, knots = knots, initial = "EMgroup2",
                   lam = 0.23,maxiter = 50,tolabs = 1e-4,tolrel = 1e-2)


res2 = FDAsubgroup(ind = dat1$ind,tm = dat1$time,y = dat1$obs,P = 2,
                   betam0 = betam0, knots = knots, initial = "EMcov", K0 = 10,
                   lam = 0.23,maxiter = 50,tolabs = 1e-4,tolrel = 1e-2)
