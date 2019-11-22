###### test with common coefficients ######


source("simdatv2.R")
source("initialv2.R")
source("FDAsubgroup2.R")

funlist21 = list(Vectorize(function(x){4*(x-0.5)^2 + 1}),
                 Vectorize(function(x){2.5*exp(-25*(x-0.25)^2) + 2*exp(-50*(x-0.75)^2)}))


eigenlist21 = list(Vectorize(function(x){sqrt(2)*sin(pi*x)}),
                   Vectorize(function(x){ sqrt(2)*cos(pi*x)}))

grids = seq(0,1, by = 0.001)

datv2 = simdatv2(xlist = list(nx = 2, meanx =0, sdx = 1, eta = c(1,2)),
                 sig2 = 0.04,lamj = c(0.1,0.05),mvec = c(10,30),ncl = 50,
                 funlist = funlist21, eigenlist = eigenlist21, 
                 grids = grids, seed = 20 + 4452)

ind = datv2$ind
tm = datv2$time
x = as.matrix(datv2[,c("x1","x2")])
y = datv2$obs
knots = seq(0,1,length.out = 5)[2:4]

xm = as.matrix(datv2[,c("x1","x2")])
inilist = initialv2(ind = datv2$ind,x = xm,tm = datv2$time,
                    y = datv2$obs,knots = knots,lam = 0.001 )

betam0 = inilist$betam
eta0 = inilist$eta


res1= FDAsubgroupv2(ind = datv2$ind,x = xm, tm = datv2$time,y = datv2$obs,P = 2,
                     betam0 = betam0, eta0 = eta0, knots = knots, initial = "EMgroup", K0 = 10,
                     lam = 0.24,maxiter = 50,tolabs = 1e-4,tolrel = 1e-2)


res2 = FDAsubgroupv2(ind = datv2$ind,x = xm, tm = datv2$time,y = datv2$obs,P = 2,
                   betam0 = betam0, eta0 = eta0,knots = knots, initial = "EMcov", K0 = 10,
                   lam = 0.24,maxiter = 50,tolabs = 1e-4,tolrel = 1e-2)
