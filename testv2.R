###### test with common coefficients ######


source("simdatv2.R")
source("initialv2.R")

funlist21 = list(Vectorize(function(x){4*(x-0.5)^2 + 1}),
                 Vectorize(function(x){2.5*exp(-25*(x-0.25)^2) + 2*exp(-50*(x-0.75)^2)}))


eigenlist21 = list(Vectorize(function(x){sqrt(2)*sin(pi*x)}),
                   Vectorize(function(x){ sqrt(2)*cos(pi*x)}))

grids = seq(0,1, by = 0.001)

datv2 = simdatv2(xlist = list(nx = 2, meanx =0, sdx = 1, eta = c(1,2)),
                 sig2 = 0.1,lamj = c(0.2,0.1),mvec = c(10,20),ncl = 50,
                 funlist = funlist21, eigenlist = eigenlist21, 
                 grids = grids, seed = mm + 4452)

ind = datv2$ind
tm = datv2$time
x = as.matrix(datv2[,c("x1","x2")])
y = datv2$obs
knots = seq(0,1,length.out = 4)[2:3]

xm = as.matrix(datv2[,c("x1","x2")])
inilist = initialv2(ind = datv2$ind,x = xm,tm = datv2$time,
                    y = datv2$obs,knots = knots,lam = 0.0001 )

betam0 = inilist$betam
eta0 = inilist$eta
