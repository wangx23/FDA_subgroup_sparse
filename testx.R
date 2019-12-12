
source("simdatx.R")
source("initial.R")
source("refitFDA.R")
source("FDAsubgroup.R")


# funlist21 = list(Vectorize(function(x){4*(x-0.5)^2 + 1}),
#                Vectorize(function(x){2.5*exp(-25*(x-0.25)^2) + 2*exp(-50*(x-0.75)^2)}))

funlist21 = list(Vectorize(function(x){4*(x/50)^(0.2)}),
               Vectorize(function(x){-4*(x/50)^(0.2)}))

eigenlist21 = list(Vectorize(function(x){sqrt(2)*sin(pi*x)}),
                 Vectorize(function(x){ sqrt(2)*cos(pi*x)}))

xlist = list(nx = 2, meanx = 0, sdx = 1, etag = matrix(c(-2,-2,2,2),nrow = 2))
sig200 = 0.01
lamj00 = c(0.1,0.05)

datx = simdatx(xlist = xlist,
                 sig2 = sig200,lamj = lamj00,mvec = c(20,30),ncl = 50,
                 funlist = funlist21, eigenlist = eigenlist21, seed = 20 + 4452)

ind = datx$ind
tm = datx$time
x = as.matrix(datx[,c("x1","x2")])
y = datx$obs
knots = seq(0,1,length.out = 5)[2:4]

group00 = unique(datx[,c("group","ind")])[,1]


res00 = refitFDAX(ind = datx$ind,tm = datx$time,x = x, y = datx$obs,P = 2,group0 = group00, 
                    knots = knots)
betam00 = res00$alpha[group00,]


betam01 = initialcoefx(ind = datx$ind,x = x,tm = datx$time,
                    y = datx$obs,knots = knots,lam = 0.01)

betam02 = betam00 + rnorm(nrow(betam00)*ncol(betam00))*0.5


res1= FDAXsubgroup(ind = datx$ind,x = x, tm = datx$time,y = datx$obs,P = 2,
                    betam0 = betam02,knots = knots, K0 = 10,
                    lam = 0.25,maxiter = 50,tolabs = 1e-4,tolrel = 1e-2)

betam0 = initialcoef2(ind = datx$ind,tm = datx$time,y = datx$obs - datx$meanx,
                      knots = knots,lam = 0.01)

res2 = FDAsubgroup(ind = datx$ind, tm = datx$time,y = datx$obs - datx$meanx, P= 2,
                   betam0 = betam01[,-c(1:2)], knots = knots, K0 = 10,
                   lam = 0.24,maxiter = 50,tolabs = 1e-4,tolrel = 1e-2)
table(res2$groupest)



########### no covariance ####
library(Spgr)
library(orthogonalsplinebasis)

xlist = list(nx = 2, meanx = 0, sdx = 1, etag = matrix(c(-2,-2,2,2),nrow = 2))

datx0 = simdatx(xlist = xlist,
               sig2 =0.01,lamj = c(0,0),mvec = c(50,50),ncl = 50,
               funlist = funlist21, eigenlist = eigenlist21, seed = 20 + 4452)


ntotal = length(datx0$obs)
knots = seq(0,1,length.out = 5)[2:4]
knotsall = c(rep(0,4),knots, rep(1,4))
obasisobj = OBasis(knotsall)
Bm = evaluate(obasisobj,datx0$time)[,-1]  ## orthogonal and include intercept 



x = as.matrix(datx0[,c("x1","x2")])

wts = rep(1, 100*(100-1)/2)

betam001 = cal_initialrx(indexy = datx0$ind,y = datx0$obs - datx0$meanx,x = Bm)
betam0011 = cal_initialrx(indexy = datx0$ind,y = datx0$obs,x = Bm[,-1]) ###
betam0012 = cal_initialrx(indexy = datx0$ind,y = datx0$obs,x = cbind(1, Bm[,-1])) ###

mean12 = rowSums(Bm[,-1] * betam0012[datx0$ind,-1])
residuals12 = datx0$obs - mean12

plot(betam0012[group00==1,1],unique(datx0$meanx)[group00==1] )
abline(0,1)
plot(betam0012[group00==2,1],unique(datx0$meanx)[group00==2] )
abline(0,1)

betam0012[group00==1,1] - unique(datx0$meanx)[group00==1]
betam0012[group00==2,1] - unique(datx0$meanx)[group00==2]

bx0 = cal_initialrx(indexy = datx0$ind,y = residuals12,x = x)

bx01 = cal_initialrx(indexy = datx0$ind,y = residuals12,x = cbind(1,x))



### estimate of x

betam0011[,1]


betam0013 = cal_initialrx(indexy = datx0$ind,y = datx0$obs - betam0011[datx0$ind,1],x = Bm) ###

res03 = Spgrrx(indexy = datx0$ind,y = datx0$obs,x = cbind(x,Bm),weights = wts,
               betam0 = cbind(xlist$etag[group00,] ,betam001),lam = 0.8)
plot(res03$group)


betam0012 = cal_initialrx(indexy = datx0$ind,y = datx0$obs - betam0011[datx0$ind,1],
                          x = Bm)



lm(datx0$obs[1:50] ~ Bm[1:50,-1])

plot(rowSums(unique(x * xlist$etag[datx0$group,])), betam0011[,1])

bx0 = cal_initialx(y = betam0011[,1],x = unique(x))


res01 = Spgrrx(indexy = datx0$ind,y = datx0$obs - datx0$meanx,x = Bm[,-1],weights = wts,
               betam0 = betam001,lam = 0.2)
plot(res01$group)
plot(res01$beta)

betam002 = cal_initialx(y = datx0$obs - datx0$mean,x = x)
res02 = Spgrx(y = datx0$obs - datx0$mean,x = x,weights = wts,
               betam0 = betam002,lam = 0.2)

group00 = unique(datx0[,c("group","ind")])[,1]
indsum = as.numeric(table(datx0$ind))

res03 = Spgrrx(indexy = datx0$ind,y = datx0$obs,x = cbind(x,Bm),weights = wts,
                betam0 = cbind(xlist$etag[group00,] ,betam0012),lam = 0.4)
plot(res03$group)


betam0011 = cal_initialrx(indexy = datx0$ind,y = datx0$obs,x = cbind(1,Bm))

res031 = Spgrrx(indexy = datx0$ind,y = datx0$obs,x = cbind(x,Bm),weights = wts,
               betam0 = cbind(bx0,betam0012),lam = 0.)
plot(res031$group)





betam000 = cbind(xlist$etag[group00,] + rnorm(100*2)*0.5,betam001)
betam004 = initialcoefx(ind = datx0$ind,tm = datx0$time,x = x,y = datx0$obs,knots = knots,lam = 0.001)

betam005 = initialcoef(ind = datx0$ind,tm = datx0$time,y = datx0$obs - lm(datx0$obs~x)$fitted.values,knots = knots,lam = 0.001)

res04 = Spgrrx(indexy = datx0$ind,y = datx0$obs,x = cbind(x,Bm),weights = wts,
               betam0 = betam004,lam = 0.5)


plot(res04$group)
plot(res04$beta)

grk = kmeans(betam005,8)$cluster
betam0041 = refitINDX(ind = datx0$ind,tm = datx0$time,x = x,y = datx0$obs, knots = knots,group0 = grk)


res05 = Spgrrx(indexy = datx0$ind,y = datx0$obs,x = cbind(x,Bm),weights = wts,
               betam0 = betam0041,lam = 0.6)
plot(res05$group)

betam006 = initialcoef(ind = datx0$ind,tm = datx0$time,y = datx0$obs,knots = knots,lam = 0.01)

res06 = Spgrrx(indexy = datx0$ind,y = datx0$obs,x = Bm,weights = wts,
               betam0 = betam006,lam = 0.4)

refit06 = refitINDX(ind = datx0$ind, y = datx0$obs, tm = datx0$time,x = x,group0 = res06$group, knots =knots) 


res07 = Spgrrx(indexy = datx0$ind,y = datx0$obs,x = cbind(x,Bm),weights = wts,
               betam0 = refit06,lam = 0.6)



#### overall model 
residuals = datx0$obs - lm(datx0$obs ~ 0 + cbind(x, Bm))$fitted

### median for each ind

plot(group00, as.numeric(by(residuals,datx0$ind,median)))

boxplot(as.numeric(by(residuals,datx0$ind,median)) ~ group00)

residm = as.numeric(by(residuals,datx0$ind,median))
### initial groups based on residuals ###
group1 = as.numeric(as.factor(cut(residm, quantile(residm, prob = 0:10 / 10, names = FALSE), include = TRUE)))

betam008 = refitINDX(ind = datx0$ind,tm = datx0$time,x = x,y = datx0$obs,group0 = group1, knots = knots)

res08 = Spgrrx(indexy = datx0$ind,y = datx0$obs,x = cbind(x,Bm),weights = wts,
               betam0 = betam008,lam = 1)



#### random ? 

coef00 = coef(lm(datx0$obs ~ 0 + cbind(x, Bm)))
betam009 = matrix(coef00,100,9, byrow = TRUE) + matrix(rnorm(100*9),100,9)*0.5
res09 = Spgrrx(indexy = datx0$ind,y = datx0$obs,x = cbind(x,Bm),weights = wts,
               betam0 = betam009,lam = 0.4)



##### each group has its own spline part coefficients but the x part assuming the same?
