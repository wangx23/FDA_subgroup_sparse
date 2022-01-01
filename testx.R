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

fun1 = function(x){funlist21[[1]](x)}
fun2 = function(x){funlist21[[2]](x)}

curve(fun1,0,1, ylim = c(-1.5,1.5))
curve(fun2,0,1, add = TRUE)

eigenlist21 = list(Vectorize(function(x){sqrt(2)*sin(pi*x)}),
                 Vectorize(function(x){ sqrt(2)*cos(pi*x)}))

xlist = list(nx = 2, meanx = 0, sdx = 1, etag = matrix(c(-2,-2,2,2),nrow = 2))
sig200 = 0.01
lamj00 = c(0.1,0.05)

datx = simdatx(xlist = xlist,
                 sig2 = sig200,lamj = lamj00,mvec = c(100,100),ncl = 50,
                 funlist = funlist21, eigenlist = eigenlist21, seed = 20 + 4452)

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
                   betam0 = betatd01,knots = knots, K0 = 10,
                   lam = 0.25,maxiter = 50,tolabs = 1e-4,tolrel = 1e-2)


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


betam003 = initial_indiv(ind = datx$ind, y = datx$obs, xm = cbind(1, Bm[,-1]), lam = 0) ## each has its own coefficients


### note: traditional method does not work

########### no covariance#######
### cal_initialrx
library(Spgr)
library(orthogonalsplinebasis)

xlist = list(nx = 2, meanx = 0, sdx = 1, etag = matrix(c(-2,-2,2,2),nrow = 2)) # two groups of x

datx0 = simdatx(xlist = xlist,
               sig2 =0.01,lamj = c(0,0),mvec = c(200,200),ncl = 50,
               funlist = funlist21, eigenlist = eigenlist21, seed = 20 + 4452)
group00 = unique(datx0[,c("group","ind")])[,1] ## true group 

ntotal = length(datx0$obs)
knots = seq(0,1,length.out = 6)[2:5]
knotsall = c(rep(0,4),knots, rep(1,4))
obasisobj = OBasis(knotsall)
Bm = evaluate(obasisobj,datx0$time)  ## orthogonal and include intercept 
Bm0 = evaluate(obasisobj,0) 

x = as.matrix(datx0[,c("x1","x2")]) ## 

wts = rep(1, 100*(100-1)/2) ## weights 

library(splines)
Bm1 <- bs(datx0$time, knots = knots, intercept = T, Boundary.knots = c(0,1))
Bm10 <- bs(0, knots = knots, intercept = T, Boundary.knots = c(0,1))

### Bm includes the intercept ###
## if parametric part is known 
betam001 = initiallap_mat(ind = datx0$ind,y = datx0$obs,xm = Bm, lam = 0.001)
betam002 = initiallap_mat(ind = datx0$ind,y = datx0$obs,xm = cbind(1, Bm[,-1]), lam = 0.001) ### the mean term is an intercept 
betam003 = initial_indiv(ind = datx0$ind, y = datx0$obs, xm = Bm,lam = 0) ## each has its own coefficients

betam0031 = initial_indiv(ind = datx0$ind, y = datx0$obs, xm = Bm1,lam = 0)

x0 <- betam003 %*% t(Bm0)

mean002 = rowSums(Bm[,-1] * betam002[datx0$ind,-1])
residuals002 = datx0$obs - mean002

## intercept is the covariate part, which is not identifiable
plot(betam002[group00==1,1],unique(datx0$meanx)[group00==1] )
abline(0,1)
plot(betam002[group00==2,1],unique(datx0$meanx)[group00==2] )
abline(0,1)

plot(betam003[group00==1,1],unique(datx0$meanx)[group00==1] )
abline(0,1)
plot(betam003[group00==2,1],unique(datx0$meanx)[group00==2] )
abline(0,1)


betam002[group00==1,1] - unique(datx0$meanx)[group00==1]
betam002[group00==2,1] - unique(datx0$meanx)[group00==2]

bx0 = initiallap_mat(ind = datx0$ind,y = residuals002,xm = x)
bx01 = initiallap_mat(ind = datx0$ind,y = residuals002,xm = cbind(1,x))





### traditional method ####
betam01 = initiallapx(ind = datx0$ind,tm = datx0$time,x = x,y = datx0$obs,knots = knots,lam = 0.001)

res01 = Spgrrx(indexy = datx0$ind,y = datx0$obs,x = cbind(x,Bm),weights = wts,
               betam0 = betam01,lam = 0.65)

plot(res01$beta)
plot(res01$group)
randIndex(group00, as.numeric(res01$group))



## remove x part?
betam02 = initiallap_mat(ind = datx0$ind,y = datx0$obs - betam002[datx0$ind,1],xm = Bm) ###

# consider betam002[datx0$ind,1] as an estimate of the intercept
res02 = Spgrrx(indexy = datx0$ind,y = datx0$obs,x = cbind(x,Bm),weights = wts,
               betam0 = cbind(xlist$etag[group00,] ,betam002),lam = 0.8) # x part uses the true value as initial to see whether betam01 is a good initial of beta
plot(res02$group)



##### based on betam002[,-1] !!!!! works! #####

x1 = x + 0.05*matrix(rnorm(nrow(x)*2),nrow(x),2)

repeat{
  groupb = kmeans(betam003[,-1],centers = 2, iter.max = 20)$cluster
  if(min(table(groupb))>1){break}
}

betam021 = refitINDX(ind = datx0$ind,tm = datx0$time,x = x,y = datx0$obs,group0 = groupb, knots = knots)

res021 = Spgrrx(indexy = datx0$ind,y = datx0$obs,x = cbind(x,Bm),weights = wts,
               betam0 = betam021,lam = 1.1)

plot(res021$beta)
plot(res021$group)
randIndex(group00, as.numeric(res021$group))


### use median? ####
betam021median = apply(betam002[,-1], 1, median)
groupb2 = as.numeric(cut(betam021median,quantile(betam021median,seq(0,1,by=0.1)), include.lowest = TRUE))
betam022 = refitINDX(ind = datx0$ind,tm = datx0$time,x = x,y = datx0$obs,group0 = groupb2, knots = knots)

res022 = Spgrrx(indexy = datx0$ind,y = datx0$obs,x = cbind(x,Bm),weights = wts,
                betam0 = betam022,lam = 1.1)
plot(res022$beta)
plot(res022$group)
randIndex(group00, as.numeric(res022$group), correct = FALSE)


betam023 = betam022 + 1*matrix(rnorm(nrow(betam022)*ncol(betam022)),nrow(betam022),ncol(betam022))

res023 = Spgrrx(indexy = datx0$ind,y = datx0$obs,x = cbind(x,Bm),weights = wts,
                betam0 = betam023,lam = 1.7)
plot(res023$beta)
plot(res023$group)
randIndex(group00, as.numeric(res023$group), correct = FALSE)


#### median for each ind based on overall model####
residuals = datx0$obs - lm(datx0$obs ~ 0 + cbind(x, Bm))$fitted
plot(group00, as.numeric(by(residuals,datx0$ind,median)))
boxplot(as.numeric(by(residuals,datx0$ind,median)) ~ group00)

residm = as.numeric(by(residuals,datx0$ind,median))


### initial groups based on residuals ###
group3 = as.numeric(as.factor(cut(residm, quantile(residm, prob = 0:10 / 10, names = FALSE), include = TRUE)))
betam03 = refitINDX(ind = datx0$ind,tm = datx0$time,x = x,y = datx0$obs,group0 = group3, knots = knots)

res03 = Spgrrx(indexy = datx0$ind,y = datx0$obs,x = cbind(x,Bm),weights = wts,
               betam0 = betam03,lam = 1)

plot(res03$group)

#### random #####

coef00 = coef(lm(datx0$obs ~ 0 + cbind(x, Bm))) ## overall model coefficients
betam04 = matrix(coef00,100,9, byrow = TRUE) + matrix(rnorm(100*9),100,9)*0.5
res04 = Spgrrx(indexy = datx0$ind,y = datx0$obs,x = cbind(x,Bm),weights = wts,
               betam0 = betam04,lam = 0.4)
plot(res04$group)


##### 

#est51 = est_initial(ind = datx0$ind,tm = datx0$time,x = x,y = datx0$obs,knots = knots)
est52 = initialv2(ind = datx0$ind,x = x,tm = datx0$time,y = datx0$obs,knots = knots)


plot(datx0$obs[datx0$group==1] - est52$bmx[datx0$group==1], datx0$meanx[datx0$group==1])
abline(0,1)

  
plot(est52$bmx[datx0$group==2], datx0$meanx[datx0$group==2])
abline(0,1)

eta51 = initial_indiv(ind = datx0$ind,y = datx0$obs - est52$bmx,xm = x, lam = 0)


#### add random to x 

x1 = x + 0.05*matrix(rnorm(nrow(x)*2),nrow(x),2)

betam06 = initiallapx(ind = datx0$ind,tm = datx0$time,x = x1,y = datx0$obs,knots = knots,lam = 0.001)


#### linear regression ####
group1 = datx0$group==1
group2 = datx0$group==2

lm(datx0$obs[group1]~ 0+ cbind(x1[group1,], Bm[group1,]))

lm(datx0$obs[group2]~ 0+ cbind(x1[group2,], Bm[group2,]))


lm(datx0$obs[datx0$ind==4]~ 0+ cbind(x1[datx0$ind==4,], Bm[datx0$ind==4,]))


############ a new distance based on observations #########


