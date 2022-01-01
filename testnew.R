#### test for examples with x involved
source("simdatx.R")
source("initial.R")
source("refitFDA.R")
source("FDAsubgroup.R")
library(flexclust)
library(Spgr)
library(cluster)


# funlist21 = list(Vectorize(function(x){4*(x-0.5)^2 + 1}),
#                Vectorize(function(x){2.5*exp(-25*(x-0.25)^2) + 2*exp(-50*(x-0.75)^2)}))

funlist21 = list(Vectorize(function(x){cos(2*pi*x)} ),
                Vectorize(function(x){1-2*exp(-6*x)} )
)


fun1 = function(x){funlist21[[1]](x)}
fun2 = function(x){funlist21[[2]](x)}

curve(fun1 ,0,1, ylim = c(-3,3))
curve(fun2 ,0,1, add = TRUE)

eigenlist21 = list(Vectorize(function(x){sqrt(2)*sin(pi*x)}),
                   Vectorize(function(x){ sqrt(2)*cos(pi*x)}))

xlist = list(nx = 2, meanx = 0, sdx = 1, etag = matrix(c(-2,-2,2,2),nrow = 2))
sig200 = 0.01
lamj00 = c(0.1,0.05)


#### without covariance 

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
Bm1 <- bs(datx0$time, knots = knots, intercept = T, Boundary.knots = c(0,1))

x = as.matrix(datx0[,c("x1","x2")]) ## 

wts = rep(1, 100*(100-1)/2) ## weights 


betam001 = initiallap_mat(ind = datx0$ind,y = datx0$obs,xm = Bm, lam = 0.001)
betam002 = initiallap_mat(ind = datx0$ind,y = datx0$obs,xm = cbind(1, Bm[,-1]), lam = 0.001) ### the mean term is an intercept 
betam003 = initial_indiv(ind = datx0$ind, y = datx0$obs, xm = Bm,lam = 0.001) ## each has its own coefficients, almost the same as betam001

betam0031 = initial_indiv(ind = datx0$ind, y = datx0$obs, xm = Bm1,lam = 0)

betam0032 = initial_indiv(ind = datx0$ind, y = datx0$obs, xm = cbind(1,Bm[,-1]),lam = 0)

meanx0 <- unique(datx0$meanx) + c(1,-1)[group00]

meanx01 <- unique(datx0$meanx)

plot(betam001, betam003)
plot(betam0031[,1], meanx0)
plot(betam0032[,1], meanx0)

x1 <- betam001 %*% t(Bm0)


meanx0est <- x1

plot(x1, meanx0)
plot(x1, meanx01)



xm1 <- cbind(1,unique(x))


betax01 <- cal_initialx(y = x1,x = xm1)
betax01 <- cal_initialrx2(indexy = 1:100,y = x1,x = xm1, K0 = 10)

solve(t(xm1[1:50,])%*%xm1[1:50,]) %*% t(xm1[1:50,]) %*% x1[1:50]

solve(t(xm1[51:100,])%*%xm1[51:100,]) %*% t(xm1[51:100,]) %*% x1[51:100]


resx01 <- Spgrx(y = x1, x = xm1,weights = wts, betam0 = betax01,lam = 0.9)




#### remove mean part 

y1 <- datx0$obs - x1[datx0$ind]   ## starting from 0 now

betam004 = initiallap_mat(ind = datx0$ind, y = y1, xm = Bm,lam = 0.001)

betam004 = initial_indiv(ind = datx0$ind, y = y1, xm = Bm,lam = 0.001)

groupb <- kmeans(betam004, 10, nstart = 100)$cluster

betam03 = refitINDX(ind = datx0$ind,tm = datx0$time,x = x,y = datx0$obs,group0 = groupb, knots = knots)

plot(betam03)

########## 
knots = seq(0,1,length.out = 9)[2:8]

t1 <- runif(200)
y1 = 1 + funlist21[[2]](t1) + 0.1*rnorm(200)
bm1 <- bs(t1,knots = knots, intercept = T, Boundary.knots = c(0,1))

solve(t(bm1)%*% bm1)%*% t(bm1) %*%y1




####### with covariance matrix #####


xlist = list(nx = 2, meanx = 0, sdx = 1, etag = matrix(c(-2,-2,2,2),nrow = 2))
sig200 = 0.01
lamj00 = c(0.1,0.05)

datx = simdatx(xlist = xlist,
               sig2 = sig200,lamj = lamj00,mvec = c(20,30),ncl = 50,
               funlist = funlist21, eigenlist = eigenlist21, seed = 20 + 4452)

group00 = unique(datx[,c("group","ind")])[,1] ## true group 

library(ggplot2)
ggplot(data = datx, aes(x = time, y = obs, group = ind, color = as.factor(group)))+geom_line() + geom_point() + theme_bw()

ntotal = length(datx$obs)
knots = seq(0,1,length.out = 6)[2:5]
knotsall = c(rep(0,4),knots, rep(1,4))
obasisobj = OBasis(knotsall)
Bm = evaluate(obasisobj,datx$time)  ## orthogonal and include intercept 
Bm0 = evaluate(obasisobj,0) 

x = as.matrix(datx[,c("x1","x2")]) ## 
ind = datx$ind
y = datx$obs
tm = datx$time
xm = Bm

wts = rep(1, 100*(100-1)/2) ## weights 


betam003 = initiallap_mat(ind = datx$ind, y = datx$obs, xm = Bm,lam = 0.001) ## each has its own coefficients, almost the same as betam001

meanx0 <- unique(datx$meanx) + c(1,-1)[group00]

meanx01 <- unique(datx$meanx)

x1 <- betam003 %*% t(Bm0)

plot(x1, meanx01)
abline(0,1)

plot(x1, meanx0)
abline(0,1)


xm1 <- cbind(1,unique(x))

solve(t(xm1[1:50,])%*%xm1[1:50,]) %*% t(xm1[1:50,]) %*% x1[1:50]

solve(t(xm1[51:100,])%*%xm1[51:100,]) %*% t(xm1[51:100,]) %*% x1[51:100]


betax01 <- cal_initialx(y = x1,x = xm1)
betax01 <- cal_initialrx2(indexy = 1:100,y = x1,x = xm1, K0 = 10)

resx01 <- Spgrx(y = x1, x = xm1,weights = wts, betam0 = betax01,lam = 0.01)

plot(resx01$beta[,1])

plot(resx01$beta[,2])

plot(betax01[,3])


y1 <- datx$obs - x1[datx$ind]   ## starting from 0 now

betam004 = initiallap_mat(ind = datx$ind, y = y1, xm = Bm,lam = 0.001)

#betam004 = initial_indiv(ind = datx$ind, y = y1, xm = Bm,lam = 0.001)

groupb <- kmeans(betam004, 10, nstart = 100)$cluster
table(groupb)

groupb1 = pam(x = betam004, k = 10,diss = FALSE)$clustering
table(groupb1)
groupb1[groupb1==10] <-  9

betam03 = refitINDX(ind = datx$ind,tm = datx$time,x = x,y = datx$obs,group0 = groupb1, knots = knots)

plot(betam03[,1])


resx01 <- Spgrrx(indexy = datx$ind, y = datx$obs - x1[datx$ind], x = Bm,weights = wts, betam0 = betam004,lam = 0.5)

groupb = resx01$group
table(groupb)

t1 <- Sys.time()
res1= FDAXsubgroup(ind = datx$ind,x = x, tm = datx$time,y = datx$obs,P = 2,
                   betam0 = betam03,knots = knots, K0 = 10, group0 = groupb1,
                   lam = 1,maxiter = 50,tolabs = 1e-4,tolrel = 1e-2, max.step = 5)

t2 <- Sys.time()

plot(res1$betam)

randIndex(res1$groupest, group00)
library(scclust)
my_dist <- distances(betam004)
my_clustering <- sc_clustering(my_dist, 6)
groupb2 <- as.vector(my_clustering) + 1

betam03 = refitINDX(ind = datx$ind,tm = datx$time,x = x,y = datx$obs,group0 = groupb2, knots = knots)


t1 <- Sys.time()
res2= FDAXsubgroup(ind = datx$ind,x = x, tm = datx$time,y = datx$obs,P = 2,
                   betam0 = betam03,knots = knots, K0 = 10, group0 = groupb2,
                   lam = 1.3,maxiter = 50,tolabs = 1e-4,tolrel = 1e-2, 
                   max.step = 5)

t2 <- Sys.time()

plot(res2$betam)
table(res2$groupest)

randIndex(res2$groupest, group00)
