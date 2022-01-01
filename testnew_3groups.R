##### three groups ####
source("BICvalue.R")
source("simdatx.R")
source("initial.R")
source("refitFDA.R")
source("FDAsubgroup.R")
source("FDAXsubgroup.R")
library(flexclust)



#### mean functions set 1 ####
# funlist1 = list(Vectorize(function(x){cos(2*pi*x)}),
#                 Vectorize(function(x){1-2*exp(-6*x)}),
#                 Vectorize(function(x){-1.5*x})
# )


funlist1 = list(Vectorize(function(x){-2 + sin(4*pi*x)}),
                Vectorize(function(x){2 + exp(-10*(x)^2) -1}),
                Vectorize(function(x){0.5 + 1.5*x})
)

curve(funlist1[[1]](x),0,1,ylim = c(-2,2))
curve(funlist1[[2]](x),0,1, add = TRUE, col = "red")
curve(funlist1[[3]](x),0,1, add = TRUE,col = "blue")

eigenlist = list(Vectorize(function(x){sqrt(2)*sin(pi*x)}),
                 Vectorize(function(x){ sqrt(2)*cos(pi*x)}))

xlist1 = list(nx = 2, meanx = 0, sdx = 1, 
              etag = matrix(c(-2,-2,2,2,0.5,0.5),nrow = 2))
sig200 = 0.04
lamj00 = c(0.2,0.1)
mvec00 = c(50,50)
ncl00 = 50
lamvec = seq(0.2,2,by = 0.1)

datx = simdatx(xlist = xlist1,
               sig2 = sig200,lamj = lamj00,mvec = mvec00,ncl = ncl00,
               funlist = funlist1, eigenlist = eigenlist, seed = 1 + 4452)
group00 = unique(datx[,c("group","ind")])[,1] ## true group 

ggplot(data = datx, aes(x = time, y = obs, group = ind, color = as.factor(group)))+geom_line()  + theme_bw()


ntotal = length(datx$obs)
knots = seq(0,1,length.out = 6)[2:5]
knotsall = c(rep(0,4),knots, rep(1,4))
obasisobj = OBasis(knotsall)
Bm = evaluate(obasisobj,datx$time)  ## orthogonal and include intercept 
Bm0 = evaluate(obasisobj,0) 


x = as.matrix(datx[,c("x1","x2")]) ## 
betam003 = initiallap_mat(ind = datx$ind, y = datx$obs, xm = Bm,lam = 0.001) #

meanx0 <- unique(datx$meanx) + c(-2,2,0.5)[group00]
meanx01 <- unique(datx$meanx)

x1 <- betam003 %*% t(Bm0)

plot(x1, meanx01)
abline(0,1)

plot(x1, meanx0)
abline(0,1)


xm1 <- cbind(1,unique(x))

solve(t(xm1[group00==1,])%*%xm1[group00==1,]) %*% t(xm1[group00==1,]) %*% x1[group00==1]

solve(t(xm1[group00==2,])%*%xm1[group00==2,]) %*% t(xm1[group00==2,]) %*% x1[group00==2]

solve(t(xm1[group00==3,])%*%xm1[group00==3,]) %*% t(xm1[group00==3,]) %*% x1[group00==3]


y1 <- datx$obs - x1[datx$ind]   ## starting from 0 now

betam004 = initiallap_mat(ind = datx$ind, y = y1, xm = Bm,lam = 0.001)


my_dist <- distances(betam004)
my_clustering <- sc_clustering(my_dist, 5)
groupb2 <- as.vector(my_clustering) + 1
table(groupb2)

betam03 = refitINDX(ind = datx$ind,tm = datx$time,x = x,y = datx$obs,group0 = groupb2, knots = knots)

plot(betam03)


res2= FDAXsubgroup(ind = datx$ind,x = x, tm = datx$time,y = datx$obs,P = 2,
                   betam0 = betam03,knots = knots, K0 = 10, group0 = groupb2,
                   lam = 0.5,maxiter = 100,tolabs = 1e-4,tolrel = 1e-2, 
                   max.step = 5)


plot(res2$betam)

table(res2$groupest)

randIndex(res2$groupest, group00, correct = F)
randIndex(res2$groupest, group00)
res2$sig2
res2$lamj
res2$niteration


res2r <- refitFDAX(ind = datx$ind, x = x, tm = datx$time, y = datx$obs, P = 2,
                   group0 = res2$groupest,knots = knots, maxiter = 100)
