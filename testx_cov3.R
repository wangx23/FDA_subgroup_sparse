###### testx with covariance, three groups #####
source("BICvalue.R")
source("simdatx.R")
source("initial.R")
source("refitFDA.R")
source("FDAsubgroup.R")
library(flexclust)



#### mean functions set 1 ####
funlist1 = list(Vectorize(function(x){cos(2*pi*x)}),
                 Vectorize(function(x){1-2*exp(-6*x)}),
                 Vectorize(function(x){-1.5*x})
)

curve(funlist1[[1]](x),0,1,ylim = c(-2,2))
curve(funlist1[[2]](x),0,1, add = TRUE, col = "red")
curve(funlist1[[3]](x),0,1, add = TRUE,col = "blue")

eigenlist = list(Vectorize(function(x){sqrt(2)*sin(pi*x)}),
                   Vectorize(function(x){ sqrt(2)*cos(pi*x)}))

xlist1 = list(nx = 2, meanx = 0, sdx = 1, etag = matrix(c(-0.5,-0.5,0.5,0.5,1,1),nrow = 2))
sig200 = 0.04
lamj00 = c(0.2,0.1)
mvec00 = c(10,20)
ncl00 = 100
lamvec = seq(0.2,2,by = 0.1)

datx = simdatx(xlist = xlist1,
               sig2 = sig200,lamj = lamj00,mvec = mvec00,ncl = ncl00,
               funlist = funlist1, eigenlist = eigenlist, seed = 1 + 4452)



ind = datx$ind
tm = datx$time
x = as.matrix(datx[,c("x1","x2")])
y = datx$obs
knots = seq(0,1,length.out = 5)[2:4]
group00 = unique(datx[,c("group","ind")])[,1] ## true group 



#### true group ###
res00 = refitFDAX(ind = datx$ind,tm = datx$time,x = x, y = datx$obs,P = 2,group0 = group00, 
                  knots = knots)
betatd00 = res00$alpha[group00,] ### true eta and beta 



betam02 = betatd00 + 0.5*rnorm(nrow(betatd00)*ncol(betatd00))
res2= FDAXsubgroup(ind = datx$ind,x = x, tm = datx$time,y = datx$obs,P = 3,
                   betam0 = betam02,knots = knots, K0 = 10, max.step = 10,
                   lam = 0.5,maxiter = 50,tolabs = 1e-4,tolrel = 1e-2)
randIndex(res2$groupest, group00)




###### initial values ###

library(Spgr)
library(orthogonalsplinebasis)

ntotal = length(datx$obs)
knots = seq(0,1,length.out = 5)[2:4]
knotsall = c(rep(0,4),knots, rep(1,4))
obasisobj = OBasis(knotsall)
Bm = evaluate(obasisobj,datx$time)  ## orthogonal and include intercept 

betam002 = initiallap_mat(ind = datx$ind,y = datx$obs,xm = cbind(1, Bm[,-1]), lam = 0.001) 
plot(betam002[group00==1,1],unique(datx$meanx)[group00==1] )
abline(0,1)
plot(betam002[group00==2,1],unique(datx$meanx)[group00==2] )
abline(0,1)

betam002[group00==1,1] - unique(datx$meanx)[group00==1]
betam002[group00==2,1] - unique(datx$meanx)[group00==2]

library(cluster)

repeat{
  groupb = kmeans(betam002[,-1],centers = 15, iter.max = 20)$cluster
  if(min(table(groupb))>1){break}
}


groupb1 = pam(x = betam002[,-1], k = 15,diss = FALSE)$clustering

betam021 = refitFDAX(ind = datx$ind,tm = datx$time,x = x,y = datx$obs,group0 = groupb, knots = knots)$alpha[groupb,]

betam0211 = refitFDAX(ind = datx$ind,tm = datx$time,x = x,y = datx$obs,group0 = groupb1, knots = knots)$alpha[groupb1,]

betam0212 = refitINDX(ind = datx$ind, tm = datx$time, x = x,y = datx$obs, group0 = groupb, knots = knots)


res021 = FDAXsubgroup(ind = datx$ind,x = x, tm = datx$time,y = datx$obs,P = 2,
                      betam0 = betam021,
                      group0 = groupb,knots = knots, K0 = 10,max.step = 5,
                      lam = 0.4,maxiter = 50,tolabs = 1e-4,tolrel = 1e-2)

randIndex(res021$groupest, group00)
res0211 = FDAXsubgroup(ind = datx$ind,x = x, tm = datx$time,y = datx$obs,P = 2,
                       betam0 = betam0211,
                       group0 = groupb1,knots = knots, K0 = 10,max.step = 5,
                       lam = 0.4,maxiter = 50,tolabs = 1e-4,tolrel = 1e-2)

randIndex(res0211$groupest, group00)

res022 = FDAXsubgroup(ind = datx$ind,x = x, tm = datx$time,y = datx$obs,P = 2,
                       betam0 = betam0212,
                       group0 = groupb,knots = knots, K0 = 10,max.step = 5,
                       lam = 0.35,maxiter = 50,tolabs = 1e-4,tolrel = 1e-2)

randIndex(res022$groupest, group00)

#### use median ###
betam021median = apply(betam002[,-1], 1, median)
groupb2 = as.numeric(cut(betam021median,quantile(betam021median,seq(0,1,by=0.1)), include.lowest = TRUE))
res022 = refitFDAX(ind = datx$ind,tm = datx$time,x = x,y = datx$obs,group0 = groupb2, knots = knots)
betam022 = res022$alpha[groupb2,]

res022= FDAXsubgroup(ind = datx$ind,x = x, tm = datx$time,y = datx$obs,P = 2,
                     betam0 = betam022,group0 = groupb2,knots = knots, K0 = 10,max.step = 5,
                     lam = 0.3,maxiter = 50,tolabs = 1e-4,tolrel = 1e-2)
randIndex(res022$groupest, group00)


#### simulation report likevalue, estimated number of groups and lambda 
funlist1 = list(Vectorize(function(x){cos(2*pi*x)}),
                Vectorize(function(x){1-2*exp(-6*x)}),
                Vectorize(function(x){-1.5*x})
)

eigenlist = list(Vectorize(function(x){sqrt(2)*sin(pi*x)}),
                 Vectorize(function(x){ sqrt(2)*cos(pi*x)}))

xlist1 = list(nx = 2, meanx = 0, sdx = 1, etag = matrix(c(-0.5,-0.5,0.5,0.5,1,1),nrow = 2))
sig200 = 0.04
lamj00 = c(0.2,0.1)
mvec00 = c(10,20)
ncl00 = 100
lamvec = seq(0.05,1,by = 0.05)




subfunx_all = function(mm, sig200, lam00, mvec00, ncl00,
                   lamvec, funlist, eigenlist, xlist, K0 = 10)
{
  
  datx = simdatx(xlist = xlist,
                 sig2 = sig200,lamj = lamj00,mvec = mvec00,ncl = ncl00,
                 funlist = funlist, eigenlist = eigenlist, seed = mm+ 4452)
  
  group0 = rep(1:length(funlist),each = ncl00)
  ng0 = length(unique(group0))
  
  
  x = as.matrix(datx[,c("x1","x2")])
  
  ntotal = length(datx$obs)
  knots = seq(0,1,length.out = 5)[2:4]
  knotsall = c(rep(0,4),knots, rep(1,4))
  obasisobj = OBasis(knotsall)
  Bm = evaluate(obasisobj,datx$time)  ## orthogonal and include intercept 
  
  betam002 = initiallap_mat(ind = datx$ind,y = datx$obs,xm = cbind(1, Bm[,-1]), lam = 0.001) 
  
 # if(method == "median")
  #{
    betam021median = apply(betam002[,-1], 1, median)
    groupb21 = as.numeric(cut(betam021median,quantile(betam021median,seq(0,1,by=0.1)), include.lowest = TRUE))
  #}
  
 # if(method == "kmeans")
  #{
    groupb22 = initialgroup_ks(betam002[,-1],K0 = K0)

  #}
  
  #if(method == "pam")
  #{
    groupb23 = initialgroup_pam(betam002[,-1],K0 = K0, diss = FALSE)
  #}
  
  #if(method == "ydist")
  #{
    ny = length(unique(datx$ind))
    distmaty = matrix(0, ny, ny)
    xm = unique(as.matrix(datx[,c("x1","x2")]))
    t1 = Sys.time()
    for(i in 1:(ny-1))
    {
      tsi = sort(datx$time[datx$ind==i],index.return = TRUE)
      tmi = tsi$x
      yi = datx$obs[datx$ind==i][tsi$ix]
      
      sfuni = stepfun(tmi, c(yi,yi[length(yi)]),f=1, right = TRUE)
      for(j in (i+1):ny)
      {
        tsj = sort(datx$time[datx$ind==j],index.return = TRUE)
        tmj = tsj$x
        yj = datx$obs[datx$ind==j][tsj$ix]
        
        sfunj = stepfun(tmj, c(yj,yj[length(yj)]),f=1, right = TRUE)
        
        funij = function(x)
        {
          (sfuni(x) - sfunj(x))^2
        }
        funij = Vectorize(funij)
        
        value1 = stats::integrate(funij,0,1, subdivisions = 500)$value
        valueij = sqrt(value1) + sqrt(sum((xm[i,] - xm[j,])^2))
        distmaty[i,j] = valueij
      }
    }
    t2 = Sys.time()
    
    distmaty = distmaty + t(distmaty)
    groupb24 = pam(distmaty, K0, diss = TRUE)$clustering
 # }
  

  subestfun = function(groupbb)
  {
    
    nllmat = matrix(0, length(lamvec),3)
    arim = matrix(0, length(lamvec),3)
    ngest = matrix(0, length(lamvec),3)
    lamest1 = lamest2 = lamest3 = matrix(0, length(lamvec),3)
    
    #betam022 = refitINDX(ind = datx$ind, tm = datx$time, x=x, y = datx$obs, group0 = groupbb, knots = knots)
    
    betam022 = refitFDAX(ind = datx$ind,tm = datx$time,x = x,y = datx$obs,group0 = groupbb, knots = knots)$alpha[groupbb,]
    
    for(Pv in 1:3)
    {
      for(j in 1:length(lamvec))
      {
        resi = try(FDAXsubgroup(ind = datx$ind,x = x, tm = datx$time,y = datx$obs,P = Pv,
                                betam0 = betam022,group0 = groupbb,knots = knots, K0 = 10,
                                max.step = 5,
                                lam = lamvec[j],maxiter = 50,tolabs = 1e-4,tolrel = 1e-2))
        errori = inherits(resi,"try-error")
        if(errori)
        {
          nllmat[j,Pv] = 9999
        }else{
          nllmat[j,Pv] = resi$likevalue
          arim[j, Pv] = randIndex(resi$groupest,group0)
          lamest1[j,Pv]= resi$lamj[1]
          lamest2[j,Pv] = resi$lamj[2]
          lamest3[j,Pv] = resi$lamj[3]
          ngest[j,Pv] = length(unique(resi$groupest))
        }
      }
    }
    
    
    outlist = list(nllmat = nllmat, arim = arim, lamest1= lamest1, lamest2 = lamest2,
                   lamest3 = lamest3, ngest = ngest)
    return(outlist)
  }
  

  res1 = subestfun(groupbb = groupb21)
  res2 = subestfun(groupbb = groupb22)
  res3 = subestfun(groupbb = groupb23)
  res4 = subestfun(groupbb = groupb24)


  output = list(res_median = res1, res_kmeans = res2, res_pam = res3,
                res_ydist = res4)
  return(output)
}


#### set 2 ####
funlist2 = list(Vectorize(function(x){4*(x - 0.5)^2+1}),
                Vectorize(function(x){2.5*exp(-25*(x-0.25)^2) + 2*exp(-50*(x-0.75)^2)})
)

curve(funlist2[[1]](x),0,1,ylim = c(0,3))
curve(funlist2[[2]](x),0,1, add = TRUE, col = "red")



#### set 3 ####

funlist3 = list(Vectorize(function(x){sqrt(2)*sin(4*pi*x)}),
                Vectorize(function(x){exp(-10*(x-0.25)^2)})
)
curve(funlist3[[1]](x),0,1,ylim = c(-2,2))
curve(funlist3[[2]](x),0,1, add = TRUE, col = "red")
