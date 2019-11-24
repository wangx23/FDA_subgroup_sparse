###### a function use Peng, J. and Paul, D. (2009) newton method and ADMM algorithm ###

## ind is the subject ID
## tm is the time
## y is the observation
## knots is the interior knots
## boundary, boundary knots
## lam is the tuning parameter in ADMM
## nu and gam: two fixed values in ADMM
## maxiter is the maximum number of iterations in the algorithm 
## tolabs and tolrel are two tolerance values in ADMM
## maxiterem, tolem, number of iterations and tolerance in EM (calculating initial values), K0 is the number of groups in Kmeans (calculating initial values)
## lamv is used to calculate initial values of beta

## sl.v: shrinkage steps for Newton, e.g., sl.v<-c(0.5,0.5,0.5) means the first three steps are 0.5
## max.step: max number of iterations of Newton

##tolnt: tolerance level to determine convergence in Newton in terms of the l_2 norm of the gradient 
##condtol: tolerance to determine singularity in Newton in terms of condition number of the Hessian matrix


source("initial.R")
source("fpcafunctions.R")
source("EMgroup.R")
source("revisedadmm.R")
source("getgroup.R")


# boundary = c(0,1)
# nu = 1
# gam = 3
# maxiter = 500
# tolabs = 1e-4
# tolrel = 1e-2
# maxiterem = 50
# tolem = 1e-3
# K0 = 10
# lamv = seq(0,20,by = 0.5)[-1]
# sl.v=rep(0.5,10)
# max.step = 20
# tolnt=1e-3
# condtol=1e+10
# seed = 2118
# P =2

#"EMgroup" or "EMcov"
FDAsubgroup = function(ind, tm, y, P = 2, betam0, knots, boundary = c(0,1),
                       lam = 0.5, nu = 1, gam = 3, maxiter = 500, initial ="EMcov",
                       tolabs = 1e-5, tolrel = 1e-3,
                       maxiterem = 50, tolem = 1e-3, K0 = 10, 
                       sl.v=rep(0.5,10), max.step = 20, 
                       tolnt=1e-3,condtol=1e+10,
                       seed = 2118)
{
  
  ### some values of the data set ###
  ntotal = length(y)
  knotsall = c(rep(boundary[1],4),knots, rep(boundary[2],4))
  obasisobj = OBasis(knotsall)
  Bm = evaluate(obasisobj,tm)  ## orthogonal
  
  p = ncol(Bm)
  uind = unique(ind)
  n = length(uind)
  npair = n*(n-1)/2
  

  Bmlist = list()
  ylist = list()
  for(i in 1:n)
  {
    indi = ind == uind[i]
    Bmlist[[i]] = Bm[indi,]
    ylist[[i]] = y[indi]
  }
  
  if(length(sl.v)<max.step){
    sl.v<-c(sl.v,rep(1,max.step-length(sl.v)))
  }else{
    sl.v<-sl.v[1:max.step]
  }
  
  
  set.seed(seed)
  
  if(initial == "EMgroup")
  {
    ### initial value of theta, lamj and sig2
    repeat{
      group0 = kmeans(betam0,K0)$cluster
      if(min(table(group0))>1){break}
    }
    
    #group0 = group
    res0 = EMgroup(ind = ind,tm = tm, y = y, knots= knots, 
                   group0 = group0, P = P, betam0 = betam0, method ="lm",boundary = boundary,
                   maxiter = maxiterem, tol = tolem)
    
    alpm = res0$alpm
    
    sig2 = res0$sig2
    theta = res0$theta
    lamj = res0$lamj
    betam = res0$alpm[group0,]
  }
  
  if(initial == "EMgroup2")  ###
  {
    ### initial value of theta, lamj and sig2

    group0 = kmeans(betam0,K0)$cluster
    #group0 = group
    res0 = EMgroup(ind = ind,tm = tm, y = y, knots= knots, 
                   group0 = group0, P = P, betam0 = betam0, method = "mean",boundary = boundary,maxiter = maxiterem, tol = tolem)
    
    alpm = res0$alpm
    
    sig2 = res0$sig2
    theta = res0$theta
    lamj = res0$lamj
    betam = res0$alpm[group0,]
  }
  
  
  
  if(initial == "EMcov")
  {
    
    group0 = kmeans(betam0,K0)$cluster
    alpm0 = do.call("rbind",by(betam0, group0, colMeans,simplify = TRUE))
    
    yresid = rep(0, ntotal)
    for(i in 1:n)
    {
      indi = ind == uind[i]
      yresid[indi] = ylist[[i]] - Bmlist[[i]]%*%alpm0[group0[i],]
    }
    res0 = EMcov(ind = ind,tm = tm,yresid = yresid, group0 = group0, knots = knots,P = P,
                 boundary = boundary)
    sig2 = res0$sig2
    theta = res0$theta
    lamj = res0$lamj
    betam = alpm0[group0,]
  }
  
  
  
  #### initial and matrix in ADMM ####
  
  ### difference matrix
  Dmat = matrix(0,n*(n-1)/2,n)
  for(j in 1:(n-1))
  {
    indexj = (n-1 + n-j+1)*(j-1)/2
    indexvj = indexj + (1:(n-j))
    Dmat[indexvj,j] = 1
    Dmat[cbind(indexvj,(j+1):n)] = -1
  }
  
  
  deltamold = t(Dmat %*% betam)
  vm = matrix(0,p, npair)
  
  flag = 0
  niteration = 0
  for(m in 1:maxiter)
  {
    if(P ==1)
    {
      lamj = matrix(lamj,1,1)
    }
    
    ## revised admm 
    resadmm = revisedadmm(Bmlist,ylist,n,p, npair,Dmat,deltamold, vm,
                   lam, nu, gam,
                   theta, lamj, sig2)
    betam = resadmm$betam
    deltam = resadmm$deltam
    betadiff = resadmm$betadiff
    vm = resadmm$vm

    
    ### remove mean part ###
    bsmean = rep(0, ntotal)
    for(i in 1:n)
    {
      indi = ind == uind[i]
      bsmean[indi] = Bmlist[[i]] %*% betam[i,]
    }

    #### newton
    datasub = data.frame(ind = ind, obs = y - bsmean, time = tm) ## substract mean
    datalist = fpca.format(datasub)

    ## auxillary of newton
    phiaux = apply(matrix(1:n),MARGIN=1,Phi.aux,Bmlist = Bmlist,datalist=datalist)


    newton.result = Newton.New(theta,phiaux,sqrt(sig2),lamj,datalist,n,sl.v,max.step,tolnt,condtol)


    step = newton.result[[7]]
    likevalue<-newton.result[[1]][step] # -2*loglikelihood over all observations
    theta = as.matrix(newton.result[[2]][,,step],p,P)
    lamj = newton.result[[3]][step,]
    sig.up = newton.result[[4]][step]
    # gradB<-newton.result[[5]][,,step]
    # gradL<-newton.result[[6]][step,]
    # cond<-newton.result[[8]]
    # error.c<-newton.result[[9]]
    
    #cv.score = CV(theta,phiaux,sig.up,lamj,datalist,n)

    sig2 = sig.up^2
    
    norm1 = norm(betadiff,"F")
    norm2 = norm(deltam,"F")

    tolpri = tolabs*sqrt(npair*p) + tolrel*max(norm1,norm2)
    toldual = tolabs*sqrt(n * p) + tolrel * norm(vm %*% Dmat, "F")

    rm = norm(betadiff - deltam, "F")
    sm = nu * norm((deltam - deltamold)%*%Dmat, "F")

    deltamold = deltam

    niteration = niteration + 1

    if(rm <= tolpri & sm <= toldual)
    {
      break
    }

  }
  
  if(niteration == maxiter)
  {
    flag = 1
  }
  
  
  
  groupest = getgroup(deltam, n)
  ugroupest = unique(groupest)
  ng = length(ugroupest)
  
  betaavg = matrix(0,ng,p)
  
  for(j in 1:ng)
  {
    indexj = groupest == ugroupest[j]
    nj = sum(indexj)
    betaavg[j,] = colMeans(betam[indexj,,drop = FALSE])
  }
  
  betaest = betaavg[groupest,]
  
  
  #### mean function estimate #####
  meanfunest = rep(0, ntotal)
  for(i in 1:n)
  {
    indi = ind == uind[i]
    meanfunest[indi] = Bm[indi,] %*% betaest[i,]
  }
  
  
  #### re order lam and theta
  indlam = sort(lamj,index.return= TRUE, decreasing = TRUE)$ix
  
  lamj = lamj[indlam]
  theta = theta[,indlam,drop = FALSE]
  
  res =  list(betam = betam, betaest = betaest, betaavg = betaavg, knots = knotsall,
              obasisobj = obasisobj,
              meanfunest = meanfunest, groupest = groupest, likevalue = likevalue,
          sig2 = sig2, theta = theta, lamj = lamj,
          deltam = deltam, rm = rm, sm = sm,
          tolpri = tolpri, toldual = toldual, niteration = niteration, flag = flag, ntotal = ntotal)
  return(res)
}