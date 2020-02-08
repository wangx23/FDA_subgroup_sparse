##### refit based on netwon algorithm with known group information ###

# refitFDA = function()
# {
#   
# }




##### refitFDA with x grouping #####

## group0 is the given group information 
refitFDAX = function(ind, tm, x, y, P = 2, group0, knots, boundary = c(0,1),
                     maxiter = 500, tol = 1e-3,
                     sl.v=rep(0.5,10), max.step = 20, 
                     tolnt=1e-3,condtol=1e+10,
                     seed = 2118)
{
  
  ntotal = length(y)
  knotsall = c(rep(boundary[1],4),knots, rep(boundary[2],4))
  obasisobj = OBasis(knotsall)
  Bm = evaluate(obasisobj,tm) ## orthogonal Bspline basis 
  Xtilde = cbind(x, Bm) ## x + orthogonal Bspline basis 
  
  
  p = ncol(Xtilde)
  uind = unique(ind)
  n = length(uind)
  
  ng = length(unique(group0))
  Ux = matrix(0, ntotal, ng*p)
  
  ylist = list()
  Bmlist = list()
  
  for(i in 1:n)
  {
    groupi = group0[i]
    indi = ind == uind[i]
    Ux[indi,((groupi - 1)*p+1):(groupi*p)] = Xtilde[indi,]
    ylist[[i]] = y[indi]
    Bmlist[[i]] = Bm[indi,]
  }
  
  est = solve(t(Ux)%*%Ux) %*% t(Ux) %*%y  #### initial value without covaraince 
  
  alpm = matrix(est, ng, p, byrow = TRUE)
  
  yresid = rep(0, ntotal)
  for(i in 1:n)
  {
    indi = ind == uind[i]
    yresid[indi] = ylist[[i]] - Xtilde[indi,]%*%alpm[group0[i],]
  }
  
  if(length(sl.v)<max.step){
    sl.v<-c(sl.v,rep(1,max.step-length(sl.v)))
  }else{
    sl.v<-sl.v[1:max.step]
  }
  
  ### intitial based on EM
  resem0 = EMcov(ind = ind,tm = tm,yresid = yresid,group0 = group0,knots = knots,P = 2)
  theta = resem0$theta
  sig2 = resem0$sig2
  lamj = resem0$lamj
  
  
  
  for(m in 1:maxiter)
  {
    #### update covariance parameters####
    for(i in 1:n)
    {
      indi = ind == uind[i]
      yresid[indi] = ylist[[i]] - Xtilde[indi,]%*%alpm[group0[i],]
    }
    
    datasub = data.frame(ind = ind, obs = yresid, time = tm) ## substract mean
    datalist = fpca.format(datasub)
    phiaux = apply(matrix(1:n),MARGIN=1,Phi.aux,Bmlist = Bmlist,datalist=datalist)
    
    newton.result = Newton.New(theta,phiaux,sqrt(sig2),lamj,datalist,n,
                               sl.v,max.step,tolnt,condtol)
    
    step = newton.result[[7]]
    likevalue<-newton.result[[1]][step] # -2*loglikelihood over all observations
    theta = as.matrix(newton.result[[2]][,,step],p,P)
    lamj = newton.result[[3]][step,]
    sig.up = newton.result[[4]][step]
    sig2 = sig.up^2
    
    #### update new coefficients #####
    UtU = matrix(0, ng*p, ng*p)
    Uty = rep(0, ng*p)
    for(i in 1:n)
    {
      indi = ind == uind[i]
      Bmi = Bmlist[[i]]
      Uxi = Ux[indi,]
      ni = nrow(Bmi)
      Bthetai = Bmi %*% theta
      Sigmainvi  = solve( diag(ni) + 1/sig2*Bthetai %*% diag(lamj) %*% t(Bthetai))
      
      UtU = UtU + t(Uxi) %*% Sigmainvi %*% Uxi
      
      Uty = Uty + t(Uxi) %*% Sigmainvi %*% ylist[[i]]
    }
    
    alpmnew = matrix(solve(UtU) %*%Uty, ng, p, byrow = TRUE)
    diffalp = max(abs(alpmnew - alpm))
    
    alpm = alpmnew
    ###### stopping ####
    
    if(diffalp <= tol)
    {break}
    
 
  }
  
  ##### covariance matrix ##### 
  
  covmat = sig2*solve(UtU)
  
  SEmat = matrix(sqrt(diag(covmat)), ng, p, byrow = TRUE)
  
  
  outlist = list(alpha = alpm, theta = theta, sig2 = sig2, lamj = lamj, SEmat = SEmat,likevalue = likevalue)
  return(outlist)
  
  
}




#### refit x independent ####
refitINDX = function(ind, tm, x, y, group0, knots, boundary = c(0,1))
{
  
  ntotal = length(y)
  knotsall = c(rep(boundary[1],4),knots, rep(boundary[2],4))
  obasisobj = OBasis(knotsall)
  Bm = evaluate(obasisobj,tm) ## orthogonal Bspline basis 
  Xtilde = cbind(x, Bm) ## x + orthogonal Bspline basis 
  
  
  p = ncol(Xtilde)
  uind = unique(ind)
  n = length(uind)
  
  ng = length(unique(group0))
  Ux = matrix(0, ntotal, ng*p)
  
  for(i in 1:n)
  {
    groupi = group0[i]
    indi = ind == uind[i]
    Ux[indi,((groupi - 1)*p+1):(groupi*p)] = Xtilde[indi,]
  }
  
  est = solve(t(Ux)%*%Ux) %*% t(Ux) %*%y  #### initial value without covaraince 
  
  alpm = matrix(est, ng, p, byrow = TRUE)
  
  betam = alpm[group0,]
  
  return(betam)
  
}
