###### IMSE based on grids #####

##### FDA ####
#### ISE of mean functions (over grids, not all simulations) ####

ISEFDAmean = function(obj, group0,grids, funlist)
{
  obasisobj = obj$obasisobj
  betaest = obj$betaest
  nobs = length(group0)
  
  bmgrid = evaluate(res$obasisobj,grids)
  
  meanfungroup = bmgrid %*% t(obj$betaavg)
  
  ISEvec = rep(0, nobs)
  for(i in 1:nobs)
  {
    groupi = group0[i]
    meani = funlist[[groupi]](grids) ## true mean 
    ISEvec[i] = mean((meani - meanfungroup[,groupi])^2)
  }
  
  return(ISEvec)

}



### eigenfunctions ###
### match the sign of the max abs value in the true
### If the true number of components is correct ####
ISEFDAeig = function(obj, grids, funlist)
{
  thetaest = obj$theta
  
  obasisobj = obj$obasisobj
  bmgrid = evaluate(res$obasisobj,grids)
  
  eigenfunest = bmgrid%*%thetaest  
  eigenfunmat = sapply(1:length(funlist),function(x){eigenlist[[x]](grids)})  ### true
  indmax = apply(abs(eigenfunmat),2,which.max)
  signmax = sapply(eigenfunmat[cbind(indmax,1:2)],sign)
  
  signmaxest = sapply(eigenfunest[cbind(indmax,1:2)],sign)
  
  eigenfunest = sapply(1:ncol(eigenfunest), function(x){eigenfunest[,x]*signmax[x]*signmaxest[x]})
  
  ISEvec = colMeans((eigenfunest - eigenfunmat)^2)
  return(ISEvec)
}



#### ISE mean for IND #####
ISEINDmean = function(betaest, group0, grids, funlist, nknots =4, order = 4)
{
  bmx = bsplineS(grids, knots_eq3(timerange, k = order, m = nknots), norder = order)
  nobs = length(group0)
  
  ISEvec = rep(0, nobs)
  for(i in 1:nobs)
  {
    groupi = group0[i]
    meani = funlist[[groupi]](grids)
    meanesti = bmx %*% betaest[i,]
    ISEvec[i] = mean((meani - meanesti)^2)
  }
  return(ISEvec)
}



#### ISE for JS ######

ISEJSmean = function(obj, group0, grids, funlist)
{
  centers = fclust.curvepred(fit=obj)$meancurves
  nobs = length(group0)
  
  ISEvec = rep(0, nobs)
  for(i in 1:nobs)
  {
    groupi = group0[i]
    meani = funlist[[groupi]](grids) ## true mean 
    ISEvec[i] = mean((meani - centers[,groupi])^2)
  }
  return(ISEvec)
}