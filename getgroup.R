##### obtain group based on deltam ###

getgroup = function(deltam, n, tol = 1e-2)
{
  p = nrow(deltam)
  b2value =sqrt(colMeans(deltam^2))
  b2value[b2value <= tol] = 0
  
  d2 = matrix(0, n, n)
  for(j in 1:(n-1))
  {
    indexj1 = (2*n -j)*(j-1)/2 + 1
    indexj2 = indexj1 + n - j - 1
    d2[(n - indexj2 + indexj1):n,j] = b2value[indexj1:indexj2]
  }
  d2 = t(d2) + d2
  
  
  ngj = 1:n
  groupest = rep(0,n)
  j = 0
  
  while (length(ngj) >0)
  {
    j = j + 1
    gj = (1:n)[d2[ngj[1],] ==0]
    indexj = ngj %in% gj
    gj = ngj[indexj]
    ngj = ngj[!indexj]
    groupest[gj] = j * rep(1,length(gj))
  }


  return(groupest)
  
}




