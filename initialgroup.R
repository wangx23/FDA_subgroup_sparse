##### a function to conduct kmeans such that if one group only has one subject, it will be assigned to the nearest group. 
## x is the original observations ##
## K0 is the original number of clusters
initialgroup_ks = function(x, K0, ming = 1)
{
  fitks = kmeans(x,centers = K0)
  groupks = fitks$cluster
  centersks = fitks$centers
  distmat = as.matrix(dist(centersks))
  diag(distmat) = 9999
  
  ngobs = table(groupks)
  grouporder = 1:K0
  
  ### for cluster with only one observation, replace it by the closest group ###
  while(min(ngobs)<=ming)
  {
    index1 = as.numeric(names(ngobs)[ngobs<=ming])[1]
    grouporder[index1] = grouporder[which.min(distmat[index1,])]
    
    groupks[groupks==index1] = grouporder[index1]
    ngobs = table(groupks)
    
  }
  groupks = as.numeric(as.factor(groupks))
  return(groupks)
  
}
  
  
  
### x can be the original observations or the dissimilarity matrix
initialgroup_pam = function(x, K0, diss = TRUE, ming = 1)
{
  if(diss)
  {
    fitpam = pam(x, K0, diss = TRUE)
    grouppam = fitpam$clustering
    id.med = fitpam$id.med
    distmat = x
    diag(distmat)  = 9999
  }else{
    fitpam = pam(x, K0, diss = FALSE)
    grouppam = fitpam$clustering
    id.med = fitpam$id.med
    distmat = as.matrix(dist(x))
    diag(distmat)  = 9999
  }
  
  ngobs = table(grouppam)
  grouporder = 1:K0
  
  while(min(ngobs) <= ming)
  {
    index1 = as.numeric(names(ngobs)[ngobs<=ming])[1]
    
    grouporder[index1] = grouporder[which.min(distmat[grouppam==index1,id.med])]
    grouppam[grouppam==index1] = grouporder[index1]
    ngobs = table(grouppam)
    
  }
  
  grouppam = as.numeric(as.factor(grouppam))
  
 return(grouppam)
}