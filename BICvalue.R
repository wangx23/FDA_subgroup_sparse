###### BIC ####

## 
BICvalue = function(obj, c0 = 1)
{
  ntotal = obj$ntotal
  n = nrow(obj$betam)
  p = ncol(obj$betam)
  group = obj$groupest
  ng = length(unique(group))
  
  Cn = c0 * log(log(n*p))
  
  Pn = ncol(obj$theta) ### number of components
  nconstraints = Pn*(Pn+1)/2
  
  ### obj$likevalue is -2L average over n
  BICvalue = n*obj$likevalue + Cn*log(ntotal)*(ng*p) + 
    log(ntotal)*(Pn*p - nconstraints)
  
  if(obj$likevalue==-99){
    BICvalue = 99999
  }
  
  return(BICvalue)
}


BICvalue2 = function(obj, c0 = 1)
{
  ntotal = obj$ntotal
  n = nrow(obj$betam)
  p = ncol(obj$betam)
  group = obj$groupest
  ng = length(unique(group))
  
  Cn = c0 * log(log(n*p))
  
  Pn = ncol(obj$theta) ### number of components
  nconstraints = Pn*(Pn+1)/2
  
  ### obj$likevalue is -2L average over n
  BICvalue = n*obj$likevalue + Cn*log(n)*(ng*p) + 
    log(n)*(Pn*p - nconstraints)
  
  if(obj$likevalue==-99){
    BICvalue = 99999
  }
  
  return(BICvalue)
}



BICvalue3 = function(obj, c0 = 1)
{
  ntotal = obj$ntotal
  n = nrow(obj$betam)
  p = ncol(obj$betam)
  group = obj$groupest
  ng = length(unique(group))
  
  Cn = c0 * log(log(n*p))
  
  Pn = ncol(obj$theta) ### number of components
  nconstraints = Pn*(Pn+1)/2
  
  ### obj$likevalue is -2L average over n
  BICvalue = n*obj$likevalue + Cn*log(ntotal)*(ng*p) + 
    n*(Pn*p - nconstraints)
  
  if(obj$likevalue==-99){
    BICvalue = 99999
  }
  
  return(BICvalue)
}