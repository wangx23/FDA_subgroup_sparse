load("~/Research/FDA_subgroup_sparse/result/result_x_3cluster_ncl100lams2.RData")
names(result_x_3cluster_ncl100lams2[[1]])

nres100lams2 = length(result_x_3cluster_ncl100lams2)
aribest100lams2 = matrix(0, nres100lams2,4)
aring3P2 = matrix(0,nres100lams2,4 )
for(i in 1:nres100lams2)
{
  aribest100lams2[i,1] = max(result_x_3cluster_ncl100lams2[[i]]$res_median$arim)
  aribest100lams2[i,2] = max(result_x_3cluster_ncl100lams2[[i]]$res_kmeans$arim)
  aribest100lams2[i,3] = max(result_x_3cluster_ncl100lams2[[i]]$res_pam$arim)
  aribest100lams2[i,4] = max(result_x_3cluster_ncl100lams2[[i]]$res_ydist$arim)
  
  aring3P2[i,1] = result_x_3cluster_ncl100lams2[[i]]$res_median$arim[ which(result_x_3cluster_ncl100lams2[[i]]$res_median$ngest[,2]==3)[1],2]
  
  aring3P2[i,2] = result_x_3cluster_ncl100lams2[[i]]$res_kmeans$arim[ which(result_x_3cluster_ncl100lams2[[i]]$res_kmeans$ngest[,2]==3)[1],2]
  
  aring3P2[i,3] = result_x_3cluster_ncl100lams2[[i]]$res_pam$arim[ which(result_x_3cluster_ncl100lams2[[i]]$res_pam$ngest[,2]==3)[1],2]
  
  aring3P2[i,4] = result_x_3cluster_ncl100lams2[[i]]$res_ydist$arim[ which(result_x_3cluster_ncl100lams2[[i]]$res_ydist$ngest[,2]==3)[1],2]
}
boxplot(aring3P2)

### select a certain number of groups, and number of P and corresponding ARI
getari32 = function(obj, ng = 3, nP = 2)
{
  arivalue = obj$arim[obj$ngest[,nP]==3,nP][1]
  return(arivalue)
}


ari32 = matrix(0, nres100lams2, 4)
for(i in 1:nres100lams2)
{
 ari32[i,1] = getari32(result_x_3cluster_ncl100lams2[[i]]$res_median)
 ari32[i,2] = getari32(result_x_3cluster_ncl100lams2[[i]]$res_kmeans)
 ari32[i,3] = getari32(result_x_3cluster_ncl100lams2[[i]]$res_pam)
 ari32[i,4] = getari32(result_x_3cluster_ncl100lams2[[i]]$res_ydist)
}

boxplot(ari32)


bicgroupfun = function(obj, p, Cnpen = "log", n = 300, Ppen = "n", combined = FALSE, colsel = 1:3)
{
  
  nllmat = obj$nllmat
  nllmat[nllmat==-99] = 9999
  
  if(Cnpen == "log" )
  {
    Cn = log(n*p)
  }
  if(Cnpen =="loglog")
  {
    Cn = log(log(n*p))
  }
  
  if(Ppen == "n")
  {
    Pn = n
  }
  
  if(Ppen == "logn")
  {
    Pn = log(n)
  }
  
  logn = log(n)
  
  nconstraints = (1:3)*((1:3)+1)/2
  bicvalue1 = t(t(n*nllmat) + t(Cn* logn*(obj$ngest*p)) + Pn*((1:3)*(p-2) - nconstraints))
  
  if(combined)
  {
    if(Cnpen == "log" )
    {
      Cn = log(n*p + (1:3)*(p-2) - nconstraints)
    }
    if(Cnpen =="loglog")
    {
      Cn = log(log(n*p + (1:3)*(p-2) - nconstraints))
    }
    
    bicvalue1 = n*nllmat + t(Cn* logn*(t(obj$ngest*p) + (1:3)*(p-2) - nconstraints))
  }
  
  
  bicvalue1 = bicvalue1[,colsel, drop = FALSE]
  index = which(bicvalue1 == min(bicvalue1), arr.ind = TRUE)[1,]
  
  ngest0 = obj$ngest[,colsel, drop = FALSE]
  ariest0 = obj$arim[,colsel, drop = FALSE]
  
  ngest = ngest0[index[1],index[2]]
  Pest = colsel[index[2]]
  ariest = ariest0[index[1],index[2]]
  aribest = max(ariest0)
  
  out = c(ngest, Pest, ariest, aribest)
  return(out)
}

bic100 = matrix(0,nres100lams2 ,4)
for(i in 1:nres100lams2 )
{
  bic100[i,] = bicgroupfun(result_x_3cluster_ncl100lams2[[i]]$res_kmeans,
                           p = result_x_3cluster_ncl100lams2[[i]]$p,
                           Cnpen = "log", Ppen = "n", colsel = 1:3)
}

sum(bic100[,1]==3)
sum(bic100[,2]==2)



load("~/Research/FDA_subgroup_sparse/result/result_x_3cluster_ncl50lams2.RData")

nres50lams2 = length(result_x_3cluster_ncl50lams2[[1]])
aribest50lams2 = matrix(0, nres50lams2,4)
for(i in 1:nres50lams2)
{
  aribest50lams2[i,1] = max(result_x_3cluster_ncl50lams2[[i]]$res_median$arim)
  aribest50lams2[i,2] = max(result_x_3cluster_ncl50lams2[[i]]$res_kmeans$arim)
  aribest50lams2[i,3] = max(result_x_3cluster_ncl50lams2[[i]]$res_pam$arim)
  aribest50lams2[i,4] = max(result_x_3cluster_ncl50lams2[[i]]$res_ydist$arim)
}
