#### x 1 and - 1

load("~/Research/FDA_subgroup_sparse/result/result_x_2cluster.RData")
load("~/Research/FDA_subgroup_sparse/result/result_x_2cluster_ydist.RData")
load("~/Research/FDA_subgroup_sparse/result/result_x_2cluster_ncl50.RData")
load("~/Research/FDA_subgroup_sparse/result/result_x_2cluster_ydist_ncl50.RData")



ari_x1_ncl100 = ari_x1_best_ncl100 = ari_x1_ncl50 = ari_x1_best_ncl50 = matrix(0,100,2)
lamest_x1_ncl100 = lamest_x1_ncl50 = array(0, dim = c(100,2,2))
Pest_x1_ncl100 = Pest_x1_ncl50 = array(0, dim = c(100,3,2))


for(i in 1:100)
{
  ari_x1_ncl100[i,1] = result_x_2cluster[[i]]$ari
  ari_x1_ncl100[i,2] = result_x_2cluster_ydist[[i]]$ari
  ari_x1_best_ncl100[i,1] = result_x_2cluster[[i]]$aribest
  ari_x1_best_ncl100[i,2] = result_x_2cluster_ydist[[i]]$aribest
  
  lamest_x1_ncl100[i,,1] = result_x_2cluster[[i]]$lamest
  lamest_x1_ncl100[i,,2] = result_x_2cluster_ydist[[i]]$lamest
  
  ari_x1_ncl50[i,1] = result_x_2cluster_ncl50[[i]]$ari
 # ari_x1_ncl50[i,2] = result_x_2cluster_ydist_ncl50[[i]]$ari
  
  lamest_x1_ncl50[i,,1] = result_x_2cluster_ncl50[[i]]$lamest
 # lamest_x1_ncl50[i,,2] = result_x_2cluster_ydist_ncl50[[i]]$lamest
}


summary(ari_x1_ncl100)
summary(ari_x1_best_ncl100)
summary(lamest_x1_ncl100[,,1])
summary(lamest_x1_ncl100[,,2])
summary(ari_x1_ncl50)
summary(lamest_x1_ncl50[,,1])


#### x 0.5 and - 0.5

load("~/Research/FDA_subgroup_sparse/result/result_x_2cluster_xvalue05.RData")
load("~/Research/FDA_subgroup_sparse/result/result_x_2cluster_xvalue05_ydist.RData")
load("~/Research/FDA_subgroup_sparse/result/result_x_2cluster_xvalue05_kmeans.RData")
load("~/Research/FDA_subgroup_sparse/result/result_x_2cluster_xvalue05_pam.RData")
load("~/Research/FDA_subgroup_sparse/result/result_x_2cluster_xvalue05_ncl50.RData")
load("~/Research/FDA_subgroup_sparse/result/result_x_2cluster_xvalue05_ydist_ncl50.RData")


ari_x05_ncl50 = ari_x05_best_ncl50 = matrix(0,100,2)
ari_x05_ncl100 = ari_x05_best_ncl100 = matrix(0,100,4)
lamest_x05_ncl50 = array(0, dim = c(100,2,2))
Pest_x05_ncl100 = Pest_x05_ncl50 = array(0, dim = c(100,3,2))
lamest_x05_ncl100 = array(0, dim = c(100,2,4))


for(i in 1:90)
{
  ari_x05_ncl100[i,1] = result_x_2cluster05[[i]]$ari
  ari_x05_ncl100[i,2] = result_x_2cluster05_ydist[[i]]$ari
  ari_x05_ncl100[i,3] = result_x_2cluster05_kmeans[[i]]$ari
  ari_x05_ncl100[i,4] = result_x_2cluster05_pam[[i]]$ari
  
  
  ari_x05_best_ncl100[i,1] = result_x_2cluster05[[i]]$aribest
  ari_x05_best_ncl100[i,2] = result_x_2cluster05_ydist[[i]]$aribest
  ari_x05_best_ncl100[i,3] = result_x_2cluster05_kmeans[[i]]$aribest
  ari_x05_best_ncl100[i,4] = result_x_2cluster05_pam[[i]]$aribest
  
  
  lamest_x05_ncl100[i,,1] = result_x_2cluster05[[i]]$lamest
  lamest_x05_ncl100[i,,2] = result_x_2cluster05_ydist[[i]]$lamest
  lamest_x05_ncl100[i,,3] = result_x_2cluster05_kmeans[[i]]$lamest
  lamest_x05_ncl100[i,,4] = result_x_2cluster05_pam[[i]]$lamest
  
  ari_x05_ncl50[i,1] = result_x_2cluster05_ncl50[[i]]$ari
  #ari_x05_ncl50[i,2] = result_x_2cluster05_ydist_ncl50[[i]]$ari
  
  lamest_x05_ncl50[i,,1] = result_x_2cluster05_ncl50[[i]]$lamest
  #lamest_x05_ncl50[i,,2] = result_x_2cluster05_ydist_ncl50[[i]]$lamest
}

summary(ari_x05_ncl100)
summary(ari_x05_best_ncl100)
summary(lamest_x05_ncl100[1:90,,1])
summary(lamest_x05_ncl100[1:90,,2])
summary(lamest_x05_ncl100[1:90,,3])
summary(lamest_x05_ncl100[1:90,,4])

xtable(summary(ari_x05_ncl100), digits = 4)
xtable(summary(cbind(lamest_x05_ncl100[,,1],lamest_x05_ncl100[,,2])), digits = 4)


summary(ari_x05_ncl50)
summary(ari_x05_ncl50)



#### 3cluster 
load("../result/result_x_3cluster_ncl100.RData")

arimbest = matrix(0,length(result_x_3cluster_ncl100),4)
arim2 = matrix(0,length(result_x_3cluster_ncl100),4)
for(i in 1:length(result_x_3cluster_ncl100))
{
  arimbest[i,1] = max(result_x_3cluster_ncl100[[i]]$res_median$arim)
  arimbest[i,2] = max(result_x_3cluster_ncl100[[i]]$res_kmeans$arim)
  arimbest[i,3] = max(result_x_3cluster_ncl100[[i]]$res_pam$arim)
  arimbest[i,4] = max(result_x_3cluster_ncl100[[i]]$res_ydist$arim)
  
  
  arim2[i,1] = result_x_3cluster_ncl100[[i]]$res_median$arim[result_x_3cluster_ncl100[[i]]$res_median$ngest[,2] == 3,2][1]
  arim2[i,2] = result_x_3cluster_ncl100[[i]]$res_kmeans$arim[result_x_3cluster_ncl100[[i]]$res_kmeans$ngest[,2] == 3,2][1]
  arim2[i,3] = result_x_3cluster_ncl100[[i]]$res_pam$arim[result_x_3cluster_ncl100[[i]]$res_pam$ngest[,2] == 3,2][1]
  arim2[i,4] = result_x_3cluster_ncl100[[i]]$res_ydist$arim[result_x_3cluster_ncl100[[i]]$res_ydist$ngest[,2] == 3,2][1]
  
}


#####  a function to generate different BIC and corresponding arim 
## Cn is the term for # of groups, Pn is the term for the number of components ##
## logn is a index for logn or logntotal 
bicgroupfun = function(obj, Cnpen = "log", n = 300, p =9, Ppen = "n", combined = FALSE)
{
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
  bicvalue1 = t(t(n*obj$nllmat) + t(Cn* logn*(obj$ngest*p)) + Pn*((1:3)*(p-2) - nconstraints))
  
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
    
    bicvalue1 = n*obj$nllmat + t(Cn* logn*(t(obj$ngest*p) + (1:3)*(p-2) - nconstraints))
  }
  index = which(bicvalue1 == min(bicvalue1), arr.ind = TRUE)
  
  ngest = obj$ngest[index]
  Pest = index[2]
  ariest = obj$arim[index]
  aribest = max(obj$arim)
  
  out = c(ngest, Pest, ariest, aribest)
  return(out)
}



BICdef = function(Cnpen, Ppen, combined = FALSE)
{
  fun1 = function(i)
  {
    out1 = bicgroupfun(result_x_3cluster_ncl100[[i]]$res_median, Cnpen = Cnpen, Ppen = Ppen,combined = combined)
    out2 = bicgroupfun(result_x_3cluster_ncl100[[i]]$res_kmeans,Cnpen =Cnpen, Ppen = Ppen,combined = combined)
    out3 = bicgroupfun(result_x_3cluster_ncl100[[i]]$res_pam,Cnpen =Cnpen, Ppen = Ppen,combined = combined)
    out4 = bicgroupfun(result_x_3cluster_ncl100[[i]]$res_ydist,Cnpen =Cnpen, Ppen = Ppen,combined = combined)
    
    out = cbind(out1, out2, out3, out4)
    return(out)
  }
  
  res = lapply(1:length(result_x_3cluster_ncl100),fun1)
  ngm = Pestm = ariestm = aribestm = matrix(0, length(result_x_3cluster_ncl100),4)
  
  for(i in 1:length(result_x_3cluster_ncl100))
  {
    ngm[i,] = res[[i]][1,]
    Pestm[i,] = res[[i]][2,]
    ariestm[i,] = res[[i]][3,]
    aribestm[i,] = res[[i]][4,]
  }
  
  
  outlist = list(ng = ngm, Pest = Pestm, ariest = ariestm, aribest = aribestm)
  return(outlist)

  
}

res1 = BICdef(Cnpen = "log", Ppen ="n")

summary(res1$ng)
apply(res1$ng,2,sd)
colSums(res1$ng==3)

summary(res1$Pest)
apply(res1$Pest,2,sd)
colSums(res1$Pest==2)

summary(res1$ariest)
summary(res1$aribest)


res2 = BICdef(Cnpen = "loglog", Ppen ="n")

summary(res2$ng)
apply(res2$ng,2,sd)
colSums(res2$ng==3)

summary(res2$Pest)
apply(res2$Pest,2,sd)
colSums(res2$Pest==2)

summary(res2$ariest)

res3 = BICdef(Cnpen = "log", Ppen ="log")

summary(res3$ng)
apply(res3$ng,2,sd)
colSums(res3$ng==3)

summary(res3$Pest)
apply(res3$Pest,2,sd)
colSums(res3$Pest==2)

res4 = BICdef(Cnpen = "loglog", Ppen ="log")

summary(res4$ng)
apply(res4$ng,2,sd)
colSums(res4$ng==3)

summary(res4$Pest)
apply(res4$Pest,2,sd)
colSums(res4$Pest==2)


res5 = BICdef(Cnpen = "loglog", Ppen = "n", combined = TRUE)

summary(res5$ng)
apply(res5$ng,2,sd)
colSums(res5$ng==3)

summary(res5$Pest)
apply(res5$Pest,2,sd)
colSums(res5$Pest==2)


######  best ##### 




