#### x 1 and - 1

load("~/Research/FDA_subgroup_sparse/result/result_x_2cluster.RData")
load("~/Research/FDA_subgroup_sparse/result/result_x_2cluster_ydist.RData")


ari_x1 = matrix(0,100,2)
lamest_x1 = array(0, dim = c(100,2,2))
Pest_x1 = array(0, dim = c(100,3,2))

for(i in 1:100)
{
  ari_x1[i,1] = result_x_2cluster[[i]]$ari
  ari_x1[i,2] = result_x_2cluster_ydist[[i]]$ari

  lamest_x1[i,,1] = result_x_2cluster[[i]]$lamest
  lamest_x1[i,,2] = result_x_2cluster_ydist[[i]]$lamest
}


load("~/Research/FDA_subgroup_sparse/result/result_x_2cluster_ncl50.RData")
load("~/Research/FDA_subgroup_sparse/result/result_x_2cluster_ydist_ncl50.RData")


load("~/Research/FDA_subgroup_sparse/result/result_x_2cluster_xvalue05.RData")

ari_xvalue1 = rep(0,100)
ari_xvalue05 = rep(0,100)
lamest_xvalue1 = lamest_xvalue05 = matrix(0,100,2)
Pest_xvalue1 = Pest_xvalue05 = matrix(0,100,3)

for(i in 1:100)
{
  ari_xvalue1[[i]] = result_x_2cluster[[i]]$ari
  ari_xvalue05[[i]] = result_x_2cluster05[[i]]$ari
  
  lamest_xvalue1[i,] = result_x_2cluster[[i]]$lamest
  lamest_xvalue05[i,] = result_x_2cluster05[[i]]$lamest
  
  Pest_xvalue1[i,] = result_x_2cluster[[i]]$Pest
  Pest_xvalue05[i,] = result_x_2cluster05[[i]]$Pest
}

summary(lamest_xvalue1)
summary(lamest_xvalue05)

summary(Pest_xvalue1)
summary(Pest_xvalue05)

summary(ari_xvalue05)
summary(ari_xvalue1)

