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