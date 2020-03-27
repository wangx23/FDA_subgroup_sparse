load("~/Research/FDA_subgroup_sparse/result/result_x_3cluster_ncl100lams2.RData")
names(result_x_3cluster_ncl100lams2[[1]])

nres100lams2 = length(result_x_3cluster_ncl100lams2[[1]])
aribest100lams2 = matrix(0, nres100lams2,4)
for(i in 1:nres100lams2)
{
  aribest100lams2[i,1] = max(result_x_3cluster_ncl100lams2[[i]]$res_median$arim)
  aribest100lams2[i,2] = max(result_x_3cluster_ncl100lams2[[i]]$res_kmeans$arim)
  aribest100lams2[i,3] = max(result_x_3cluster_ncl100lams2[[i]]$res_pam$arim)
  aribest100lams2[i,4] = max(result_x_3cluster_ncl100lams2[[i]]$res_ydist$arim)
}


load("~/Research/FDA_subgroup_sparse/result/result_x_2cluster_ncl100lam15.RData")

nres100lams15 = length(result_x_3cluster_ncl100lams15[[1]])
aribest100lams15 = matrix(0, nres100lams15,4)
for(i in 1:nres100lams15)
{
  aribest100lams15[i,1] = max(result_x_3cluster_ncl100lams15[[i]]$res_median$arim)
  aribest100lams15[i,2] = max(result_x_3cluster_ncl100lams15[[i]]$res_kmeans$arim)
  aribest100lams15[i,3] = max(result_x_3cluster_ncl100lams15[[i]]$res_pam$arim)
  aribest100lams15[i,4] = max(result_x_3cluster_ncl100lams15[[i]]$res_ydist$arim)
}
