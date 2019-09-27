load("~/Research/FDA_subgroup_sparse/result/resultcompgrids1.RData")

load("~/Research/FDA_subgroup_sparse/result/resultcompgrids3.RData")

output_array3 = output_array2 = array(0, dim = c(5,3,100))
imsemat3 = matrix(0,150,3)
imse_eign_mat3= matrix(0,150,2)
mse_lam_mat3 = matrix(0,150,2)

imsemat2 = matrix(0,100,3)
imse_eign_mat2= matrix(0,100,2)
mse_lam_mat2 = matrix(0,100,2)

maxvalue = rep(0,100)

for(i in 1:100)
{

  output_array3[,,i] = resultcomp3[[i]]$outputmat
  imsemat3 = imsemat3 + resultcomp3[[i]]$isemeanmat
  imse_eign_mat3[i,]= resultcomp3[[i]]$ise_eig
  mse_lam_mat3[i,] = resultcomp3[[i]]$mse_lamj
  
  output_array2[,,i] = resultcomp1[[i]]$outputmat
  imsemat2 = imsemat2 + resultcomp1[[i]]$isemeanmat
  imse_eign_mat2[i,]= resultcomp1[[i]]$ise_eig
  mse_lam_mat2[i,] = resultcomp1[[i]]$mse_lamj
  
  maxvalue[i] = max(resultcomp1[[i]]$isemeanmat[,1])
}



apply(output_array3,c(1,2),mean)
apply(output_array3,c(1,2),sd)


boxplot(imsemat3/100)
colMeans(imse_eign_mat3[rowSums(imse_eign_mat3)!=0,])
sqrt(colMeans(mse_lam_mat3[rowSums(mse_lam_mat3)!=0,]))



apply(output_array2,c(1,2),mean)
apply(output_array2,c(1,2),sd)


boxplot(imsemat2/100)
colMeans(imse_eign_mat2[rowSums(imse_eign_mat2)!=0,])
sqrt(colMeans(mse_lam_mat2[rowSums(mse_lam_mat2)!=0,]))


