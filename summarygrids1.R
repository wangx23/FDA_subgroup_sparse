#######  

load("../result/resultcompgrids1.RData")

resultcomp1[[1]]

output_array = array(0, dim = c(4,3,100))
imsemat = matrix(0,100,3)
imse_eign_mat= matrix(0,100,2)
mse_lam_mat = matrix(0,100,2)


for(m in 1:100)
{
  output_array[,,m] = resultcomp1[[m]]$outputmat
  imsemat = imsemat + resultcomp1[[m]]$isemeanmat
  imse_eign_mat[m,]= resultcomp1[[m]]$ise_eig
  mse_lam_mat[m,] = resultcomp1[[m]]$mse_lamj
}

apply(output_array,c(1,2), mean)
boxplot(imsemat/100)
colMeans(imse_eign_mat[rowSums(imse_eign_mat)!=0,])
sqrt(colMeans(mse_lam_mat[rowSums(mse_lam_mat)!=0,]))

