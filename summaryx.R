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