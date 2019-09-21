##### test for server #####

fun1 = function(x)
{
  set.seed(x)
  mean(rnorm(100000))
}

cl <- makeCluster(24)  
registerDoParallel(cl)  
test <- foreach(mm=1:100) %dopar%  fun1(mm)
stopCluster(cl) 
save(test,file = "../result/test.RData")