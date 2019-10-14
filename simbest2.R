source("subfun_best.R")

lam00 = c(0.1,0.05)
sig2 = 0.05
mvec = c(5,20)
ncl= 50

lamvec = exp(seq(-1.5,-0.9, length.out = 50))
lamvec2 = exp(seq(-0.5,0.3,length.out = 50))


# t1 = Sys.time()
# res1 = subfun_best(mm = 77,sig2 = sig2,lam00 = lam00,mvec = mvec,ncl = ncl,lamvec = lamvec,lamvec2 = lamvec2)
# t2 = Sys.time()

subfun2 = function(mm)
{
  subfun_best(mm,sig2 = sig2,lam00 = lam00,mvec = mvec,ncl = ncl,lamvec = lamvec,lamvec2 = lamvec2)
}

cl <- makeCluster(24)  
registerDoParallel(cl)  
resultcomp3_best2 <- foreach(mm=1:120,
                            .packages=c("flexclust","orthogonalsplinebasis","plyr","fda","Lclust","igraph","mcclust"),.errorhandling = "remove") %dopar%  subfun2(mm)
stopCluster(cl) 
save(resultcomp3_best2,file = "../result/resultcompgrids3_best2.RData")