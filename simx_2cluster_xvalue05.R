
source("simx_2cluster.R")

funlist = list(Vectorize(function(x){4*(x/50)^(0.2)}),
               Vectorize(function(x){-4*(x/50)^(0.2)}))

eigenlist = list(Vectorize(function(x){sqrt(2)*sin(pi*x)}),
                 Vectorize(function(x){ sqrt(2)*cos(pi*x)}))

xlist = list(nx = 2, meanx = 0, sdx = 1, etag = matrix(c(-0.5,-0.5,0.5,0.5),nrow = 2))
sig200 = 0.04
lamj00 = c(0.2,0.1)
mvec00 = c(10,20)
ncl00 = 50
lamvec = seq(0.2,2,by = 0.1)


#### ncl = 50 #####
# subfunx2cluster = function(mm){
#   subfunx(mm, sig200 = sig200, lam00 = lamj00, mvec00 = mvec00,
#           ncl00 = 50, lamvec = lamvec,
#           funlist = funlist, eigenlist = eigenlist,xlist= xlist)
# }
# 
# cl <- makeCluster(24)
# registerDoParallel(cl)
# result_x_2cluster05_ncl50 <- foreach(mm=1:100,
#                              .packages=c("flexclust","orthogonalsplinebasis","plyr"),
#                              .errorhandling = "remove") %dopar%subfunx2cluster(mm)
# stopCluster(cl)
# save(result_x_2cluster05_ncl50,file = "../result/result_x_2cluster_xvalue05_ncl50.RData")
# 
# 
# subfunx2cluster2 = function(mm){
#   subfunx(mm, sig200 = sig200, lam00 = lamj00, mvec00 = mvec00,
#           ncl00 = 50, lamvec = lamvec, 
#           funlist = funlist, eigenlist = eigenlist,xlist= xlist, method = "ydist")
# }
# 
# cl <- makeCluster(24)  
# registerDoParallel(cl)  
# result_x_2cluster05_ydist_ncl50 <- foreach(mm=1:100,
#                                .packages=c("flexclust","orthogonalsplinebasis","plyr","cluster"),
#                                .errorhandling = "remove") %dopar%subfunx2cluster2(mm)
# stopCluster(cl) 
# save(result_x_2cluster05_ydist_ncl50,file = "../result/result_x_2cluster_xvalue05_ydist_ncl50.RData")




#### ncl = 100 #####
# subfunx2cluster = function(mm){
#   subfunx(mm, sig200 = sig200, lam00 = lamj00, mvec00 = mvec00,
#           ncl00 = 100, lamvec = lamvec,
#           funlist = funlist, eigenlist = eigenlist,xlist= xlist)
# }
# 
# cl <- makeCluster(24)
# registerDoParallel(cl)
# result_x_2cluster05 <- foreach(mm=1:100,
#                                      .packages=c("flexclust","orthogonalsplinebasis","plyr"),
#                                      .errorhandling = "remove") %dopar%subfunx2cluster(mm)
# stopCluster(cl)
# save(result_x_2cluster05,file = "../result/result_x_2cluster_xvalue05.RData")


# subfunx2cluster2 = function(mm){
#   subfunx(mm, sig200 = sig200, lam00 = lamj00, mvec00 = mvec00,
#           ncl00 = 100, lamvec = lamvec, 
#           funlist = funlist, eigenlist = eigenlist,xlist= xlist, method = "ydist")
# }
# 
# cl <- makeCluster(24)  
# registerDoParallel(cl)  
# result_x_2cluster05_ydist <- foreach(mm=1:100,
#                                            .packages=c("flexclust","orthogonalsplinebasis","plyr","cluster"),
#                                            .errorhandling = "remove") %dopar%subfunx2cluster2(mm)
# stopCluster(cl) 
# save(result_x_2cluster05_ydist,file = "../result/result_x_2cluster_xvalue05_ydist.RData")



####### kmeans ###
subfunx2cluster3 = function(mm){
  subfunx(mm, sig200 = sig200, lam00 = lamj00, mvec00 = mvec00,
          ncl00 = 100, lamvec = lamvec,
          funlist = funlist, eigenlist = eigenlist,xlist= xlist, method = "kmeans")
}

cl <- makeCluster(24)
registerDoParallel(cl)
result_x_2cluster05_kmeans <- foreach(mm=1:100,
                               .packages=c("flexclust","orthogonalsplinebasis","plyr"),
                               .errorhandling = "remove") %dopar%subfunx2cluster3(mm)
stopCluster(cl)
save(result_x_2cluster05_kmeans,file = "../result/result_x_2cluster_xvalue05_kmeans.RData")


subfunx2cluster4 = function(mm){
  subfunx(mm, sig200 = sig200, lam00 = lamj00, mvec00 = mvec00,
          ncl00 = 100, lamvec = lamvec, 
          funlist = funlist, eigenlist = eigenlist,xlist= xlist, method = "pam")
}

cl <- makeCluster(24)  
registerDoParallel(cl)  
result_x_2cluster05_pam <- foreach(mm=1:100,
                                     .packages=c("flexclust","orthogonalsplinebasis","plyr","cluster"),
                                     .errorhandling = "remove") %dopar%subfunx2cluster4(mm)
stopCluster(cl) 
save(result_x_2cluster05_pam,file = "../result/result_x_2cluster_xvalue05_pam.RData")