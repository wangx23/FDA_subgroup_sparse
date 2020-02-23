
source("subfunx_all.R")

funlist1 = list(Vectorize(function(x){cos(2*pi*x)}),
                Vectorize(function(x){1-2*exp(-6*x)}),
                Vectorize(function(x){-1.5*x})
)

eigenlist = list(Vectorize(function(x){sqrt(2)*sin(pi*x)}),
                 Vectorize(function(x){ sqrt(2)*cos(pi*x)}))

xlist1 = list(nx = 2, meanx = 0, sdx = 1, etag = matrix(c(-0.5,-0.5,0.5,0.5,1,1),nrow = 2))
sig200 = 0.04
# lamj00 = c(0.2,0.1)
lamj00 = c(0.15, 0.1)
mvec00 = c(10,20)
lamvec = seq(0.2,0.8,by = 0.025)


subfunx3cluster = function(mm){
  subfunx_all(mm, sig200 = sig200, lam00 = lamj00, mvec00 = mvec00,
          ncl00 = 100, lamvec = lamvec,
          funlist = funlist1, eigenlist = eigenlist,xlist= xlist1, 
          K0 = 15,ming = 1, max.step = 5)
}

# subfunx3cluster = function(mm){
#   subfunx_all(mm, sig200 = sig200, lam00 = lamj00, mvec00 = mvec00,
#               ncl00 = 50, lamvec = lamvec,
#               funlist = funlist1, eigenlist = eigenlist,xlist= xlist1, 
#               K0 = 12, ming = 2, max.step = 5)
# }



# t1 = Sys.time()
# res1 = subfunx3cluster(1)
# t2 = Sys.time()
# 
# t1 = Sys.time()
# res2 = subfunx3cluster(30)
# t2 = Sys.time()


library(doParallel)
cl <- makeCluster(24)
registerDoParallel(cl)
result_x_3cluster_ncl100lams2 <- foreach(mm=1:150,
                             .packages=c("flexclust","orthogonalsplinebasis","plyr","cluster"),
                             .errorhandling = "remove") %dopar%subfunx3cluster(mm)
stopCluster(cl)
save(result_x_3cluster_ncl100lams2,file = "../result/result_x_3cluster_ncl100lams2.RData")



