
source("subfunx_all.R")

funlist = list(Vectorize(function(x){4*(x -0.5)^2+1}),
               Vectorize(function(x){2.5*exp(-25*(x-0.25)^2) + 2*exp(-50*(x-0.75)^2)})
)


eigenlist = list(Vectorize(function(x){sqrt(2)*sin(pi*x)}),
                 Vectorize(function(x){ sqrt(2)*cos(pi*x)}))

xlist = list(nx = 2, meanx = 0, sdx = 1, etag = matrix(c(-0.5,-0.5,0.5,0.5),nrow = 2))

sig200 = 0.04
lamj00 = c(0.1, 0.05)
mvec00 = c(10,20)

lamvec1 = seq(0.3,0.7, by = 0.0)
#lamvec2 = seq(0.2,0.6, by = 0.02)
lamvec2 = c(seq(0.2,0.25, by = 0.0025),seq(0.26,0.6, by = 0.1))
#lamvec3 = c(seq(0.15,0.2,by = 0.005),seq(0.25,0.7, by = 0.05))
lamvec3 = seq(0.15,0.7, by = 0.02)
#lamvec4 = c(seq(0.3,0.5, by = 0.02), seq(0.501,0.55,by = 0.001))
lamvec4 = seq(0.2,0.6, by = 0.02)

lamvec1 = 0.2
#lamvec2 = 0.1
lamvec3 = 0.1
lamvec4 = 0.1



# subfunx3cluster = function(mm){
#   subfunx_all(mm, sig200 = sig200, lam00 = lamj00, mvec00 = mvec00,
#           ncl00 = 100, lamvec = lamvec,
#           funlist = funlist1, eigenlist = eigenlist,xlist= xlist1, K0 = 15)
# }

subfunx2cluster = function(mm){
  subfunx_all(mm, sig200 = sig200, lam00 = lamj00, mvec00 = mvec00,
              ncl00 = 50, 
              lamvec1 = lamvec1, lamvec2 = lamvec2, lamvec3 = lamvec3, lamvec4 = lamvec4,
              funlist = funlist, eigenlist = eigenlist,xlist= xlist, 
              K0 = 15, ming = 2, max.step = 5)
}




res1 = subfunx2cluster(1)

library(doParallel)
cl <- makeCluster(24)
registerDoParallel(cl)
result_x_2cluster_ncl100lam15 <- foreach(mm=1:150,
                                         .packages=c("flexclust","orthogonalsplinebasis","plyr","cluster"),
                                         .errorhandling = "remove") %dopar%subfunx2cluster(mm)
stopCluster(cl)
save(result_x_2cluster_ncl100lam15,file = "../result/result_x_2cluster_ncl100lam15.RData")
