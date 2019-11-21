##### test ###

#setwd("Google Drive/Research/FDA_subgroup_sparse/code/")
source("simdat.R")
source("initial.R")
source("EMgroup.R")
source("FDAsubgroup.R")
source("BICvalue.R")

sig2 = 0.1
lamj = c(0.2,0.1)
mvec = c(5,20)
ncl = 50

#### another easier test simulated based on knots ####
beta1 = rep(-1,6)
beta2 = rep(1,6)
betam12 = rbind(beta1,beta2)
knots = seq(0,1,length.out = 4)[2:3]



dat1 = simdat1(knots, betam12, sig2 = sig2,lamj = lamj,mvec = mvec,ncl = ncl)

ind = dat1$ind
tm = dat1$time
y = dat1$obs
group = unique(dat1[,1:2])[,1]
lamv = seq(0,20, by = 0.5)[-1]

res1 = FDAsubgroup(ind = ind,tm = tm,y = y,P = 2,knots = knots,
                   lam = 1,maxiter = 50)



###### curves simulation ###
dat = simdat(sig2 = 0.1,lamj = c(0.2,0.1),mvec = c(5,20),ncl = 70)

knots = seq(0,1,length.out = 6)[2:5]
boundary = c(0,1)

ind = dat$ind
tm = dat$time
y = dat$obs
group = unique(dat[,1:2])[,1]
lamv = seq(0,20, by = 0.5)[-1]

res1 = FDAsubgroup(ind = dat$ind,tm = dat$time,y = dat$obs,P = 2,
                  knots = seq(0,1,length.out = 6)[2:5],
                   lam = 1,maxiter = 50)
res2 = FDAsubgroup(ind = dat$ind,tm = dat$time,y = dat$obs,P = 2,
                  knots = seq(0,1,length.out = 6)[2:5],
                  lam = 0.5,maxiter = 50)

res3 = FDAsubgroup(ind = dat$ind,tm = dat$time,y = dat$obs,P = 1,
                   knots = seq(0,1,length.out = 6)[2:5],
                   lam = 0.5,maxiter = 50)


res4 = FDAsubgroup(ind = dat$ind,tm = dat$time,y = dat$obs,P = 3,
                   knots = seq(0,1,length.out = 6)[2:5],
                   lam = 0.5,maxiter = 50)




res$groupest

BICvalue(res1)
BICvalue(res2)
BICvalue(res3)
BICvalue(res4)

lamv = seq(0.2,1,by = 0.05)
BICvec = rep(0,length(lamv))
for(j in 1:length(lamv))
{
  resi = FDAsubgroup(ind = dat$ind,tm = dat$time,y = dat$obs,P = 2,
                     knots = seq(0,1,length.out = 6)[2:5],
                     lam = lamv[j],maxiter = 50)
  BICvec[j] = BICvalue(resi)
  
}



res = FDAsubgroup(ind = dat$ind,tm = dat$time,y = dat$obs,P = 2,
                   knots = seq(0,1,length.out = 6)[2:5],
                   lam = lamv[which.min(BICvec)],maxiter = 50)



#### test for both lambda and P ####

lamv = seq(0.2,1,by = 0.05)
BICm = matrix(0,length(lamv), 3)
for(j in 1:length(lamv))
{
  for(Pv in 1:3)
  {
    resi = FDAsubgroup(ind = dat$ind,tm = dat$time,y = dat$obs,P = Pv,
                       knots = seq(0,1,length.out = 6)[2:5],
                       lam = lamv[j],maxiter = 50)
    BICm[j,Pv] = BICvalue(resi)
  }
}

inds = which(BICm == min(BICm), arr.ind = TRUE)

res = FDAsubgroup(ind = dat$ind,tm = dat$time,y = dat$obs,P = inds[2],
                  knots = seq(0,1,length.out = 6)[2:5],
                  lam = lamv[inds[1]],maxiter = 50)








