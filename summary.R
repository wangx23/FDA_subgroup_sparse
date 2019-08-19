load("~/Research/FDA_subgroup_sparse/result/resulttest1.rdata")
result1 = do.call("rbind",resulttest1)

summary(result1)

sqrt(mean((result1[,2] - 0.1)^2))


load("~/Research/FDA_subgroup_sparse/result/resulttestP1.rdata")
resultP1 = do.call("rbind",resulttestP1)
summary(resultP1)

sqrt(mean((resultP1[,2] - 0.1)^2))

summary(resultP1[,3])
sum(resultP1[,3]==2)
summary(resultP1[,1])
sd(resultP1[,1])
