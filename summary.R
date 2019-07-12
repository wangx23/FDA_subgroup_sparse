load("~/Research/FDA_subgroup_sparse/result/resulttest1.rdata")
result1 = do.call("rbind",resulttest1)

summary(result1)

sqrt(mean((result1[,2] - 0.1)^2))
