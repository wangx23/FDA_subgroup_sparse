library(gridExtra)
load("Google Drive File Stream/My Drive/Research/FDA_subgroup_sparse/result/resultcomp1.rdata")
matcomp1 = do.call("rbind",resultcomp1)

result1 = matcomp1[seq(1,to = 200,by = 2),]
result2 = matcomp1[seq(2,to = 200,by = 2),]

hist(result1[,1])
hist(result2[,1])


pdf("Google Drive File Stream/My Drive/Research/FDA_subgroup_sparse/docs/figures/aricomp1.pdf",width = 6,height = 3)
par(mfrow = c(1,2))
boxplot(cbind(FDA = result1[,1], IND = result2[,1]))
boxplot(cbind(FDA = result1[,1], IND = result2[,1])[result2[,1]!=0,])
dev.off()

pdf("Google Drive File Stream/My Drive/Research/FDA_subgroup_sparse/docs/figures/rmsecomp.pdf",width = 6,height = 5)
par(mfrow = c(1,1))
boxplot(cbind(FDA = result1[,2],IND = result2[,2]))
dev.off()

pdf("Google Drive File Stream/My Drive/Research/FDA_subgroup_sparse/docs/figures/numgroupscomp.pdf",width = 6,height = 5)
par(mfrow = c(1,1))
boxplot(cbind(FDA = result1[,5],IND = result2[,5]))
dev.off()








load("~/Research/FDA_subgroup_sparse/result/resultcomp2.rdata")

matcomp2 = do.call("rbind",resultcomp1)
result21 = matcomp2[seq(1,to = 200,by = 2),]
result22 = matcomp2[seq(2,to = 200,by = 2),]

hist(result21[,1])
hist(result22[,1])

which(result22[,1]==0)

boxplot(cbind(result21[,1],result22[,1])[-c(which(result22[,1]==0)),])


sum(result21[,1]>=0.95)
sum(result22[,1]>=0.95)
