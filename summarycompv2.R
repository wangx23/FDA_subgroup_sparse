##### output for resultnew_v2 #######
library(ggplot2)
library(gridExtra)
setwd("Research/FDA_subgroup/resultnew_v2/")

resultsim3m10 = read.csv("resultsim3m10.csv",header = FALSE)

round(colMeans(resultsim3m10[,1:4]),2)
round(apply(resultsim3m10[,1:4],2,sd),3)

round(colMeans(resultsim3m10[,5:8]),2)
round(apply(resultsim3m10[,5:8],2,sd),2)
colMeans(resultsim3m10[,5:8]==3)


mean(resultsim3m10[,17])
sd(resultsim3m10[,17])
mean(resultsim3m10[,17]==2)

rmse3m10 = data.frame(model = c(rep("IND",100), rep("FDA",100)), 
                      RMSE = c(as.matrix(resultsim3m10[,c(13,15)])))

plot10 = qplot(x = model, y = RMSE, data = rmse3m10, main = "m=10") + 
  geom_boxplot() + theme_bw() +  theme(plot.title = element_text(hjust = 0.5))





resultsim3m20 = read.csv("resultsim3m20.csv",header = FALSE)

round(colMeans(resultsim3m20[,1:4]),2)
round(apply(resultsim3m20[,1:4],2,sd),3)

round(colMeans(resultsim3m20[,5:8]),2)
round(apply(resultsim3m20[,5:8],2,sd),2)
colMeans(resultsim3m20[,5:8]==3)


mean(resultsim3m20[,17])
sd(resultsim3m20[,17])
mean(resultsim3m20[,17]==2)

rmse3m20 = data.frame(model = c(rep("IND",100), rep("FDA",100)), 
                      RMSE = c(as.matrix(resultsim3m20[,c(13,15)])))

plot20 = qplot(x = model, y = RMSE, data = rmse3m20, main = "m=20") + 
  geom_boxplot() + theme_bw() +  theme(plot.title = element_text(hjust = 0.5))


