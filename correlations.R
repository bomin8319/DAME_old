Xvectors = list()
for (i in 1:5) {
	Xvectors[[i]] = sapply(1:32, function(t){UNdatafull$X[t,,,i+1][upper.tri(UNdatafull$X[t,,,i+1])]})
}

cor(c(Xvectors[[1]]), c(Xvectors[[1]]), use = "complete.obs")
cor(c(Xvectors[[1]]), c(Xvectors[[2]]), use = "complete.obs")
cor(c(Xvectors[[1]]), c(Xvectors[[3]]), use = "complete.obs")
cor(c(Xvectors[[1]]), c(Xvectors[[4]]), use = "complete.obs")
cor(c(Xvectors[[1]]), c(Xvectors[[5]]), use = "complete.obs")

cor(c(Xvectors[[2]]), c(Xvectors[[2]]), use = "complete.obs")
cor(c(Xvectors[[2]]), c(Xvectors[[3]]), use = "complete.obs")
cor(c(Xvectors[[2]]), c(Xvectors[[4]]), use = "complete.obs")
cor(c(Xvectors[[2]]), c(Xvectors[[5]]), use = "complete.obs")

cor(c(Xvectors[[3]]), c(Xvectors[[3]]), use = "complete.obs")
cor(c(Xvectors[[3]]), c(Xvectors[[4]]), use = "complete.obs")
cor(c(Xvectors[[3]]), c(Xvectors[[5]]), use = "complete.obs")

cor(c(Xvectors[[4]]), c(Xvectors[[4]]), use = "complete.obs")
cor(c(Xvectors[[4]]), c(Xvectors[[5]]), use = "complete.obs")