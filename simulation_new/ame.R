library(amen)
ame_sim = list()

for (tp in 1:10) {
	ame_sim[[tp]] = ame(Ys[tp,,], symmetric = TRUE, R=2, nscan = 1000, burn = 100, odens = 5)	
	
}

sapply(1:10, function(x) var(c(ame_sim[[x]]$ULUPM[upper.tri(diag(20))])) / mean(ame_sim[[x]]$VC[,1]))


library(igraph)

par(mfrow = c(2,5))
for (tp in 1:10) {
	net <- graph_from_adjacency_matrix(Ys[tp,,], mode = "undirected", diag = FALSE)
	plot(net)
	
}


sapply(1:5, function(x) sum((UDU[[x]][upper.tri(UDU[[x]])]- ame_sim[[x]]$ULUPM[upper.tri(UDU[[x]])])^2))


sapply(1:10, function(x) sum((UDU[[x]][upper.tri(UDU[[x]])]- M3$UDU[[x]][upper.tri(UDU[[x]])])^2))

sapply(1:10, function(x) cor(UDU[[x]][upper.tri(UDU[[x]])], M3$UDU[[x]][upper.tri(UDU[[x]])]))
sapply(1:10, function(x) cor(UDU[[x]][upper.tri(UDU[[x]])], M1$UDU[[x]][upper.tri(UDU[[x]])]))




sapply(1:10, function(x) sum((Ys[tp,,][upper.tri(UDU[[x]])]- M1$YPM[[x]][upper.tri(UDU[[x]])])^2))
sapply(1:10, function(x) sum((Ys[tp,,][upper.tri(UDU[[x]])]-M3$YPM[[x]][upper.tri(UDU[[x]])])^2))
