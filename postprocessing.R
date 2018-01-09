#postprocessing on U
load("UN_full.RData")
Time = 32
N = 97
#construct eigenvector-ized U's
Unew = list()
Dnew = list()
R = 2
for (tp in 1:32) {
	Unew[[tp]] = matrix(NA, nrow = nrow(UN$U[[tp]]), ncol = 2*500)
	Dnew[[tp]] = matrix(NA, nrow = 2, ncol = 500)
	rownames(UN$U[[tp]])[which(rownames(UN$U[[tp]]) %in% "GFR")] = "GMY"
	rownames(Unew[[tp]]) = rownames(UN$U[[tp]])
	for (iter in 1:500){
		if (sum(is.na(UN$UPS[[tp]][,1]))>0){
		Uold = UN$UPS[[tp]][-which(is.na(UN$UPS[[tp]][,1])),(2*iter-1):(2*iter)] 
		} else {
			Uold = UN$UPS[[tp]][,(2*iter-1):(2*iter)] 
		}
		UDU = Uold %*% diag(UN$DPS[[tp]][iter,]) %*% t(Uold)
		eULU = eigen(UDU)
		Unew[[tp]][, (2*iter-1):(2*iter)] = eULU$vec[, seq(1, R, length = R), drop = FALSE]	
		Dnew[[tp]][, iter] =  eULU$val[which(rank(-abs(eULU$val), ties.method = "first") <= R)]
		}
}


#1. Do a procrustes for the U's (from eigendecomposition of UDU) in each iteration and get their interval. 
#fix USA and JPN's position but not fix those to be same
library(MCMCpack)
Upost1 <- lapply(1:32, function(tp) Unew[[tp]][,1])
Upost2 <- lapply(1:32, function(tp) Unew[[tp]][,2]
)
for (tp in 1:32) {
	Xstar = matrix(0, nrow = nrow(UN$U[[tp]]), ncol = 2)
	rownames(Xstar) =rownames(UN$U[[tp]])
	Xstar = Unew[[tp]][,1:2] %*% diag(sqrt(Dnew[[tp]][,1]))
	for (iter in 2:500){
	Utrans= procrustes(Unew[[tp]][,(2*iter-1):(2*iter)], Xstar)$X.new %*% diag(sqrt(Dnew[[tp]][,iter]))
	Upost1[[tp]] = cbind(Upost1[[tp]], Utrans[,1])
	Upost2[[tp]] = cbind(Upost2[[tp]], Utrans[,2])
	}
}

yearsummary = list()
for (tp in 1:32) {
	yearsummary[[tp]] = cbind(summary(mcmc(t(Upost1[[tp]])))[[2]][,c(1,3, 5)], summary(mcmc(t(Upost2[[tp]])))[[2]][,c(1,3, 5)])
	rownames(yearsummary[[tp]]) = rownames(UN$U[[tp]])
}
library(plotrix)
par(mfrow = c(4,8), oma = c(3,3,3,3), mar = c(2,1,1,1))
for (tp in 1:32) {
plot(yearsummary[[tp]][,2], yearsummary[[tp]][,5], pch = 16, col = 1:97)
draw.ellipse(yearsummary[[tp]][,2],yearsummary[[tp]][,5], a = (yearsummary[[tp]][,3]-yearsummary[[tp]][,1])/2, b = (yearsummary[[tp]][,6]-yearsummary[[tp]][,4])/2, lty = 2)	
}
library(ggplot2)
library(FastGP)
library(mvtnorm)
library(fields)
library(reshape)
library(MCMCpack)
library(expm)
library(igraph)
library(DLFM2)
library(coda)
library(ggplot2)
library(gridExtra)
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}
number_ticks <- function(n) {function(limits) pretty(limits, n)}
countrynames = ggplotColours(N)
#UD plots_full
setwd('/Users/bomin8319/Desktop/UDU')
years = c(1983:2014)
plots= list()
data2 = list()
sums = list()
for (t in 1:Time) {
  limits = summary(c(unlist(Upost1[[t]]), unlist(Upost2[[t]])))
data2[[t]] = data.frame(Country = c(sapply(rownames(Upost1[[t]]), function(x){rep(x, 500)})), r1 = c(sapply(1:nrow(Upost1[[t]]), function(x){Upost1[[t]][x,]})), r2 = c(sapply(1:nrow(Upost2[[t]]), function(x){Upost2[[t]][x,]})))
p <- ggplot(data2[[t]], aes(x = r1, y = r2, colour = Country, label = Country))
sums[[t]] <- data.frame(r1 = yearsummary[[t]][,2], r2 = yearsummary[[t]][,5], Country = rownames(yearsummary[[t]]))
plots[[t]] = p+ stat_ellipse(show.legend = F, type = "norm")+ ggtitle(years[t]) + theme_minimal()+ theme(plot.title = element_text(hjust = 0.5)) + scale_x_continuous(limits =c(limits[[1]]-0.05, limits[[6]]+0.05), breaks=number_ticks(5)) + scale_y_continuous(limits = c(limits[[1]]-0.05, limits[[6]]+0.05), breaks=number_ticks(5))+geom_text(data = sums[[t]], size = 3, show.legend = F)+scale_color_manual(values = countrynames)
mname = paste0("plot", years[t], "full.png")
print(plots[[t]])
ggsave(filename = mname, width = 10, height = 10)
}
marrangeGrob(plots[c(4,8,12,16,20,24,28,Time)], nrow = 2, ncol = 4, top = NULL)



#UD plots_reduced
data3 = list()
sums2 = list()
for (t in 1:Time) {
countrynames2 = countrynames[which(sort(unique(data2[[t]]$Country)) %in% c("USA", "CHN", "IND", "ROK","PRK","IRQ","RUS","GRG","UKR","UKG", "FRN", "GMY", "TUR", "JPN", "ISR", "SYR", "LEB", "SUD", "IRN", "AUL", "PAK", "EGY","AFG"))]
    data3[[t]] = data2[[t]][data2[[t]]$Country %in% c("USA", "CHN", "IND", "ROK","PRK","IRQ","RUS","GRG","UKR","UKG", "FRN", "GMY", "TUR", "JPN", "ISR", "SYR", "LEB", "SUD", "IRN", "AUL", "PAK", "EGY","AFG"),]
    sums2[[t]] =  sums[[t]][sums[[t]]$Country %in% c("USA", "CHN", "IND", "ROK","PRK","IRQ","RUS","GRG","UKR","UKG", "FRN", "GMY", "TUR", "JPN", "ISR", "SYR", "LEB", "SUD", "IRN", "AUL", "PAK", "EGY","AFG"),]
}
for (t in 1:Time) {
  limits = summary(c(data3[[t]]$r1, data3[[t]]$r2))
p <- ggplot(data3[[t]], aes(x = r1, y = r2, colour = Country, label = Country))
plots[[t]] = p+ stat_ellipse(show.legend = F, type = "norm")+ ggtitle(years[t]) + theme_minimal()+ theme(plot.title = element_text(hjust = 0.5)) + scale_x_continuous(limits = c(limits[[1]]-0.05, limits[[6]]+0.05), breaks=number_ticks(5)) + scale_y_continuous(limits = c(limits[[1]]-0.05, limits[[6]]+0.05),  breaks=number_ticks(5))+geom_text(data = sums2[[t]], size = 3, show.legend = F)+scale_color_manual(values = countrynames2)
mname = paste0("plot", years[t], "reduced.png")
print(plots[[t]])
ggsave(filename = mname,width = 7.5, height = 7.5)

}

data3 = list()
sums2 = list()
for (t in 1:Time) {
	countrynames2 = countrynames[which(sort(unique(data2[[t]]$Country)) %in% c("USA", "CHN", "IND", "ROK","PRK","IRQ","RUS","GRG","UKR","UKG", "FRN", "GMY", "TUR", "JPN", "ISR", "SYR", "LEB", "SUD", "IRN", "AUL", "PAK", "EGY","AFG"))]
    data3[[t]] = data2[[t]][data2[[t]]$Country %in% c("USA", "CHN", "IND", "ROK","PRK","IRQ","RUS","GRG","UKR","UKG", "FRN", "GMY", "TUR", "JPN", "ISR", "SYR", "LEB", "SUD", "IRN", "AUL", "PAK", "EGY","AFG"),]
    sums2[[t]] =  sums[[t]][sums[[t]]$Country %in% c("USA", "CHN", "IND", "ROK","PRK","IRQ","RUS","GRG","UKR","UKG", "FRN", "GMY", "TUR", "JPN", "ISR", "SYR", "LEB", "SUD", "IRN", "AUL", "PAK", "EGY","AFG"),]
    limits = summary(c(data3[[t]]$r1, data3[[t]]$r2))
  p <- ggplot(data3[[t]], aes(x = r1, y = r2, colour = Country, label = Country))
plots[[t]] = p+ stat_ellipse(show.legend = F, type = "norm")+ ggtitle(years[t]) + theme_minimal()+ theme(plot.title = element_text(hjust = 0.5)) + scale_x_continuous(limits = c(limits[[1]]-0.05, limits[[6]]+0.05),  breaks=number_ticks(5)) + scale_y_continuous(limits = c(limits[[1]]-0.05, limits[[6]]+0.05),  breaks=number_ticks(5))+geom_text(data = sums2[[t]], size = 2, show.legend = F)+scale_color_manual(values = countrynames2)
mname = paste0("plot", years[t], "reduced.png")
print(plots[[t]])
#ggsave(filename = mname)
} 
marrangeGrob(plots[c(4,8,12,16,20,24,28,Time)], nrow = 2, ncol = 4, top = NULL)




#2. #eigenvector's product (be unit length) arccosine transformation which will give degrees
#UDU -> UU^T -> interval for pairwise angles
#NOT DISTANCE
angle <- function(x,y){
  dot.prod <- x%*%y 
  norm.x <- norm(x,type="2")
  norm.y <- norm(y,type="2")
  theta <- acos(dot.prod / (norm.x * norm.y))
  if (is.na(theta)) theta <- 0 
  as.numeric(theta)
}
angles = list()
countrypair = list()
for (tp in 1:32) {
	N = nrow(Unew[[tp]])
	countrypair[[tp]] =  t(combn(rownames(Unew[[tp]]), 2))
	angles[[tp]] = matrix(NA, nrow = N * (N-1) / 2, ncol = 500)
	for (iter in 1:500){
		Unow = Unew[[tp]][,(2*iter-1):(2*iter)]
		for (i in 1:nrow(angles[[tp]])) {
			angles[[tp]][i, iter] = angle(Unow[which(rownames(Unow) %in% countrypair[[tp]][i,1]),], Unow[which(rownames(Unow) %in% countrypair[[tp]][i,2]),])
		}
	}
}
save(angles, file = "angles.RData")
load("/Users/bomin8319/Desktop/UDU/angles.RData")
library(MCMCpack)
countrypair2 = list()
for (tp in 1:32) {
	countrypair2[[tp]] = matrix(NA, nrow = nrow(countrypair[[tp]]), ncol = 4)
 for (i in 1:nrow(countrypair[[tp]])) {
 	summary = summary(mcmc(angles[[tp]][i,]))[[2]][c(1,3, 5)]
 	countrypair2[[tp]][i, ] =  c(summary, (summary[1] <= 1 & summary[3] >= 1)) 
 }
 	countrypair[[tp]] = data.frame(countrypair[[tp]], countrypair2[[tp]])
} 

similar = list()
for (tp in 1:32){
	similar[[tp]] =countrypair[[tp]][which(countrypair[[tp]][,6] == 1), 1:2]	
}
