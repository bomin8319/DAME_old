#full
library(devtools)
setwd('/Users/bomin8319/Desktop/DAME/pkg/R')
load_all()
load("/Users/bomin8319/Desktop/DAME/UNdatafull.RData")
attach(UNdatafull)
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


#reduced

reduced <-which(dimnames(Y)[[2]] %in% c("USA", "CHN", "IND", "UKG", "FRN", "GFR", "TUR", "JPN", "ISR", "SYR", "LEB", "SUD", "IRN", "AUL", "PAK", "EGY", "AFG", "PRK", "RUS", "GRG", "UKR", "ROK", "IRQ"))
  
Y = Y[, reduced, reduced]
X = X[, reduced, reduced, 1:6]

# not existing countries -> all missing values imputed using model (biased)
avail1 = matrix(1, 32, 23)
avail1[1:8, c(12, 20)] =0 #North and South Korea did not joined UN voting until 1990
avail1[1:10, 21] = 0 #GRG no voted until 1992
avail1[1:9, 22] = 0 #RUS X variables not existed until 1991
avail1[1:9, 23] = 0 #UKR not existed until 1991
avail1[13:21, 9] = 0 #IRQ under sanction

UN = DLFM_MH(Y, X, RE = c("additive", "multiplicative"), R = 2, avail = avail1, burn = 1000, nscan = 5000, odens = 10, plot = FALSE)
save(UN, file = "UN_reduced.RData")
UN2 = DLFM_MH(Y, X, RE = c("additive"), R = 2, avail = avail1, burn = 1000, nscan = 5000, odens = 10, plot = FALSE)
save(UN, file = "UN_reduced2.RData")
UN3 = DLFM_MH(Y, X, RE = c(), R = 2, avail = avail1, burn = 1000, nscan = 5000, odens = 10, plot = FALSE)
save(UN, file = "UN_reduced3.RData")


## side by side 
hi = factor(1983:2014)
hello = list()
hello2 = list()
hello3 = list()
hellonew = list()
mean = list()
countrynames = sort(rownames(UN$U[[32]]))

n2 = 1
for 	(n in 1:23){
	hello[[n]] = matrix(NA, nrow = 0, ncol = 3)
	mean[[n]] = matrix(NA, nrow = 0, ncol = 3)
	for (d in 1:32){
	hello[[n]]  = rbind(hello[[n]], cbind(UN$Degree[[d]][,n2], rep(hi[d], length(UN$Degree[[d]][,n2])), rep("DAME", length(UN$Degree[[d]][,n2]))))
	diag(Y[d,, ])= 0
	Y[d, which(avail1[d,]==0), ] = 0
	Y[d, , which(avail1[d,]==0)] = 0
	mean[[n]] = rbind(mean[[n]], c(rowSums(Y[d,,], na.rm=TRUE)[n2], hi[d], "DAME"))
}
colnames(hello[[n]]) = c("Degree", "Year", "Model")
hello[[n]] = as.data.frame(hello[[n]])
hello[[n]]$Year = factor(sort(as.numeric(hello[[n]]$Year)), labels = c(1983:2014))

	hello2[[n]] = matrix(NA, nrow = 0, ncol = 3)
	for (d in 1:32){
	hello2[[n]]  = rbind(hello2[[n]], cbind(UN2$Degree[[d]][,n2], rep(hi[d], length(UN2$Degree[[d]][,n2])), rep("AE", length(UN2$Degree[[d]][,n2]))))
	}
colnames(hello2[[n]]) = c("Degree", "Year", "Model")
hello2[[n]] = as.data.frame(hello2[[n]])
hello2[[n]]$Year = factor(sort(as.numeric(hello2[[n]]$Year)), labels = c(1983:2014))

hello3[[n]] = matrix(NA, nrow = 0, ncol = 3)
	for (d in 1:32){
	hello3[[n]]  = rbind(hello3[[n]], cbind(UN3$Degree[[d]][,n2], rep(hi[d], length(UN3$Degree[[d]][,n2])), rep("NO", length(UN3$Degree[[d]][,n2]))))
	}
colnames(hello3[[n]]) = c("Degree", "Year", "Model")
hello3[[n]] = as.data.frame(hello3[[n]])
hello3[[n]]$Year = factor(sort(as.numeric(hello3[[n]]$Year)), labels = c(1983:2014))


hellonew[[n]] = as.data.frame(rbind(hello[[n]], hello2[[n]], hello3[[n]]))
colnames(mean[[n]]) =  c("Degree", "Year", "Model")
mean[[n]] = as.data.frame(mean[[n]])
n2 = n2 + 1
}

countryname = (rownames(UN$U[[32]]))
countryname2 = (rownames(UN$U[[32]]))
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}
setwd('/Users/bomin8319/Desktop/results')
countrynames = ggplotColours(3)
p10 = list()
n2 = 1
for (n in 1:23){
	hellonew[[n]]$Degree = as.numeric(as.character(hellonew[[n]]$Degree))
	mean[[n]]$Degree = as.numeric(as.character(mean[[n]]$Degree))
	mean[[n]]$Year = hi
p10[[n]] = ggplot(hellonew[[n]], aes(x = Year, y = Degree, color= Model, fill =Model)) + geom_boxplot(outlier.size = 0.5, position = position_dodge()) + theme_minimal() +scale_fill_manual(values = alpha(countrynames,0.5), name = countryname[n])+scale_color_manual(values = countrynames,name = countryname[n]) + geom_line(data = mean[[n]], color = "grey50", size = 1, group = 1)
#annotate("text", 30, 10, size = 5, label = countryname2[n], colour = countrynames[n])
n2 = n2 + 1

mname = paste0(countryname2[n], n, "reduced.png")
print(p10[[n]])
ggsave(filename = mname, width = 12, height = 6)
}


beta = lapply(1:32, function(t){summary(mcmc(UN$BETA[[t]][1:500,]))[[2]]})
betas = list()
for (i in 1:6) {betas[[i]] = sapply(1:32, function(t){beta[[t]][i,]})}
betacols= ggplotColours(6)
number_ticks <- function(n) {function(limits) pretty(limits, n)}
plots = list()
i= 1
years = c(1983:2014)
data = data.frame(cbind(years,t(betas[[i]])))
 colnames(data)[4] = "beta"
 f = ggplot(data, aes(x = years))
plots[[1]]=f + geom_line(aes(y = beta), colour=betacols[1]) + geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), alpha = 0.1) + scale_x_continuous(breaks=number_ticks(8))+ylab("Intercept") + theme_minimal()

color =betacols[-1]
varname = c("log(distance)", "Polity", "Alliance", "Lower trade-to-GDP ratio", "Common Language")
for (i in 2:6){
years = c(1983:2014)
data = data.frame(cbind(years,t(betas[[i]])))
 colnames(data)[4] = "beta"
 f = ggplot(data, aes(x = years))
 
 plots[[i]] <- f + geom_line(aes(y = beta), colour=color[i-1]) + geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), alpha = 0.1) + ylab(varname[i-1]) + scale_x_continuous(breaks=number_ticks(8)) + geom_hline(yintercept = 0) + theme_minimal()

 }
 marrangeGrob(plots[1:6], nrow = 2, ncol = 3, top = NULL)

#thetaplot
colors = sort(rownames(UN$U[[32]]))
colors[which(colors == "GFR")] = "GMY"
rownames(UN$U[[32]])[22] = "GMY"
	thetanew = t(sapply(1:32, function(t){colMeans(UN$theta[[t]])}))
	orders = sapply(1:23, function(n){which(colors[n]== rownames(UN$U[[32]]))})
	thetanew = thetanew[,orders]
	data3 = data.frame(years = years, theta = thetanew)
	colnames(data3)[-1] = colors
	data3new = melt(data3, id = "years")
	colnames(data3new)[3] = "theta"
	f <- ggplot(data3new, aes(years, theta, colour = variable, label = variable))
	f + geom_line() + scale_x_continuous(breaks=number_ticks(8)) + scale_colour_discrete(name = "countries") + theme_minimal()

#UD plots
Xstar = matrix(0, nrow = 23, ncol = 2)
rownames(Xstar) = sort(c("USA", "CHN", "IND", "ROK","PRK","IRQ","RUS","GRG","UKR","UKG", "FRN", "GMY", "TUR", "JPN", "ISR", "SYR", "LEB", "SUD", "IRN", "AUL", "PAK", "EGY","AFG"))
Xstar[23, ]= c(1,1)
Xstar[15, ] = c(-1,-1)
#Xstar[19, ] = c(0, -1)
##Xstar[20, ] = c(0.5, 0.5)
#Xstar[11, ] = c(1, 0)
#Xstar[16,] =c(0.5, 1)
plots= list()
data2 = list()
colors = sort(rownames(Xstar))
for (t in 1:31) {
UDmat = matrix(NA, 23, 2)
rownames(UDmat) = colors
rownames(UN$U[[t]])[which(rownames(UN$U[[t]]) == "GFR")] = "GMY"
UDmat[which(rownames(UDmat) %in% names(UN$U[[t]][which(rownames(UN$U[[t]]) %in% c("USA", "CHN", "IND", "ROK","PRK","IRQ","RUS","GRG","UKR","UKG", "FRN", "GMY", "TUR", "JPN", "ISR", "SYR", "LEB", "SUD", "IRN", "AUL", "PAK", "EGY","AFG")),1] * UN$D[[t]][1] )),1] =  as.numeric(UN$U[[t]][which(rownames(UN$U[[t]]) %in% c("USA", "CHN", "IND", "ROK","PRK","IRQ","RUS","GRG","UKR","UKG", "FRN", "GMY", "TUR", "JPN", "ISR", "SYR", "LEB", "SUD", "IRN", "AUL", "PAK", "EGY","AFG")),1] * UN$D[[t]][1] )
UDmat[which(rownames(UDmat) %in% names(UN$U[[t]][which(rownames(UN$U[[t]]) %in% c("USA", "CHN", "IND", "ROK","PRK","IRQ","RUS","GRG","UKR","UKG", "FRN", "GMY", "TUR", "JPN", "ISR", "SYR", "LEB", "SUD", "IRN", "AUL", "PAK", "EGY","AFG")),2] * UN$D[[t]][2] )),2] =  UN$U[[t]][which(rownames(UN$U[[t]]) %in% c("USA", "CHN", "IND", "ROK","PRK","IRQ","RUS","GRG","UKR","UKG", "FRN", "GMY", "TUR", "JPN", "ISR", "SYR", "LEB", "SUD", "IRN", "AUL", "PAK", "EGY","AFG")),2] * UN$D[[t]][2] 
UDmat[which(rownames(UDmat) %in% names(UN$U[[t]][which(rownames(UN$U[[t]]) %in% c("USA", "CHN", "IND", "ROK","PRK","IRQ","RUS","GRG","UKR","UKG", "FRN", "GMY", "TUR", "JPN", "ISR", "SYR", "LEB", "SUD", "IRN", "AUL", "PAK", "EGY","AFG")),1] * UN$D[[t]][1])), ] = procrustes(UDmat[which(rownames(UDmat) %in% rownames(UN$U[[t]])),], Xstar[which(rownames(Xstar) %in% rownames(UN$U[[t]])),])$X.new
UDmat = UDmat - UDmat[which(rownames(UDmat) =="USA"),]
data2[[t]] = data.frame(UDmat)
colnames(data2[[t]])[1:2] = c("r1", "r2")
Xstar = UDmat
Xstar[is.na(Xstar)] = rep(0, 2)
}
for (t in 1:32) {
UDmat = matrix(NA, 23, 2)
rownames(UDmat) = colors
rownames(UN$U[[t]])[which(rownames(UN$U[[t]]) == "GFR")] = "GMY"
orders = sapply(1:23, function(n){which(rownames(UDmat)[n]== rownames(UN$U[[32]]))})
UDmat[,1] =  as.numeric(UN$U[[t]][which(rownames(UN$U[[t]]) %in% c("USA", "CHN", "IND", "ROK","PRK","IRQ","RUS","GRG","UKR","UKG", "FRN", "GMY", "TUR", "JPN", "ISR", "SYR", "LEB", "SUD", "IRN", "AUL", "PAK", "EGY","AFG")),1] * UN$D[[t]][1] )[orders]
UDmat[,2] =  as.numeric(UN$U[[t]][which(rownames(UN$U[[t]]) %in% c("USA", "CHN", "IND", "ROK","PRK","IRQ","RUS","GRG","UKR","UKG", "FRN", "GMY", "TUR", "JPN", "ISR", "SYR", "LEB", "SUD", "IRN", "AUL", "PAK", "EGY","AFG")),2] * UN$D[[t]][2] )[orders]
UDmat = procrustes(UDmat[which(rownames(UDmat) %in% rownames(UN$U[[t]])),], Xstar[which(rownames(Xstar) %in% rownames(UN$U[[t]])),])$X.new
UDmat = UDmat - UDmat[which(rownames(UDmat) =="USA"),]
data2[[t]] = data.frame(UDmat)
colnames(data2[[t]])[1:2] = c("r1", "r2")
Xstar = UDmat
Xstar[is.na(Xstar)] = rep(0, 2)
}
rangex = summary(unlist(sapply(1:32, function(t){data2[[t]][,1][!is.na(data2[[t]][,1])]})))
rangey = summary(unlist(sapply(1:32, function(t){data2[[t]][,2][!is.na(data2[[t]][,2])]})))

for (t in 1:32) {
colors.t = colors[which(colors %in% rownames(UN$U[[t]]))]
p <- ggplot(data2[[t]], aes(x = r1, y = r2, colour = colors, label = rownames(data2[[t]])))

plots[[t]] = p+ geom_text(size = 10, show.legend = F)+ ggtitle(years[t]) + theme_minimal()+ theme(plot.title = element_text(hjust = 0.5)) + xlim(rangex[1], rangex[6]) + ylim(rangey[1], rangey[6]) 
mname = paste0("plot", t, "reduced.png")
print(plots[[t]])
ggsave(filename = mname)
}

rangex = summary(unlist(sapply(c(4,8,12,16,20,24,28,32), function(t){data2[[t]][,1][!is.na(data2[[t]][,1])]})))
rangey = summary(unlist(sapply(c(4,8,12,16,20,24,28,32), function(t){data2[[t]][,2][!is.na(data2[[t]][,2])]})))
for (t in c(4,8,12,16,20,24,28,32)) {
colors.t = colors[which(colors %in% rownames(UN$U[[t]]))]
p <- ggplot(data2[[t]], aes(x = r1, y = r2, colour = colors, label = rownames(data2[[t]])))

plots[[t]] = p+ geom_text(size = 5, show.legend = F)+ ggtitle(years[t]) + theme_minimal()+ theme(plot.title = element_text(hjust = 0.5)) + xlim(rangex[1]-0.25, rangex[6]+0.25) + ylim(rangey[1], rangey[6]) 
}
marrangeGrob(plots[c(4,8,12,16,20,24,28,32)], nrow = 2, ncol = 4, top = NULL)

