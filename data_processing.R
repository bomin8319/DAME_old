library(gdata)
ideal <- read.table('/Users/bomin8319/Desktop/Old Research/DaveMaggie/dataverse_files/Idealpoints.tab', sep="\t", header=TRUE)
description <- read.csv('/Users/bomin8319/Desktop/Old Research/DaveMaggie/dataverse_files/descriptionsnew.csv',header=TRUE)

important <- description[description$importantvote==1,-9:-10]
time<-strptime(important[,7], "%m/%d/%y")
years <-unique(time$year+1900)

ideal <- ideal[ideal$Year %in% years,]

countrycode <-unique(ideal[,c(2,3,6)])
countrycode<-countrycode[order(countrycode[,1]),]

##########################################################
#first take the subset of countries in the response & covariates

load('/Users/bomin8319/Desktop/Old Research/DaveMaggie/agreedata.RData')
years = 1983:2014
Time =32
countrylist = sapply(1:Time, function(t){rownames(agreedata[[t]]$network)})
countrylist_all = countrylist[[1]]
for (i in 2:32){countrylist_all = union(countrylist_all, countrylist[[i]])}

polity <- read.csv('/Users/bomin8319/Desktop/DLFM/UNfit/polity.csv', header = TRUE, sep = ",")
polity <- polity[, -(4:7)]

polcountries = unique(polity$scode)
countrylist_all = intersect(unique(polity[polity$year %in% years,]$scode), countrylist_all)
#result in 163 countries remaining
countrylist_all2 = countrycode[which(countrycode[,2] %in% countrylist_all),3]


GDP <- read.csv('/Users/bomin8319/Desktop/DAME/GDP.csv')

countrylist_all2 = intersect(countrylist_all2, GDP[,1])
#result in 144 countries remaining
#GDP[-which(GDP[,1] %in% countrylist_all2),1]

countrycode_all = unique(countrycode[which(countrycode[,3] %in% countrylist_all2),])
countrylist_all = countrycode_all[,2]

source('/Users/bomin8319/Downloads/Trade_DOT/DOT_EX_MX/DOTX8005.R')
export <- Matrix.X
for (tp in 1:26) {
	rownames(export[[tp]])[which(rownames(export[[tp]])=="GMY")] = "GFR"
    rownames(export[[tp]])[which(rownames(export[[tp]])=="RUM")] = "ROM"
    rownames(export[[tp]])[which(rownames(export[[tp]])=="YEM")] = "YPR"
}
a = rownames(export[[1]])
for (i in 2:26) { a = union(a, rownames(export[[i]]))}
countrylist_all = intersect(a, countrylist_all)

#result in 138 countries remaining
countrycode_all = countrycode_all[which(countrycode_all[,2] %in% countrylist_all),]


#distance overlaps with all 138 countries
dist <- read.csv('/Users/bomin8319/Downloads/Trade_DOT/DOT_EX_MX/distance.csv', header=TRUE)


alliance <- read.csv('/Users/bomin8319/Downloads/version4.1_csv/alliance_v4.1_by_dyad_yearly.csv')
alliance <- alliance[which(alliance$year %in% years), c(2:5, 18)]
alliancelist <- unique(union(alliance$ccode1, alliance$ccode2))
countrycode_all <- countrycode_all[which(countrycode_all[,1] %in% alliancelist),]
countrycode_all <- countrycode_all[-33,] #deleted GFR twice
#result in 119 countries
countrylist_all <- countrycode_all[,2]
#################################################################
#delete countries with large NA's

#nothing deleted with distance
dist <- dist[,-1]
rownames(dist)<- colnames(dist) 

dist2 <- matrix(NA, length(countrylist_all), length(countrylist_all))
rownames(dist2)<- colnames(dist2) <- countrylist_all
for (i in 1:length(countrylist_all)) {
	for (j in 1:length(countrylist_all)) {
	where1 <- which(rownames(dist) %in% countrylist_all[i])
	where2 <- which(colnames(dist) %in% countrylist_all[j])
	if (length(where1) >0 & length(where2) > 0) {
	dist2[i, j] <- dist[where1, where2]
	}
}
}

#polity delte 2 countries
polity <- polity[polity$year %in% years,]
polity[polity$scode == "GMY",]$scode = "GFR"

politylist <-list()
for (tp in 1:Time){ 
  politylist[[tp]] = matrix(NA, length(countrylist_all), length(countrylist_all))
  diag(politylist[[tp]]) = 0
  rownames(politylist[[tp]]) <- colnames(politylist[[tp]]) <- countrylist_all
  for (i in 1:length(countrylist_all)) {
    for (j in 1:length(countrylist_all)) {
      polityi = unique(polity$polity2[polity$scode %in% countrylist_all[i] & polity$year == years[tp]])
      if (length(polityi) > 0 && is.na(polityi)==TRUE) {polityi = 0}
      polityj = unique(polity$polity2[polity$scode %in% countrylist_all[j] & polity$year == years[tp]])
      if (length(polityj) > 0 && is.na(polityj)==TRUE) {polityj = 0}
        if (length(polityi) + length(polityj) >= 2)
        politylist[[tp]][i,j] = abs(polityi - polityj)
    }
  }
}

missing1 <- sapply(1:32, function(t){which(is.na(politylist[[t]][,1]))})
#delete ETH and YPR since no polity score over 10 years
countrycode_all <- countrycode_all[-which(countrycode_all[,2] %in% c("YPR", "ETH")),] #deleted GFR twice
#result in 117 countries
countrylist_all <- countrycode_all[,2]

#reprocess earlier ones
dist2 <- matrix(NA, length(countrylist_all), length(countrylist_all))
rownames(dist2)<- colnames(dist2) <- countrylist_all
for (i in 1:length(countrylist_all)) {
	for (j in 1:length(countrylist_all)) {
	where1 <- which(rownames(dist) %in% countrylist_all[i])
	where2 <- which(colnames(dist) %in% countrylist_all[j])
	if (length(where1) >0 & length(where2) > 0) {
	dist2[i, j] <- dist[where1, where2]
	}
}
}

politylist <-list()
for (tp in 1:Time){ 
  politylist[[tp]] = matrix(NA, length(countrylist_all), length(countrylist_all))
  diag(politylist[[tp]]) = 0
  rownames(politylist[[tp]]) <- colnames(politylist[[tp]]) <- countrylist_all
  for (i in 1:length(countrylist_all)) {
    for (j in 1:length(countrylist_all)) {
      polityi = unique(polity$polity2[polity$scode %in% countrylist_all[i] & polity$year == years[tp]])
      if (length(polityi) > 0 && is.na(polityi)==TRUE) {polityi = 0}
      polityj = unique(polity$polity2[polity$scode %in% countrylist_all[j] & polity$year == years[tp]])
      if (length(polityj) > 0 && is.na(polityj)==TRUE) {polityj = 0}
        if (length(polityi) + length(polityj) >= 2)
        politylist[[tp]][i,j] = abs(polityi - polityj)
    }
  }
}
missing1 <- sapply(1:32, function(t){which(is.na(politylist[[t]][,1]))})

#nothing deleted with alliance (treat none as 0)
allilist <-list()
for (tp in 1:30){ 
  allilist[[tp]] = matrix(0, length(countrylist_all), length(countrylist_all))
  rownames(allilist[[tp]]) <- colnames(allilist[[tp]]) <- countrylist_all
  alliyear = alliance[which(alliance$year == years[tp]),]
  for (i in 1:nrow(alliyear)){
  	allilist[[tp]][which(countrycode_all[,1] == alliyear[i,1]), which(countrycode_all[,1] == alliyear[i,3])] = 1
  	allilist[[tp]][which(countrycode_all[,1] == alliyear[i,3]), which(countrycode_all[,1] == alliyear[i,1])] = 1
  }
   diag(allilist[[tp]]) = 1
   } 	
allilist[[31]] = allilist[[32]] = allilist[[30]]

#process GDP
GDPlist <- list()
for (tp in 1:Time){
  GDPlist[[tp]]<- matrix(NA, nrow=length(countrylist_all), ncol=length(countrylist_all))
  rownames(GDPlist[[tp]]) <- countrylist_all
  colnames(GDPlist[[tp]]) <- countrylist_all
  GDPyeart <- GDP[,c(1,tp+1)]
  for (i in 1:length(countrylist_all)){
  	GDPi <-  as.numeric(GDPyeart[which(GDPyeart[,1] %in% countrycode_all[i,3]),2])
    GDPlist[[tp]][i, ] <- ifelse(length(GDPi) > 0, GDPi, NA)
  }
}
missing2 = sapply(1:32, function(t){which(is.na(GDPlist[[t]][,1])==TRUE)})

#now export
source('/Users/bomin8319/Downloads/Trade_DOT/DOT_EX_MX/DOTX8005.R')
export <- Matrix.X
export <- export[4:26]
for (tp in 1:23) {
	rownames(export[[tp]])[which(rownames(export[[tp]])=="GMY")] = "GFR"
}
exportlist <- list()
for (tp in 1:23){
  exportlist[[tp]] <- matrix(NA, nrow=length(countrylist_all), ncol=length(countrylist_all))
  rownames(exportlist[[tp]]) <- countrylist_all
  colnames(exportlist[[tp]]) <- countrylist_all
  countrylist2 <- which(countrylist_all %in% rownames(export[[tp]]))
  countrylist3 <- sapply(1:length(countrylist2), function(x){which(rownames(export[[tp]]) %in% countrylist_all[countrylist2[x]])})
  exportlist[[tp]][countrylist2, countrylist2] <- as.matrix(export[[tp]][countrylist3,countrylist3], nrow=length(countrylist3), ncol=length(countrylist3))
  diag(exportlist[[tp]]) = 0
  }

export0609 <- read.csv('/Users/bomin8319/Downloads/Trade_DOT/exports_2006-09.csv')
export0609[(export0609 == "n.a.")]<- NA
export0609 <- as.matrix(export0609)
export0609[which(export0609[,1] %in% "Korea, Republic of"),1] = as.character("South Korea")
export0609[which(export0609[,2] %in% "Korea, Republic of"),2] = as.character("South Korea")
export0609[which(export0609[,1] %in% "United States"),1] = as.character("United States of America")
export0609[which(export0609[,2] %in% "United States"),2] = as.character("United States of America")
export0609[which(export0609[,1] %in% "Venezuela, Rep\xfablica Bolivariana de"),1] = as.character("Venezuela")
export0609[which(export0609[,2] %in% "Venezuela, Rep\xfablica Bolivariana de"),2] = as.character("Venezuela")
export0609[which(export0609[,1] %in% "China, P.R.: Mainland"),1] = as.character("China")
export0609[which(export0609[,2] %in% "China, P.R.: Mainland"),2] = as.character("China")
export0609[which(export0609[,1] %in% "Afghanistan, Islamic Republic of"),1] = as.character("Afghanistan")
export0609[which(export0609[,2] %in% "Afghanistan, Islamic Republic of"),2] = as.character("Afghanistan")
export0609[which(export0609[,1] %in% "Germany"),1] = as.character("German Federal Republic")
export0609[which(export0609[,2] %in% "Germany"),2] = as.character("German Federal Republic")
export0609[which(export0609[,1] %in% "Iran, Islamic Republic of"),1] = as.character("Iran")
export0609[which(export0609[,2] %in% "Iran, Islamic Republic of"),2] = as.character("Iran")
export0609[which(export0609[,1] %in% "Russian Federation"),1] = as.character("Russia")
export0609[which(export0609[,2] %in% "Russian Federation"),2] = as.character("Russia")
export0609[which(export0609[,1] %in% "Syrian Arab Republic"),1] = as.character("Syria")
export0609[which(export0609[,2] %in% "Syrian Arab Republic"),2] = as.character("Syria")
export0609[which(export0609[,1] %in% "Slovak Republic"),1] = as.character("Slovakia")
export0609[which(export0609[,2] %in% "Slovak Republic"),2] = as.character("Slovakia")
export0609[which(export0609[,1] %in% "Lao People's Democratic Republic"),1] = as.character("Laos")
export0609[which(export0609[,2] %in% "Lao People's Democratic Republic"),2] = as.character("Laos")
export0609[which(export0609[,1] %in% "Macedonia, FYR"),1] = as.character("Macedonia")
export0609[which(export0609[,2] %in% "Macedonia, FYR"),2] = as.character("Macedonia")
export0609[which(export0609[,1] %in% "Gambia, The"),1] = as.character("Gambia")
export0609[which(export0609[,2] %in% "Gambia, The"),2] = as.character("Gambia")
export0609[which(export0609[,1] %in% "Bahrain, Kingdom of"),1] = as.character("Bahrain")
export0609[which(export0609[,2] %in% "Bahrain, Kingdom of"),2] = as.character("Bahrain")
export0609[which(export0609[,1] %in% "Armenia, Republic of"),1] = as.character("Armenia")
export0609[which(export0609[,2] %in% "Armenia, Republic of"),2] = as.character("Armenia")
export0609[which(export0609[,1] %in% "Kyrgyz Republic"),1] = as.character("Kyrgyzstan")
export0609[which(export0609[,2] %in% "Kyrgyz Republic"),2] = as.character("Kyrgyzstan")
export0609[which(export0609[,1] %in% "Azerbaijan, Republic of"),1] = as.character("Azerbaijan")
export0609[which(export0609[,2] %in% "Azerbaijan, Republic of"),2] = as.character("Azerbaijan")
export0609[which(export0609[,1] %in% "Azerbaijan, Republic of"),1] = as.character("Azerbaijan")
export0609[which(export0609[,2] %in% "Azerbaijan, Republic of"),2] = as.character("Azerbaijan")

export0609 <- export0609[which(export0609[,1] %in% countrycode_all[,3]), -3]
export0609 <- export0609[which(export0609[,2] %in% countrycode_all[,3]), ] 

  
iter = 1
for (i in 24:27){
  exportlist[[i]] <- matrix(NA, nrow=length(countrylist_all), ncol=length(countrylist_all))
  rownames(exportlist[[i]]) <- colnames(exportlist[[i]]) <- countrylist_all
  for (j in 1:length(countrylist_all)){
    for (k in 1:length(countrylist_all)){
    	 exportijt = as.numeric(export0609[which((export0609[,1] %in% countrycode_all[j,3]) & (export0609[,2] %in% countrycode_all[k,3])), iter+2])/10^6
     exportlist[[i]][j,k] <- ifelse(length(exportijt > 0), exportijt, NA)
           }
  }
  diag(exportlist[[i]]) = 0
  iter = iter + 1
}  

setClass("num.with.commas")
setAs("character", "num.with.commas",
function(from) as.numeric(gsub(",","",from)))
export10 <- read.csv('/Users/bomin8319/Downloads/2010export.csv', skip=1, header=TRUE, colClasses=c('character', 'character', 'num.with.commas'))
export10[(export10 == "n.a.")]<- NA
export10 <- as.matrix(export10[,-4])

export10[which(export10[,1] %in% "Afghanistan, I.R. of"),1] = as.character("Afghanistan")
export10[which(export10[,2] %in% "Afghanistan, I.R. of"),2] = as.character("Afghanistan")
export10[which(export10[,1] %in% "Korea, Republic of"),1] = as.character("South Korea")
export10[which(export10[,2] %in% "Korea, Republic of"),2] = as.character("South Korea")
export10[which(export10[,1] %in% "United States"),1] = as.character("United States of America")
export10[which(export10[,2] %in% "United States"),2] = as.character("United States of America")
export10[which(export10[,1] %in% "China,P.R.: Mainland"),1] = as.character("China")
export10[which(export10[,2] %in% "China,P.R.: Mainland"),2] = as.character("China")
export10[which(export10[,1] %in% "Germany"),1] = as.character("German Federal Republic")
export10[which(export10[,2] %in% "Germany"),2] = as.character("German Federal Republic")
export10[which(export10[,1] %in% "Iran, I.R. of"),1] = as.character("Iran")
export10[which(export10[,2] %in% "Iran, I.R. of"),2] = as.character("Iran")
export10[which(export10[,1] %in% "Russian Federation"),1] = as.character("Russia")
export10[which(export10[,2] %in% "Russian Federation"),2] = as.character("Russia")
export10[which(export10[,1] %in% "Syrian Arab Republic"),1] = as.character("Syria")
export10[which(export10[,2] %in% "Syrian Arab Republic"),2] = as.character("Syria")
export10[which(export10[,1] %in% "Bahrain, Kingdom of"),1] = as.character("Bahrain")
export10[which(export10[,2] %in% "Bahrain, Kingdom of"),2] = as.character("Bahrain")
export10 <- export10[which(export10[,1] %in% countrycode_all[,3]), ]
export10 <- export10[which(export10[,2] %in% countrycode_all[,3]), ] 
for (i in 28){
  exportlist[[i]] <- matrix(NA, nrow=length(countrylist_all), ncol=length(countrylist_all))
  rownames(exportlist[[i]]) <- colnames(exportlist[[i]]) <- countrylist_all
  for (j in 1:length(countrylist_all)){
    for (k in 1:length(countrylist_all)){
    	 exportijt = as.numeric(export10[which((export10[,1] %in% countrycode_all[j,3]) & (export10[,2] %in% countrycode_all[k,3])), 3])/10^6
     exportlist[[i]][j,k] <- ifelse(length(exportijt > 0), exportijt, NA)
           }
  }
  diag(exportlist[[i]]) = 0
} 

export11 <- read.csv('/Users/bomin8319/Downloads/2011export.csv', skip=1, header=TRUE, colClasses=c('character', 'character', 'num.with.commas'))
export11[(export11 == "n.a.")]<- NA
export11 <- as.matrix(export11[,-4])
export11[which(export11[,1] %in% "Afghanistan, I.R. of"),1] = as.character("Afghanistan")
export11[which(export11[,2] %in% "Afghanistan, I.R. of"),2] = as.character("Afghanistan")
export11[which(export11[,1] %in% "Korea, Republic of"),1] = as.character("South Korea")
export11[which(export11[,2] %in% "Korea, Republic of"),2] = as.character("South Korea")
export11[which(export11[,1] %in% "United States"),1] = as.character("United States of America")
export11[which(export11[,2] %in% "United States"),2] = as.character("United States of America")
export11[which(export11[,1] %in% "China,P.R.: Mainland"),1] = as.character("China")
export11[which(export11[,2] %in% "China,P.R.: Mainland"),2] = as.character("China")
export11[which(export11[,1] %in% "Germany"),1] = as.character("German Federal Republic")
export11[which(export11[,2] %in% "Germany"),2] = as.character("German Federal Republic")
export11[which(export11[,1] %in% "Iran, I.R. of"),1] = as.character("Iran")
export11[which(export11[,2] %in% "Iran, I.R. of"),2] = as.character("Iran")
export11[which(export11[,1] %in% "Russian Federation"),1] = as.character("Russia")
export11[which(export11[,2] %in% "Russian Federation"),2] = as.character("Russia")
export11[which(export11[,1] %in% "Syrian Arab Republic"),1] = as.character("Syria")
export11[which(export11[,2] %in% "Syrian Arab Republic"),2] = as.character("Syria")
export11[which(export11[,1] %in% "Bahrain, Kingdom of"),1] = as.character("Bahrain")
export11[which(export11[,2] %in% "Bahrain, Kingdom of"),2] = as.character("Bahrain")
for (i in 29){
  exportlist[[i]] <- matrix(NA, nrow=length(countrylist_all), ncol=length(countrylist_all))
  rownames(exportlist[[i]]) <- colnames(exportlist[[i]]) <- countrylist_all
  for (j in 1:length(countrylist_all)){
    for (k in 1:length(countrylist_all)){
    	 exportijt = as.numeric(export11[which((export11[,1] %in% countrycode_all[j,3]) & (export11[,2] %in% countrycode_all[k,3])), 3])/10^6
     exportlist[[i]][j,k] <- ifelse(length(exportijt > 0), exportijt, NA)
           }
  }
  diag(exportlist[[i]]) = 0
} 

export12 <- read.csv('/Users/bomin8319/Downloads/2012export.csv', skip=1, header=TRUE, colClasses=c('character', 'character', 'num.with.commas'))
export12[(export12 == "n.a.")]<- NA
export12 <- as.matrix(export12[,-4])
export12[which(export12[,1] %in% "Afghanistan, I.R. of"),1] = as.character("Afghanistan")
export12[which(export12[,2] %in% "Afghanistan, I.R. of"),2] = as.character("Afghanistan")
export12[which(export12[,1] %in% "Korea, Republic of"),1] = as.character("South Korea")
export12[which(export12[,2] %in% "Korea, Republic of"),2] = as.character("South Korea")
export12[which(export12[,1] %in% "United States"),1] = as.character("United States of America")
export12[which(export12[,2] %in% "United States"),2] = as.character("United States of America")
export12[which(export12[,1] %in% "China,P.R.: Mainland"),1] = as.character("China")
export12[which(export12[,2] %in% "China,P.R.: Mainland"),2] = as.character("China")
export12[which(export12[,1] %in% "Germany"),1] = as.character("German Federal Republic")
export12[which(export12[,2] %in% "Germany"),2] = as.character("German Federal Republic")
export12[which(export12[,1] %in% "Iran, I.R. of"),1] = as.character("Iran")
export12[which(export12[,2] %in% "Iran, I.R. of"),2] = as.character("Iran")
export12[which(export12[,1] %in% "Russian Federation"),1] = as.character("Russia")
export12[which(export12[,2] %in% "Russian Federation"),2] = as.character("Russia")
export12[which(export12[,1] %in% "Syrian Arab Republic"),1] = as.character("Syria")
export12[which(export12[,2] %in% "Syrian Arab Republic"),2] = as.character("Syria")
export12[which(export12[,1] %in% "Bahrain, Kingdom of"),1] = as.character("Bahrain")
export12[which(export12[,2] %in% "Bahrain, Kingdom of"),2] = as.character("Bahrain")
for (i in 30){
  exportlist[[i]] <- matrix(NA, nrow=length(countrylist_all), ncol=length(countrylist_all))
  rownames(exportlist[[i]]) <- colnames(exportlist[[i]]) <- countrylist_all
  for (j in 1:length(countrylist_all)){
    for (k in 1:length(countrylist_all)){
    	 exportijt = as.numeric(export12[which((export12[,1] %in% countrycode_all[j,3]) & (export12[,2] %in% countrycode_all[k,3])), 3])/10^6
     exportlist[[i]][j,k] <- ifelse(length(exportijt > 0), exportijt, NA)
           }
  }
  diag(exportlist[[i]]) = 0
} 

export13 <- read.csv('/Users/bomin8319/Downloads/2013export.csv', skip=1, header=TRUE, colClasses=c('character', 'character', 'num.with.commas'))
export13[(export13 == "n.a.")]<- NA
export13 <- as.matrix(export13[,-4])

export13[which(export13[,1] %in% "Afghanistan, I.R. of"),1] = as.character("Afghanistan")
export13[which(export13[,2] %in% "Afghanistan, I.R. of"),2] = as.character("Afghanistan")
export13[which(export13[,1] %in% "Korea, Republic of"),1] = as.character("South Korea")
export13[which(export13[,2] %in% "Korea, Republic of"),2] = as.character("South Korea")
export13[which(export13[,1] %in% "United States"),1] = as.character("United States of America")
export13[which(export13[,2] %in% "United States"),2] = as.character("United States of America")
export13[which(export13[,1] %in% "China,P.R.: Mainland"),1] = as.character("China")
export13[which(export13[,2] %in% "China,P.R.: Mainland"),2] = as.character("China")
export13[which(export13[,1] %in% "Germany"),1] = as.character("German Federal Republic")
export13[which(export13[,2] %in% "Germany"),2] = as.character("German Federal Republic")
export13[which(export13[,1] %in% "Iran, I.R. of"),1] = as.character("Iran")
export13[which(export13[,2] %in% "Iran, I.R. of"),2] = as.character("Iran")
export13[which(export13[,1] %in% "Russian Federation"),1] = as.character("Russia")
export13[which(export13[,2] %in% "Russian Federation"),2] = as.character("Russia")
export13[which(export13[,1] %in% "Syrian Arab Republic"),1] = as.character("Syria")
export13[which(export13[,2] %in% "Syrian Arab Republic"),2] = as.character("Syria")
export13[which(export13[,1] %in% "Bahrain, Kingdom of"),1] = as.character("Bahrain")
export13[which(export13[,2] %in% "Bahrain, Kingdom of"),2] = as.character("Bahrain")
for (i in 31){
  exportlist[[i]] <- matrix(NA, nrow=length(countrylist_all), ncol=length(countrylist_all))
  rownames(exportlist[[i]]) <- colnames(exportlist[[i]]) <- countrylist_all
  for (j in 1:length(countrylist_all)){
    for (k in 1:length(countrylist_all)){
    	 exportijt = as.numeric(export13[which((export13[,1] %in% countrycode_all[j,3]) & (export13[,2] %in% countrycode_all[k,3])), 3])/10^6
     exportlist[[i]][j,k] <- ifelse(length(exportijt > 0), exportijt, NA)
           }
  }
  diag(exportlist[[i]]) = 0
}

export14 <- read.csv('/Users/bomin8319/Downloads/2014export.csv', skip=1, header=TRUE, colClasses=c('character', 'character', 'num.with.commas'))
export14[(export14 == "n.a.")]<- NA
export14 <- as.matrix(export13[,-4])
export14[which(export14[,1] %in% "Afghanistan, I.R. of"),1] = as.character("Afghanistan")
export14[which(export14[,2] %in% "Afghanistan, I.R. of"),2] = as.character("Afghanistan")
export14[which(export14[,1] %in% "Korea, Republic of"),1] = as.character("South Korea")
export14[which(export14[,2] %in% "Korea, Republic of"),2] = as.character("South Korea")
export14[which(export14[,1] %in% "United States"),1] = as.character("United States of America")
export14[which(export14[,2] %in% "United States"),2] = as.character("United States of America")
export14[which(export14[,1] %in% "China,P.R.: Mainland"),1] = as.character("China")
export14[which(export14[,2] %in% "China,P.R.: Mainland"),2] = as.character("China")
export14[which(export14[,1] %in% "Germany"),1] = as.character("German Federal Republic")
export14[which(export14[,2] %in% "Germany"),2] = as.character("German Federal Republic")
export14[which(export14[,1] %in% "Iran, I.R. of"),1] = as.character("Iran")
export14[which(export14[,2] %in% "Iran, I.R. of"),2] = as.character("Iran")
export14[which(export14[,1] %in% "Russian Federation"),1] = as.character("Russia")
export14[which(export14[,2] %in% "Russian Federation"),2] = as.character("Russia")
export14[which(export14[,1] %in% "Syrian Arab Republic"),1] = as.character("Syria")
export14[which(export14[,2] %in% "Syrian Arab Republic"),2] = as.character("Syria")
export14[which(export14[,1] %in% "Bahrain, Kingdom of"),1] = as.character("Bahrain")
export14[which(export14[,2] %in% "Bahrain, Kingdom of"),2] = as.character("Bahrain")
for (i in 32){
  exportlist[[i]] <- matrix(NA, nrow=length(countrylist_all), ncol=length(countrylist_all))
  rownames(exportlist[[i]]) <- colnames(exportlist[[i]]) <- countrylist_all
  for (j in 1:length(countrylist_all)){
    for (k in 1:length(countrylist_all)){
    	 exportijt = as.numeric(export14[which((export14[,1] %in% countrycode_all[j,3]) & (export14[,2] %in% countrycode_all[k,3])), 3])/10^6
     exportlist[[i]][j,k] <- ifelse(length(exportijt > 0), exportijt, NA)
           }
  }
  diag(exportlist[[i]]) = 0
}

missing3 = sapply(1:32, function(t){which(colSums(exportlist[[t]], na.rm = TRUE)==0)})

#delete BEL, LUX, SAF, CRO, MLD, LAT, BLR, ARM, AZE, TKM, KYR, UZB, KZK since missing >= 10
#     ccode CountryAbb  CountryName
#3854   211        BEL      Belgium
#3855   212        LUX   Luxembourg
#5296   344        CRO      Croatia
#5303   359        MLD      Moldova
#5136   367        LAT       Latvia
#3875   370        BLR      Belarus
#5140   371        ARM      Armenia
#5141   373        AZE   Azerbaijan
#5530   560        SAF South Africa
#5382   701        TKM Turkmenistan
#5384   703        KYR   Kyrgyzstan
#5737   704        UZB   Uzbekistan
#5385   705        KZK   Kazakhstan

#become 104 countries
countrycode_all = countrycode_all[-which(countrycode_all[,2] %in% c("BEL", "LUX", "SAF", "CRO", "MLD", "LAT", "BLR", "ARM", "AZE", "TKM", "KYR", "UZB", "KZK")),]
countrylist_all = countrycode_all[,2]
################
#now reprocess data with these final countries

dist2 <- matrix(NA, length(countrylist_all), length(countrylist_all))
rownames(dist2)<- colnames(dist2) <- countrylist_all
for (i in 1:length(countrylist_all)) {
	for (j in 1:length(countrylist_all)) {
	where1 <- which(rownames(dist) %in% countrylist_all[i])
	where2 <- which(colnames(dist) %in% countrylist_all[j])
	if (length(where1) >0 & length(where2) > 0) {
	dist2[i, j] <- dist[where1, where2]
	}
}
}

politylist <-list()
for (tp in 1:Time){ 
  politylist[[tp]] = matrix(NA, length(countrylist_all), length(countrylist_all))
  diag(politylist[[tp]]) = 0
  rownames(politylist[[tp]]) <- colnames(politylist[[tp]]) <- countrylist_all
  for (i in 1:length(countrylist_all)) {
    for (j in 1:length(countrylist_all)) {
      polityi = unique(polity$polity2[polity$scode %in% countrylist_all[i] & polity$year == years[tp]])
      if (length(polityi) > 0 && is.na(polityi)==TRUE) {polityi = 0}
      polityj = unique(polity$polity2[polity$scode %in% countrylist_all[j] & polity$year == years[tp]])
      if (length(polityj) > 0 && is.na(polityj)==TRUE) {polityj = 0}
        if (length(polityi) + length(polityj) >= 2)
        politylist[[tp]][i,j] = abs(polityi - polityj)
    }
  }
}
missing1 <- sapply(1:32, function(t){which(is.na(politylist[[t]][,1]))})

#nothing deleted with alliance (treat none as 0)
allilist <-list()
for (tp in 1:30){ 
  allilist[[tp]] = matrix(0, length(countrylist_all), length(countrylist_all))
  rownames(allilist[[tp]]) <- colnames(allilist[[tp]]) <- countrylist_all
  alliyear = alliance[which(alliance$year == years[tp]),]
  for (i in 1:nrow(alliyear)){
  	allilist[[tp]][which(countrycode_all[,1] == alliyear[i,1]), which(countrycode_all[,1] == alliyear[i,3])] = 1
  	allilist[[tp]][which(countrycode_all[,1] == alliyear[i,3]), which(countrycode_all[,1] == alliyear[i,1])] = 1
  }
   diag(allilist[[tp]]) = 1
   } 	
allilist[[31]] = allilist[[32]] = allilist[[30]]

#process GDP
GDPlist <- list()
for (tp in 1:Time){
  GDPlist[[tp]]<- matrix(NA, nrow=length(countrylist_all), ncol=length(countrylist_all))
  rownames(GDPlist[[tp]]) <- countrylist_all
  colnames(GDPlist[[tp]]) <- countrylist_all
  GDPyeart <- GDP[,c(1,tp+1)]
  for (i in 1:length(countrylist_all)){
  	GDPi <-  as.numeric(GDPyeart[which(GDPyeart[,1] %in% countrycode_all[i,3]),2])
    GDPlist[[tp]][i, ] <- ifelse(length(GDPi) > 0, GDPi, NA)
  }
}
missing2 = sapply(1:32, function(t){which(is.na(GDPlist[[t]][,1])==TRUE)})

#now export
source('/Users/bomin8319/Downloads/Trade_DOT/DOT_EX_MX/DOTX8005.R')
export <- Matrix.X
export <- export[4:26]
for (tp in 1:23) {
	rownames(export[[tp]])[which(rownames(export[[tp]])=="GMY")] = "GFR"
}
exportlist <- list()
for (tp in 1:23){
  exportlist[[tp]] <- matrix(NA, nrow=length(countrylist_all), ncol=length(countrylist_all))
  rownames(exportlist[[tp]]) <- countrylist_all
  colnames(exportlist[[tp]]) <- countrylist_all
  countrylist2 <- which(countrylist_all %in% rownames(export[[tp]]))
  countrylist3 <- sapply(1:length(countrylist2), function(x){which(rownames(export[[tp]]) %in% countrylist_all[countrylist2[x]])})
  exportlist[[tp]][countrylist2, countrylist2] <- as.matrix(export[[tp]][countrylist3,countrylist3], nrow=length(countrylist3), ncol=length(countrylist3))
  diag(exportlist[[tp]]) = 0
  }


export0609 <- read.csv('/Users/bomin8319/Downloads/Trade_DOT/exports_2006-09.csv')
export0609[(export0609 == "n.a.")]<- NA
export0609 <- as.matrix(export0609)
export0609[which(export0609[,1] %in% "Korea, Republic of"),1] = as.character("South Korea")
export0609[which(export0609[,2] %in% "Korea, Republic of"),2] = as.character("South Korea")
export0609[which(export0609[,1] %in% "United States"),1] = as.character("United States of America")
export0609[which(export0609[,2] %in% "United States"),2] = as.character("United States of America")
export0609[which(export0609[,1] %in% "Venezuela, Rep\xfablica Bolivariana de"),1] = as.character("Venezuela")
export0609[which(export0609[,2] %in% "Venezuela, Rep\xfablica Bolivariana de"),2] = as.character("Venezuela")
export0609[which(export0609[,1] %in% "China, P.R.: Mainland"),1] = as.character("China")
export0609[which(export0609[,2] %in% "China, P.R.: Mainland"),2] = as.character("China")
export0609[which(export0609[,1] %in% "Afghanistan, Islamic Republic of"),1] = as.character("Afghanistan")
export0609[which(export0609[,2] %in% "Afghanistan, Islamic Republic of"),2] = as.character("Afghanistan")
export0609[which(export0609[,1] %in% "Germany"),1] = as.character("German Federal Republic")
export0609[which(export0609[,2] %in% "Germany"),2] = as.character("German Federal Republic")
export0609[which(export0609[,1] %in% "Iran, Islamic Republic of"),1] = as.character("Iran")
export0609[which(export0609[,2] %in% "Iran, Islamic Republic of"),2] = as.character("Iran")
export0609[which(export0609[,1] %in% "Russian Federation"),1] = as.character("Russia")
export0609[which(export0609[,2] %in% "Russian Federation"),2] = as.character("Russia")
export0609[which(export0609[,1] %in% "Syrian Arab Republic"),1] = as.character("Syria")
export0609[which(export0609[,2] %in% "Syrian Arab Republic"),2] = as.character("Syria")
export0609[which(export0609[,1] %in% "Slovak Republic"),1] = as.character("Slovakia")
export0609[which(export0609[,2] %in% "Slovak Republic"),2] = as.character("Slovakia")
export0609[which(export0609[,1] %in% "Lao People's Democratic Republic"),1] = as.character("Laos")
export0609[which(export0609[,2] %in% "Lao People's Democratic Republic"),2] = as.character("Laos")
export0609[which(export0609[,1] %in% "Macedonia, FYR"),1] = as.character("Macedonia")
export0609[which(export0609[,2] %in% "Macedonia, FYR"),2] = as.character("Macedonia")
export0609[which(export0609[,1] %in% "Gambia, The"),1] = as.character("Gambia")
export0609[which(export0609[,2] %in% "Gambia, The"),2] = as.character("Gambia")
export0609[which(export0609[,1] %in% "Bahrain, Kingdom of"),1] = as.character("Bahrain")
export0609[which(export0609[,2] %in% "Bahrain, Kingdom of"),2] = as.character("Bahrain")
export0609[which(export0609[,1] %in% "Armenia, Republic of"),1] = as.character("Armenia")
export0609[which(export0609[,2] %in% "Armenia, Republic of"),2] = as.character("Armenia")
export0609[which(export0609[,1] %in% "Kyrgyz Republic"),1] = as.character("Kyrgyzstan")
export0609[which(export0609[,2] %in% "Kyrgyz Republic"),2] = as.character("Kyrgyzstan")
export0609[which(export0609[,1] %in% "Azerbaijan, Republic of"),1] = as.character("Azerbaijan")
export0609[which(export0609[,2] %in% "Azerbaijan, Republic of"),2] = as.character("Azerbaijan")
export0609[which(export0609[,1] %in% "Azerbaijan, Republic of"),1] = as.character("Azerbaijan")
export0609[which(export0609[,2] %in% "Azerbaijan, Republic of"),2] = as.character("Azerbaijan")

export0609 <- export0609[which(export0609[,1] %in% countrycode_all[,3]), -3]
export0609 <- export0609[which(export0609[,2] %in% countrycode_all[,3]), ] 

  
iter = 1
for (i in 24:27){
  exportlist[[i]] <- matrix(NA, nrow=length(countrylist_all), ncol=length(countrylist_all))
  rownames(exportlist[[i]]) <- colnames(exportlist[[i]]) <- countrylist_all
  for (j in 1:length(countrylist_all)){
    for (k in 1:length(countrylist_all)){
    	 exportijt = as.numeric(export0609[which((export0609[,1] %in% countrycode_all[j,3]) & (export0609[,2] %in% countrycode_all[k,3])), iter+2])/10^6
     exportlist[[i]][j,k] <- ifelse(length(exportijt > 0), exportijt, NA)
           }
  }
  diag(exportlist[[i]]) = 0
  iter = iter + 1
}  

setClass("num.with.commas")
setAs("character", "num.with.commas",
function(from) as.numeric(gsub(",","",from)))
export10 <- read.csv('/Users/bomin8319/Downloads/2010export.csv', skip=1, header=TRUE, colClasses=c('character', 'character', 'num.with.commas'))
export10[(export10 == "n.a.")]<- NA
export10 <- as.matrix(export10[,-4])

export10[which(export10[,1] %in% "Afghanistan, I.R. of"),1] = as.character("Afghanistan")
export10[which(export10[,2] %in% "Afghanistan, I.R. of"),2] = as.character("Afghanistan")
export10[which(export10[,1] %in% "Korea, Republic of"),1] = as.character("South Korea")
export10[which(export10[,2] %in% "Korea, Republic of"),2] = as.character("South Korea")
export10[which(export10[,1] %in% "United States"),1] = as.character("United States of America")
export10[which(export10[,2] %in% "United States"),2] = as.character("United States of America")
export10[which(export10[,1] %in% "China,P.R.: Mainland"),1] = as.character("China")
export10[which(export10[,2] %in% "China,P.R.: Mainland"),2] = as.character("China")
export10[which(export10[,1] %in% "Germany"),1] = as.character("German Federal Republic")
export10[which(export10[,2] %in% "Germany"),2] = as.character("German Federal Republic")
export10[which(export10[,1] %in% "Iran, I.R. of"),1] = as.character("Iran")
export10[which(export10[,2] %in% "Iran, I.R. of"),2] = as.character("Iran")
export10[which(export10[,1] %in% "Russian Federation"),1] = as.character("Russia")
export10[which(export10[,2] %in% "Russian Federation"),2] = as.character("Russia")
export10[which(export10[,1] %in% "Syrian Arab Republic"),1] = as.character("Syria")
export10[which(export10[,2] %in% "Syrian Arab Republic"),2] = as.character("Syria")
export10[which(export10[,1] %in% "Bahrain, Kingdom of"),1] = as.character("Bahrain")
export10[which(export10[,2] %in% "Bahrain, Kingdom of"),2] = as.character("Bahrain")
export10 <- export10[which(export10[,1] %in% countrycode_all[,3]), ]
export10 <- export10[which(export10[,2] %in% countrycode_all[,3]), ] 
for (i in 28){
  exportlist[[i]] <- matrix(NA, nrow=length(countrylist_all), ncol=length(countrylist_all))
  rownames(exportlist[[i]]) <- colnames(exportlist[[i]]) <- countrylist_all
  for (j in 1:length(countrylist_all)){
    for (k in 1:length(countrylist_all)){
    	 exportijt = as.numeric(export10[which((export10[,1] %in% countrycode_all[j,3]) & (export10[,2] %in% countrycode_all[k,3])), 3])/10^6
     exportlist[[i]][j,k] <- ifelse(length(exportijt > 0), exportijt, NA)
           }
  }
  diag(exportlist[[i]]) = 0
} 

export11 <- read.csv('/Users/bomin8319/Downloads/2011export.csv', skip=1, header=TRUE, colClasses=c('character', 'character', 'num.with.commas'))
export11[(export11 == "n.a.")]<- NA
export11 <- as.matrix(export11[,-4])
export11[which(export11[,1] %in% "Afghanistan, I.R. of"),1] = as.character("Afghanistan")
export11[which(export11[,2] %in% "Afghanistan, I.R. of"),2] = as.character("Afghanistan")
export11[which(export11[,1] %in% "Korea, Republic of"),1] = as.character("South Korea")
export11[which(export11[,2] %in% "Korea, Republic of"),2] = as.character("South Korea")
export11[which(export11[,1] %in% "United States"),1] = as.character("United States of America")
export11[which(export11[,2] %in% "United States"),2] = as.character("United States of America")
export11[which(export11[,1] %in% "China,P.R.: Mainland"),1] = as.character("China")
export11[which(export11[,2] %in% "China,P.R.: Mainland"),2] = as.character("China")
export11[which(export11[,1] %in% "Germany"),1] = as.character("German Federal Republic")
export11[which(export11[,2] %in% "Germany"),2] = as.character("German Federal Republic")
export11[which(export11[,1] %in% "Iran, I.R. of"),1] = as.character("Iran")
export11[which(export11[,2] %in% "Iran, I.R. of"),2] = as.character("Iran")
export11[which(export11[,1] %in% "Russian Federation"),1] = as.character("Russia")
export11[which(export11[,2] %in% "Russian Federation"),2] = as.character("Russia")
export11[which(export11[,1] %in% "Syrian Arab Republic"),1] = as.character("Syria")
export11[which(export11[,2] %in% "Syrian Arab Republic"),2] = as.character("Syria")
export11[which(export11[,1] %in% "Bahrain, Kingdom of"),1] = as.character("Bahrain")
export11[which(export11[,2] %in% "Bahrain, Kingdom of"),2] = as.character("Bahrain")
for (i in 29){
  exportlist[[i]] <- matrix(NA, nrow=length(countrylist_all), ncol=length(countrylist_all))
  rownames(exportlist[[i]]) <- colnames(exportlist[[i]]) <- countrylist_all
  for (j in 1:length(countrylist_all)){
    for (k in 1:length(countrylist_all)){
    	 exportijt = as.numeric(export11[which((export11[,1] %in% countrycode_all[j,3]) & (export11[,2] %in% countrycode_all[k,3])), 3])/10^6
     exportlist[[i]][j,k] <- ifelse(length(exportijt > 0), exportijt, NA)
           }
  }
  diag(exportlist[[i]]) = 0
} 

export12 <- read.csv('/Users/bomin8319/Downloads/2012export.csv', skip=1, header=TRUE, colClasses=c('character', 'character', 'num.with.commas'))
export12[(export12 == "n.a.")]<- NA
export12 <- as.matrix(export12[,-4])
export12[which(export12[,1] %in% "Afghanistan, I.R. of"),1] = as.character("Afghanistan")
export12[which(export12[,2] %in% "Afghanistan, I.R. of"),2] = as.character("Afghanistan")
export12[which(export12[,1] %in% "Korea, Republic of"),1] = as.character("South Korea")
export12[which(export12[,2] %in% "Korea, Republic of"),2] = as.character("South Korea")
export12[which(export12[,1] %in% "United States"),1] = as.character("United States of America")
export12[which(export12[,2] %in% "United States"),2] = as.character("United States of America")
export12[which(export12[,1] %in% "China,P.R.: Mainland"),1] = as.character("China")
export12[which(export12[,2] %in% "China,P.R.: Mainland"),2] = as.character("China")
export12[which(export12[,1] %in% "Germany"),1] = as.character("German Federal Republic")
export12[which(export12[,2] %in% "Germany"),2] = as.character("German Federal Republic")
export12[which(export12[,1] %in% "Iran, I.R. of"),1] = as.character("Iran")
export12[which(export12[,2] %in% "Iran, I.R. of"),2] = as.character("Iran")
export12[which(export12[,1] %in% "Russian Federation"),1] = as.character("Russia")
export12[which(export12[,2] %in% "Russian Federation"),2] = as.character("Russia")
export12[which(export12[,1] %in% "Syrian Arab Republic"),1] = as.character("Syria")
export12[which(export12[,2] %in% "Syrian Arab Republic"),2] = as.character("Syria")
export12[which(export12[,1] %in% "Bahrain, Kingdom of"),1] = as.character("Bahrain")
export12[which(export12[,2] %in% "Bahrain, Kingdom of"),2] = as.character("Bahrain")
for (i in 30){
  exportlist[[i]] <- matrix(NA, nrow=length(countrylist_all), ncol=length(countrylist_all))
  rownames(exportlist[[i]]) <- colnames(exportlist[[i]]) <- countrylist_all
  for (j in 1:length(countrylist_all)){
    for (k in 1:length(countrylist_all)){
    	 exportijt = as.numeric(export12[which((export12[,1] %in% countrycode_all[j,3]) & (export12[,2] %in% countrycode_all[k,3])), 3])/10^6
     exportlist[[i]][j,k] <- ifelse(length(exportijt > 0), exportijt, NA)
           }
  }
  diag(exportlist[[i]]) = 0
} 

export13 <- read.csv('/Users/bomin8319/Downloads/2013export.csv', skip=1, header=TRUE, colClasses=c('character', 'character', 'num.with.commas'))
export13[(export13 == "n.a.")]<- NA
export13 <- as.matrix(export13[,-4])

export13[which(export13[,1] %in% "Afghanistan, I.R. of"),1] = as.character("Afghanistan")
export13[which(export13[,2] %in% "Afghanistan, I.R. of"),2] = as.character("Afghanistan")
export13[which(export13[,1] %in% "Korea, Republic of"),1] = as.character("South Korea")
export13[which(export13[,2] %in% "Korea, Republic of"),2] = as.character("South Korea")
export13[which(export13[,1] %in% "United States"),1] = as.character("United States of America")
export13[which(export13[,2] %in% "United States"),2] = as.character("United States of America")
export13[which(export13[,1] %in% "China,P.R.: Mainland"),1] = as.character("China")
export13[which(export13[,2] %in% "China,P.R.: Mainland"),2] = as.character("China")
export13[which(export13[,1] %in% "Germany"),1] = as.character("German Federal Republic")
export13[which(export13[,2] %in% "Germany"),2] = as.character("German Federal Republic")
export13[which(export13[,1] %in% "Iran, I.R. of"),1] = as.character("Iran")
export13[which(export13[,2] %in% "Iran, I.R. of"),2] = as.character("Iran")
export13[which(export13[,1] %in% "Russian Federation"),1] = as.character("Russia")
export13[which(export13[,2] %in% "Russian Federation"),2] = as.character("Russia")
export13[which(export13[,1] %in% "Syrian Arab Republic"),1] = as.character("Syria")
export13[which(export13[,2] %in% "Syrian Arab Republic"),2] = as.character("Syria")
export13[which(export13[,1] %in% "Bahrain, Kingdom of"),1] = as.character("Bahrain")
export13[which(export13[,2] %in% "Bahrain, Kingdom of"),2] = as.character("Bahrain")
for (i in 31){
  exportlist[[i]] <- matrix(NA, nrow=length(countrylist_all), ncol=length(countrylist_all))
  rownames(exportlist[[i]]) <- colnames(exportlist[[i]]) <- countrylist_all
  for (j in 1:length(countrylist_all)){
    for (k in 1:length(countrylist_all)){
    	 exportijt = as.numeric(export13[which((export13[,1] %in% countrycode_all[j,3]) & (export13[,2] %in% countrycode_all[k,3])), 3])/10^6
     exportlist[[i]][j,k] <- ifelse(length(exportijt > 0), exportijt, NA)
           }
  }
  diag(exportlist[[i]]) = 0
}

export14 <- read.csv('/Users/bomin8319/Downloads/2014export.csv', skip=1, header=TRUE, colClasses=c('character', 'character', 'num.with.commas'))
export14[(export14 == "n.a.")]<- NA
export14 <- as.matrix(export13[,-4])
export14[which(export14[,1] %in% "Afghanistan, I.R. of"),1] = as.character("Afghanistan")
export14[which(export14[,2] %in% "Afghanistan, I.R. of"),2] = as.character("Afghanistan")
export14[which(export14[,1] %in% "Korea, Republic of"),1] = as.character("South Korea")
export14[which(export14[,2] %in% "Korea, Republic of"),2] = as.character("South Korea")
export14[which(export14[,1] %in% "United States"),1] = as.character("United States of America")
export14[which(export14[,2] %in% "United States"),2] = as.character("United States of America")
export14[which(export14[,1] %in% "China,P.R.: Mainland"),1] = as.character("China")
export14[which(export14[,2] %in% "China,P.R.: Mainland"),2] = as.character("China")
export14[which(export14[,1] %in% "Germany"),1] = as.character("German Federal Republic")
export14[which(export14[,2] %in% "Germany"),2] = as.character("German Federal Republic")
export14[which(export14[,1] %in% "Iran, I.R. of"),1] = as.character("Iran")
export14[which(export14[,2] %in% "Iran, I.R. of"),2] = as.character("Iran")
export14[which(export14[,1] %in% "Russian Federation"),1] = as.character("Russia")
export14[which(export14[,2] %in% "Russian Federation"),2] = as.character("Russia")
export14[which(export14[,1] %in% "Syrian Arab Republic"),1] = as.character("Syria")
export14[which(export14[,2] %in% "Syrian Arab Republic"),2] = as.character("Syria")
export14[which(export14[,1] %in% "Bahrain, Kingdom of"),1] = as.character("Bahrain")
export14[which(export14[,2] %in% "Bahrain, Kingdom of"),2] = as.character("Bahrain")
for (i in 32){
  exportlist[[i]] <- matrix(NA, nrow=length(countrylist_all), ncol=length(countrylist_all))
  rownames(exportlist[[i]]) <- colnames(exportlist[[i]]) <- countrylist_all
  for (j in 1:length(countrylist_all)){
    for (k in 1:length(countrylist_all)){
    	 exportijt = as.numeric(export14[which((export14[,1] %in% countrycode_all[j,3]) & (export14[,2] %in% countrycode_all[k,3])), 3])/10^6
     exportlist[[i]][j,k] <- ifelse(length(exportijt > 0), exportijt, NA)
           }
  }
  diag(exportlist[[i]]) = 0
}

missing3 = sapply(1:32, function(t){which(colSums(exportlist[[t]], na.rm = TRUE)==0)})

#filled out NA's in exprot manually from other data source -> or just fill in the previous years
for (tp in 2:32) {
	print(sum(is.na(exportlist[[tp]])))
	exportlist[[tp]][which(is.na(exportlist[[tp]]))] = exportlist[[tp-1]][which(is.na(exportlist[[tp]]))]
	print(sum(is.na(exportlist[[tp]])))
}

missing4 = lapply(1:32, function(t){sapply(1:104, function(c){sum(is.na(exportlist[[t]][c,]))>31})})
missing4 = lapply(1:32, function(x){rownames(exportlist[[1]])[missing4[[x]]>0]})
#those miss over 30% of export more than 10 years: EQG, DJI, GNB, MON, CHA, GRG, UKR
countrycode_all = countrycode_all[-which(countrycode_all[,2] %in% c("EQG", "DJI", "GNB", "MON", "CHA", "GRG", "UKR")),]
countrylist_all = countrycode_all[,2]

#become 97 countries
dist2 <- matrix(NA, length(countrylist_all), length(countrylist_all))
rownames(dist2)<- colnames(dist2) <- countrylist_all
for (i in 1:length(countrylist_all)) {
	for (j in 1:length(countrylist_all)) {
	where1 <- which(rownames(dist) %in% countrylist_all[i])
	where2 <- which(colnames(dist) %in% countrylist_all[j])
	if (length(where1) >0 & length(where2) > 0) {
	dist2[i, j] <- dist[where1, where2]
	}
}
}

politylist <-list()
for (tp in 1:Time){ 
  politylist[[tp]] = matrix(NA, length(countrylist_all), length(countrylist_all))
  diag(politylist[[tp]]) = 0
  rownames(politylist[[tp]]) <- colnames(politylist[[tp]]) <- countrylist_all
  for (i in 1:length(countrylist_all)) {
    for (j in 1:length(countrylist_all)) {
      polityi = unique(polity$polity2[polity$scode %in% countrylist_all[i] & polity$year == years[tp]])
      if (length(polityi) > 0 && is.na(polityi)==TRUE) {polityi = 0}
      polityj = unique(polity$polity2[polity$scode %in% countrylist_all[j] & polity$year == years[tp]])
      if (length(polityj) > 0 && is.na(polityj)==TRUE) {polityj = 0}
        if (length(polityi) + length(polityj) >= 2)
        politylist[[tp]][i,j] = abs(polityi - polityj)
    }
  }
}
missing1 <- sapply(1:32, function(t){which(is.na(politylist[[t]][,1]))})

#nothing deleted with alliance (treat none as 0)
allilist <-list()
for (tp in 1:30){ 
  allilist[[tp]] = matrix(0, length(countrylist_all), length(countrylist_all))
  rownames(allilist[[tp]]) <- colnames(allilist[[tp]]) <- countrylist_all
  alliyear = alliance[which(alliance$year == years[tp]),]
  for (i in 1:nrow(alliyear)){
  	allilist[[tp]][which(countrycode_all[,1] == alliyear[i,1]), which(countrycode_all[,1] == alliyear[i,3])] = 1
  	allilist[[tp]][which(countrycode_all[,1] == alliyear[i,3]), which(countrycode_all[,1] == alliyear[i,1])] = 1
  }
   diag(allilist[[tp]]) = 1
   } 	
allilist[[31]] = allilist[[32]] = allilist[[30]]

#process GDP
GDPlist <- list()
for (tp in 1:Time){
  GDPlist[[tp]]<- matrix(NA, nrow=length(countrylist_all), ncol=length(countrylist_all))
  rownames(GDPlist[[tp]]) <- countrylist_all
  colnames(GDPlist[[tp]]) <- countrylist_all
  GDPyeart <- GDP[,c(1,tp+1)]
  for (i in 1:length(countrylist_all)){
  	GDPi <-  as.numeric(GDPyeart[which(GDPyeart[,1] %in% countrycode_all[i,3]),2])
    GDPlist[[tp]][i, ] <- ifelse(length(GDPi) > 0, GDPi, NA)
  }
}
missing2 = sapply(1:32, function(t){which(is.na(GDPlist[[t]][,1])==TRUE)})
#NOTE: imputed Iraq's GDP over 1991 - 1996 using linear equation

#now export
source('/Users/bomin8319/Downloads/Trade_DOT/DOT_EX_MX/DOTX8005.R')
export <- Matrix.X
export <- export[4:26]
for (tp in 1:23) {
	rownames(export[[tp]])[which(rownames(export[[tp]])=="GMY")] = "GFR"
}
exportlist <- list()
for (tp in 1:23){
  exportlist[[tp]] <- matrix(NA, nrow=length(countrylist_all), ncol=length(countrylist_all))
  rownames(exportlist[[tp]]) <- countrylist_all
  colnames(exportlist[[tp]]) <- countrylist_all
  countrylist2 <- which(countrylist_all %in% rownames(export[[tp]]))
  countrylist3 <- sapply(1:length(countrylist2), function(x){which(rownames(export[[tp]]) %in% countrylist_all[countrylist2[x]])})
  exportlist[[tp]][countrylist2, countrylist2] <- as.matrix(export[[tp]][countrylist3,countrylist3], nrow=length(countrylist3), ncol=length(countrylist3))
  diag(exportlist[[tp]]) = 0
  }


export0609 <- read.csv('/Users/bomin8319/Downloads/Trade_DOT/exports_2006-09.csv')
export0609[(export0609 == "n.a.")]<- NA
export0609 <- as.matrix(export0609)
export0609[which(export0609[,1] %in% "Korea, Republic of"),1] = as.character("South Korea")
export0609[which(export0609[,2] %in% "Korea, Republic of"),2] = as.character("South Korea")
export0609[which(export0609[,1] %in% "United States"),1] = as.character("United States of America")
export0609[which(export0609[,2] %in% "United States"),2] = as.character("United States of America")
export0609[which(export0609[,1] %in% "Venezuela, Rep\xfablica Bolivariana de"),1] = as.character("Venezuela")
export0609[which(export0609[,2] %in% "Venezuela, Rep\xfablica Bolivariana de"),2] = as.character("Venezuela")
export0609[which(export0609[,1] %in% "China, P.R.: Mainland"),1] = as.character("China")
export0609[which(export0609[,2] %in% "China, P.R.: Mainland"),2] = as.character("China")
export0609[which(export0609[,1] %in% "Afghanistan, Islamic Republic of"),1] = as.character("Afghanistan")
export0609[which(export0609[,2] %in% "Afghanistan, Islamic Republic of"),2] = as.character("Afghanistan")
export0609[which(export0609[,1] %in% "Germany"),1] = as.character("German Federal Republic")
export0609[which(export0609[,2] %in% "Germany"),2] = as.character("German Federal Republic")
export0609[which(export0609[,1] %in% "Iran, Islamic Republic of"),1] = as.character("Iran")
export0609[which(export0609[,2] %in% "Iran, Islamic Republic of"),2] = as.character("Iran")
export0609[which(export0609[,1] %in% "Russian Federation"),1] = as.character("Russia")
export0609[which(export0609[,2] %in% "Russian Federation"),2] = as.character("Russia")
export0609[which(export0609[,1] %in% "Syrian Arab Republic"),1] = as.character("Syria")
export0609[which(export0609[,2] %in% "Syrian Arab Republic"),2] = as.character("Syria")
export0609[which(export0609[,1] %in% "Slovak Republic"),1] = as.character("Slovakia")
export0609[which(export0609[,2] %in% "Slovak Republic"),2] = as.character("Slovakia")
export0609[which(export0609[,1] %in% "Lao People's Democratic Republic"),1] = as.character("Laos")
export0609[which(export0609[,2] %in% "Lao People's Democratic Republic"),2] = as.character("Laos")
export0609[which(export0609[,1] %in% "Macedonia, FYR"),1] = as.character("Macedonia")
export0609[which(export0609[,2] %in% "Macedonia, FYR"),2] = as.character("Macedonia")
export0609[which(export0609[,1] %in% "Gambia, The"),1] = as.character("Gambia")
export0609[which(export0609[,2] %in% "Gambia, The"),2] = as.character("Gambia")
export0609[which(export0609[,1] %in% "Bahrain, Kingdom of"),1] = as.character("Bahrain")
export0609[which(export0609[,2] %in% "Bahrain, Kingdom of"),2] = as.character("Bahrain")
export0609[which(export0609[,1] %in% "Armenia, Republic of"),1] = as.character("Armenia")
export0609[which(export0609[,2] %in% "Armenia, Republic of"),2] = as.character("Armenia")
export0609[which(export0609[,1] %in% "Kyrgyz Republic"),1] = as.character("Kyrgyzstan")
export0609[which(export0609[,2] %in% "Kyrgyz Republic"),2] = as.character("Kyrgyzstan")
export0609[which(export0609[,1] %in% "Azerbaijan, Republic of"),1] = as.character("Azerbaijan")
export0609[which(export0609[,2] %in% "Azerbaijan, Republic of"),2] = as.character("Azerbaijan")
export0609[which(export0609[,1] %in% "Azerbaijan, Republic of"),1] = as.character("Azerbaijan")
export0609[which(export0609[,2] %in% "Azerbaijan, Republic of"),2] = as.character("Azerbaijan")

export0609 <- export0609[which(export0609[,1] %in% countrycode_all[,3]), -3]
export0609 <- export0609[which(export0609[,2] %in% countrycode_all[,3]), ] 

  
iter = 1
for (i in 24:27){
  exportlist[[i]] <- matrix(NA, nrow=length(countrylist_all), ncol=length(countrylist_all))
  rownames(exportlist[[i]]) <- colnames(exportlist[[i]]) <- countrylist_all
  for (j in 1:length(countrylist_all)){
    for (k in 1:length(countrylist_all)){
    	 exportijt = as.numeric(export0609[which((export0609[,1] %in% countrycode_all[j,3]) & (export0609[,2] %in% countrycode_all[k,3])), iter+2])/10^6
     exportlist[[i]][j,k] <- ifelse(length(exportijt > 0), exportijt, NA)
           }
  }
  diag(exportlist[[i]]) = 0
  iter = iter + 1
}  

setClass("num.with.commas")
setAs("character", "num.with.commas",
function(from) as.numeric(gsub(",","",from)))
export10 <- read.csv('/Users/bomin8319/Downloads/2010export.csv', skip=1, header=TRUE, colClasses=c('character', 'character', 'num.with.commas'))
export10[(export10 == "n.a.")]<- NA
export10 <- as.matrix(export10[,-4])

export10[which(export10[,1] %in% "Afghanistan, I.R. of"),1] = as.character("Afghanistan")
export10[which(export10[,2] %in% "Afghanistan, I.R. of"),2] = as.character("Afghanistan")
export10[which(export10[,1] %in% "Korea, Republic of"),1] = as.character("South Korea")
export10[which(export10[,2] %in% "Korea, Republic of"),2] = as.character("South Korea")
export10[which(export10[,1] %in% "United States"),1] = as.character("United States of America")
export10[which(export10[,2] %in% "United States"),2] = as.character("United States of America")
export10[which(export10[,1] %in% "China,P.R.: Mainland"),1] = as.character("China")
export10[which(export10[,2] %in% "China,P.R.: Mainland"),2] = as.character("China")
export10[which(export10[,1] %in% "Germany"),1] = as.character("German Federal Republic")
export10[which(export10[,2] %in% "Germany"),2] = as.character("German Federal Republic")
export10[which(export10[,1] %in% "Iran, I.R. of"),1] = as.character("Iran")
export10[which(export10[,2] %in% "Iran, I.R. of"),2] = as.character("Iran")
export10[which(export10[,1] %in% "Russian Federation"),1] = as.character("Russia")
export10[which(export10[,2] %in% "Russian Federation"),2] = as.character("Russia")
export10[which(export10[,1] %in% "Syrian Arab Republic"),1] = as.character("Syria")
export10[which(export10[,2] %in% "Syrian Arab Republic"),2] = as.character("Syria")
export10[which(export10[,1] %in% "Bahrain, Kingdom of"),1] = as.character("Bahrain")
export10[which(export10[,2] %in% "Bahrain, Kingdom of"),2] = as.character("Bahrain")
export10 <- export10[which(export10[,1] %in% countrycode_all[,3]), ]
export10 <- export10[which(export10[,2] %in% countrycode_all[,3]), ] 
for (i in 28){
  exportlist[[i]] <- matrix(NA, nrow=length(countrylist_all), ncol=length(countrylist_all))
  rownames(exportlist[[i]]) <- colnames(exportlist[[i]]) <- countrylist_all
  for (j in 1:length(countrylist_all)){
    for (k in 1:length(countrylist_all)){
    	 exportijt = as.numeric(export10[which((export10[,1] %in% countrycode_all[j,3]) & (export10[,2] %in% countrycode_all[k,3])), 3])/10^6
     exportlist[[i]][j,k] <- ifelse(length(exportijt > 0), exportijt, NA)
           }
  }
  diag(exportlist[[i]]) = 0
} 

export11 <- read.csv('/Users/bomin8319/Downloads/2011export.csv', skip=1, header=TRUE, colClasses=c('character', 'character', 'num.with.commas'))
export11[(export11 == "n.a.")]<- NA
export11 <- as.matrix(export11[,-4])
export11[which(export11[,1] %in% "Afghanistan, I.R. of"),1] = as.character("Afghanistan")
export11[which(export11[,2] %in% "Afghanistan, I.R. of"),2] = as.character("Afghanistan")
export11[which(export11[,1] %in% "Korea, Republic of"),1] = as.character("South Korea")
export11[which(export11[,2] %in% "Korea, Republic of"),2] = as.character("South Korea")
export11[which(export11[,1] %in% "United States"),1] = as.character("United States of America")
export11[which(export11[,2] %in% "United States"),2] = as.character("United States of America")
export11[which(export11[,1] %in% "China,P.R.: Mainland"),1] = as.character("China")
export11[which(export11[,2] %in% "China,P.R.: Mainland"),2] = as.character("China")
export11[which(export11[,1] %in% "Germany"),1] = as.character("German Federal Republic")
export11[which(export11[,2] %in% "Germany"),2] = as.character("German Federal Republic")
export11[which(export11[,1] %in% "Iran, I.R. of"),1] = as.character("Iran")
export11[which(export11[,2] %in% "Iran, I.R. of"),2] = as.character("Iran")
export11[which(export11[,1] %in% "Russian Federation"),1] = as.character("Russia")
export11[which(export11[,2] %in% "Russian Federation"),2] = as.character("Russia")
export11[which(export11[,1] %in% "Syrian Arab Republic"),1] = as.character("Syria")
export11[which(export11[,2] %in% "Syrian Arab Republic"),2] = as.character("Syria")
export11[which(export11[,1] %in% "Bahrain, Kingdom of"),1] = as.character("Bahrain")
export11[which(export11[,2] %in% "Bahrain, Kingdom of"),2] = as.character("Bahrain")
for (i in 29){
  exportlist[[i]] <- matrix(NA, nrow=length(countrylist_all), ncol=length(countrylist_all))
  rownames(exportlist[[i]]) <- colnames(exportlist[[i]]) <- countrylist_all
  for (j in 1:length(countrylist_all)){
    for (k in 1:length(countrylist_all)){
    	 exportijt = as.numeric(export11[which((export11[,1] %in% countrycode_all[j,3]) & (export11[,2] %in% countrycode_all[k,3])), 3])/10^6
     exportlist[[i]][j,k] <- ifelse(length(exportijt > 0), exportijt, NA)
           }
  }
  diag(exportlist[[i]]) = 0
} 

export12 <- read.csv('/Users/bomin8319/Downloads/2012export.csv', skip=1, header=TRUE, colClasses=c('character', 'character', 'num.with.commas'))
export12[(export12 == "n.a.")]<- NA
export12 <- as.matrix(export12[,-4])
export12[which(export12[,1] %in% "Afghanistan, I.R. of"),1] = as.character("Afghanistan")
export12[which(export12[,2] %in% "Afghanistan, I.R. of"),2] = as.character("Afghanistan")
export12[which(export12[,1] %in% "Korea, Republic of"),1] = as.character("South Korea")
export12[which(export12[,2] %in% "Korea, Republic of"),2] = as.character("South Korea")
export12[which(export12[,1] %in% "United States"),1] = as.character("United States of America")
export12[which(export12[,2] %in% "United States"),2] = as.character("United States of America")
export12[which(export12[,1] %in% "China,P.R.: Mainland"),1] = as.character("China")
export12[which(export12[,2] %in% "China,P.R.: Mainland"),2] = as.character("China")
export12[which(export12[,1] %in% "Germany"),1] = as.character("German Federal Republic")
export12[which(export12[,2] %in% "Germany"),2] = as.character("German Federal Republic")
export12[which(export12[,1] %in% "Iran, I.R. of"),1] = as.character("Iran")
export12[which(export12[,2] %in% "Iran, I.R. of"),2] = as.character("Iran")
export12[which(export12[,1] %in% "Russian Federation"),1] = as.character("Russia")
export12[which(export12[,2] %in% "Russian Federation"),2] = as.character("Russia")
export12[which(export12[,1] %in% "Syrian Arab Republic"),1] = as.character("Syria")
export12[which(export12[,2] %in% "Syrian Arab Republic"),2] = as.character("Syria")
export12[which(export12[,1] %in% "Bahrain, Kingdom of"),1] = as.character("Bahrain")
export12[which(export12[,2] %in% "Bahrain, Kingdom of"),2] = as.character("Bahrain")
for (i in 30){
  exportlist[[i]] <- matrix(NA, nrow=length(countrylist_all), ncol=length(countrylist_all))
  rownames(exportlist[[i]]) <- colnames(exportlist[[i]]) <- countrylist_all
  for (j in 1:length(countrylist_all)){
    for (k in 1:length(countrylist_all)){
    	 exportijt = as.numeric(export12[which((export12[,1] %in% countrycode_all[j,3]) & (export12[,2] %in% countrycode_all[k,3])), 3])/10^6
     exportlist[[i]][j,k] <- ifelse(length(exportijt > 0), exportijt, NA)
           }
  }
  diag(exportlist[[i]]) = 0
} 

export13 <- read.csv('/Users/bomin8319/Downloads/2013export.csv', skip=1, header=TRUE, colClasses=c('character', 'character', 'num.with.commas'))
export13[(export13 == "n.a.")]<- NA
export13 <- as.matrix(export13[,-4])

export13[which(export13[,1] %in% "Afghanistan, I.R. of"),1] = as.character("Afghanistan")
export13[which(export13[,2] %in% "Afghanistan, I.R. of"),2] = as.character("Afghanistan")
export13[which(export13[,1] %in% "Korea, Republic of"),1] = as.character("South Korea")
export13[which(export13[,2] %in% "Korea, Republic of"),2] = as.character("South Korea")
export13[which(export13[,1] %in% "United States"),1] = as.character("United States of America")
export13[which(export13[,2] %in% "United States"),2] = as.character("United States of America")
export13[which(export13[,1] %in% "China,P.R.: Mainland"),1] = as.character("China")
export13[which(export13[,2] %in% "China,P.R.: Mainland"),2] = as.character("China")
export13[which(export13[,1] %in% "Germany"),1] = as.character("German Federal Republic")
export13[which(export13[,2] %in% "Germany"),2] = as.character("German Federal Republic")
export13[which(export13[,1] %in% "Iran, I.R. of"),1] = as.character("Iran")
export13[which(export13[,2] %in% "Iran, I.R. of"),2] = as.character("Iran")
export13[which(export13[,1] %in% "Russian Federation"),1] = as.character("Russia")
export13[which(export13[,2] %in% "Russian Federation"),2] = as.character("Russia")
export13[which(export13[,1] %in% "Syrian Arab Republic"),1] = as.character("Syria")
export13[which(export13[,2] %in% "Syrian Arab Republic"),2] = as.character("Syria")
export13[which(export13[,1] %in% "Bahrain, Kingdom of"),1] = as.character("Bahrain")
export13[which(export13[,2] %in% "Bahrain, Kingdom of"),2] = as.character("Bahrain")
for (i in 31){
  exportlist[[i]] <- matrix(NA, nrow=length(countrylist_all), ncol=length(countrylist_all))
  rownames(exportlist[[i]]) <- colnames(exportlist[[i]]) <- countrylist_all
  for (j in 1:length(countrylist_all)){
    for (k in 1:length(countrylist_all)){
    	 exportijt = as.numeric(export13[which((export13[,1] %in% countrycode_all[j,3]) & (export13[,2] %in% countrycode_all[k,3])), 3])/10^6
     exportlist[[i]][j,k] <- ifelse(length(exportijt > 0), exportijt, NA)
           }
  }
  diag(exportlist[[i]]) = 0
}

export14 <- read.csv('/Users/bomin8319/Downloads/2014export.csv', skip=1, header=TRUE, colClasses=c('character', 'character', 'num.with.commas'))
export14[(export14 == "n.a.")]<- NA
export14 <- as.matrix(export13[,-4])
export14[which(export14[,1] %in% "Afghanistan, I.R. of"),1] = as.character("Afghanistan")
export14[which(export14[,2] %in% "Afghanistan, I.R. of"),2] = as.character("Afghanistan")
export14[which(export14[,1] %in% "Korea, Republic of"),1] = as.character("South Korea")
export14[which(export14[,2] %in% "Korea, Republic of"),2] = as.character("South Korea")
export14[which(export14[,1] %in% "United States"),1] = as.character("United States of America")
export14[which(export14[,2] %in% "United States"),2] = as.character("United States of America")
export14[which(export14[,1] %in% "China,P.R.: Mainland"),1] = as.character("China")
export14[which(export14[,2] %in% "China,P.R.: Mainland"),2] = as.character("China")
export14[which(export14[,1] %in% "Germany"),1] = as.character("German Federal Republic")
export14[which(export14[,2] %in% "Germany"),2] = as.character("German Federal Republic")
export14[which(export14[,1] %in% "Iran, I.R. of"),1] = as.character("Iran")
export14[which(export14[,2] %in% "Iran, I.R. of"),2] = as.character("Iran")
export14[which(export14[,1] %in% "Russian Federation"),1] = as.character("Russia")
export14[which(export14[,2] %in% "Russian Federation"),2] = as.character("Russia")
export14[which(export14[,1] %in% "Syrian Arab Republic"),1] = as.character("Syria")
export14[which(export14[,2] %in% "Syrian Arab Republic"),2] = as.character("Syria")
export14[which(export14[,1] %in% "Bahrain, Kingdom of"),1] = as.character("Bahrain")
export14[which(export14[,2] %in% "Bahrain, Kingdom of"),2] = as.character("Bahrain")
for (i in 32){
  exportlist[[i]] <- matrix(NA, nrow=length(countrylist_all), ncol=length(countrylist_all))
  rownames(exportlist[[i]]) <- colnames(exportlist[[i]]) <- countrylist_all
  for (j in 1:length(countrylist_all)){
    for (k in 1:length(countrylist_all)){
    	 exportijt = as.numeric(export14[which((export14[,1] %in% countrycode_all[j,3]) & (export14[,2] %in% countrycode_all[k,3])), 3])/10^6
     exportlist[[i]][j,k] <- ifelse(length(exportijt > 0), exportijt, NA)
           }
  }
  diag(exportlist[[i]]) = 0
}

missing3 = sapply(1:32, function(t){which(colSums(exportlist[[t]], na.rm = TRUE)==0)})

#filled out NA's in exprot manually from other data source -> or just fill in the previous years
for (tp in 2:32) {
	print(sum(is.na(exportlist[[tp]])))
	exportlist[[tp]][which(is.na(exportlist[[tp]]))] = exportlist[[tp-1]][which(is.na(exportlist[[tp]]))]
	print(sum(is.na(exportlist[[tp]])))
}

missing4 = lapply(1:32, function(t){sapply(1:97, function(c){sum(is.na(exportlist[[t]][c,]))>31})})
missing4 = lapply(1:32, function(x){rownames(exportlist[[1]])[missing4[[x]]>0]})

tradelist <- list()
for (tp in 1:32) {
    exportlist[[tp]][is.na(exportlist[[tp]])] = 0 #finally, missing export plugged in as zeros
	tradelist[[tp]] = exportlist[[tp]] + t(exportlist[[tp]])
}


load("/Users/bomin8319/Desktop/DAME/clang.RData")
rownames(clang)[36] = "GFR"
colnames(clang)[36] = "GFR"
clanglist <- matrix(NA, length(countrylist_all), length(countrylist_all))
rownames(clanglist) <- colnames(clanglist) <- countrylist_all
 for (i in 1:length(countrylist_all)) {
	for (j in 1:length(countrylist_all)) {
	where1 <- which(rownames(clang) %in% countrylist_all[i])
	where2 <- which(colnames(clang) %in% countrylist_all[j])
	if (length(where1) >0 & length(where2) > 0) {
	clanglist[i, j] <- clang[where1, where2]
	}
}
}
diag(clanglist) = 1
#Haiti, Mali, Benin, Guinea, Burkina Faso, Togo, Gabon,  Central African Republic, Congo, Burundi, Rwanda uses French
clanglist[3,] = clanglist[28,]
clanglist[,3] = clanglist[,28]
clanglist[43, ] = clanglist[28,]
clanglist[,43] = clanglist[,28]
clanglist[45, ] = clanglist[28,]
clanglist[,45] = clanglist[,28]
clanglist[48, ] = clanglist[28,]
clanglist[,48] = clanglist[,28]
clanglist[49, ] = clanglist[28,]
clanglist[,49] = clanglist[,28]
clanglist[52, ] = clanglist[28,]
clanglist[,52] = clanglist[,28]
clanglist[55, ] = clanglist[28,]
clanglist[,55] = clanglist[,28]
clanglist[56, ] = clanglist[28,]
clanglist[,56] = clanglist[,28]
clanglist[57, ] = clanglist[28,]
clanglist[,57] = clanglist[,28]
clanglist[61, ] = clanglist[28,]
clanglist[,61] = clanglist[,28]
clanglist[62, ] = clanglist[28,]
clanglist[,62] = clanglist[,28]
#El Salvador, Nicaragua uses Spanish
clanglist[10, ] = clanglist[29,]
clanglist[,10] = clanglist[,29]
clanglist[11, ] = clanglist[29,]
clanglist[,11] = clanglist[,29]
#Suriname uses Dutch
clanglist[17, ] = clanglist[31,]
clanglist[,17] = clanglist[,31]
#Algania uses Albanian
clanglist[35, ] = rep(0, 97)
clanglist[,35] = rep(0, 97)
clanglist[35, 35] = 1
#Gambia, Sierra Leone, Ghana, Uganda, Tanzania, Rwanda uses English
clanglist[42, ] = clanglist[1,]
clanglist[,42] = clanglist[,1]
clanglist[50, ] = clanglist[1,]
clanglist[,50] = clanglist[,1]
clanglist[51, ] = clanglist[1,]
clanglist[,51] = clanglist[,1]
clanglist[58, ] = clanglist[1,]
clanglist[,58] = clanglist[,1]
clanglist[60, ] = clanglist[1,]
clanglist[,60] = clanglist[,1]
clanglist[62,which(clanglist[,1]>0)] = 1
clanglist[which(clanglist[,1]>0),62] =1
#Mauritania, Algeria, Libya, Iraq uses Arabic
clanglist[46, ] = clanglist[77,]
clanglist[,46] = clanglist[,77]
clanglist[68, ] = clanglist[77,]
clanglist[,68] = clanglist[,77]
clanglist[70, ] = clanglist[77,]
clanglist[,70] = clanglist[,77]
clanglist[74, ] = clanglist[77,]
clanglist[,74] = clanglist[,77]
#Angola,  Mozambique uses Protuguise
clanglist[63,] = clanglist[30,]
clanglist[,63] = clanglist[,30]
clanglist[64,] = clanglist[30,]
clanglist[,64] = clanglist[,30]



Ylist <- list()
for (tp in 1:Time){
  Ylist[[tp]] <- matrix(NA, nrow=length(countrylist_all), ncol=length(countrylist_all))
  rownames(Ylist[[tp]]) <- countrylist_all
  colnames(Ylist[[tp]]) <- countrylist_all
  countrylist2 <- which(countrylist_all %in% rownames(agreedata[[tp]]$network))
  countrylist3 <- sapply(1:length(countrylist2), function(x){which(rownames(agreedata[[tp]]$network) %in% countrylist_all[countrylist2[x]])})
  Ylist[[tp]][countrylist2, countrylist2] <-agreedata[[tp]]$network[countrylist3, countrylist3] 
}
missing5 = sapply(1:32, function(t){which(colSums(Ylist[[t]], na.rm = TRUE)==0)})


Ylist2 <- list()
for (tp in 1:Time){
    Ylist2[[tp]] <- matrix(NA, nrow=length(countrylist_all), ncol=length(countrylist_all))
    rownames(Ylist2[[tp]]) <- countrylist_all
    colnames(Ylist2[[tp]]) <- countrylist_all
    countrylist2 <- which(countrylist_all %in% rownames(agreedata[[tp]]$jointvote))
    countrylist3 <- sapply(1:length(countrylist2), function(x){which(rownames(agreedata[[tp]]$jointvote) %in% countrylist_all[countrylist2[x]])})
    Ylist2[[tp]][countrylist2, countrylist2] <-agreedata[[tp]]$jointvote[countrylist3, countrylist3]
}


Y <- array(NA, dim=c(Time, length(countrylist_all), length(countrylist_all)))
dimnames(Y) <- list(c(1983:2014), countrylist_all, countrylist_all)
for (tp in 1:Time){
  Y[tp,,]<- Ylist[[tp]]
  diag(Y[tp,,]) = 1
}

P = 1 + 7
N = dim(Y)[2]
X <- array(NA, dim=c(Time, N, N, P))
dimnames(X) <- list(c(1983:2014), countrylist_all, countrylist_all, c("intercept", "log(distance)", "polity", "alliance", "min(Trade_ij/GDP_i, Trade_ij/GDP_j)", "language", "GDP", "Export"))

for (tp in 1:Time){
  X[tp,,,1] <- 1
  X[tp,,,2] <- as.matrix(log(dist2), nrow=dim(Y)[2], ncol=dim(Y)[2])
  diag(X[tp,,,2]) <- 0
  X[tp,,,3] <- as.matrix(politylist[[tp]], nrow=dim(Y)[2], ncol=dim(Y)[2]) 
  X[tp,,,4] <- as.matrix(allilist[[tp]], nrow=dim(Y)[2], ncol=dim(Y)[2]) 
  for (i in 1:dim(Y)[2]){
  	for (j in 1:dim(Y)[2]){
  		ratioi = tradelist[[tp]][i,j] / GDPlist[[tp]][i,1]
  		ratioj = tradelist[[tp]][i,j] / GDPlist[[tp]][j,1]
  		X[tp, i, j, 5] = min(ratioi, ratioj)
  	}
  }
  X[tp,,,6] =  as.matrix(clanglist, nrow=dim(Y)[2], ncol=dim(Y)[2])
  X[tp,,,7] <- as.matrix(GDPlist[[tp]], nrow=dim(Y)[2], ncol=dim(Y)[2])
  X[tp,,,8] <- as.matrix(exportlist[[tp]], nrow=dim(Y)[2], ncol=dim(Y)[2]) 
  	}

UNdatafull = list()
UNdatafull$Y = Y
UNdatafull$X = X

save(UNdatafull, file = "UNdatafull2.RData")
