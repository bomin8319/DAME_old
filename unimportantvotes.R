library(gdata)
dyadic<-read.table('/Users/bomin8319/Desktop/Old Research/DaveMaggie/dataverse_files/Dyadicdata.tab', sep="\t", header=TRUE)
raw <- read.table('/Users/bomin8319/Desktop/Old Research/DaveMaggie/dataverse_files/RawVotingdata.tab', sep="\t", header=TRUE)
ideal <- read.table('/Users/bomin8319/Desktop/Old Research/DaveMaggie/dataverse_files/Idealpoints.tab', sep="\t", header=TRUE)
description <- read.csv('/Users/bomin8319/Desktop/Old Research/DaveMaggie/dataverse_files/descriptionsnew.csv',header=TRUE)

important <- description[description$importantvote==0,-9:-10]
time<-strptime(important[,7], "%m/%d/%y")
years <-unique(time$year+1900)

years = c(1983:2014)
dyadic <- dyadic[dyadic$year %in% years,]
ideal <- ideal[ideal$Year %in% years,]


countrycode <-unique(ideal[,c(2,3,6)])
countrycode<-countrycode[order(countrycode[,1]),]


rawimportant <-raw[raw$rcid %in% important$rcid,]
rawimportant$vote[which(rawimportant$vote==8)]<-4
rawimportant$vote[which(rawimportant$vote==9)]<-5

countrycode2 <-unique(ideal[,c(2,3,6)])
countrycode2<-countrycode2[order(countrycode2[,1]),]



agreedata <-list()
for (t in 1:length(years)){
  print(t)
  agreedata[[t]]<- list(year=years[t])
  data <- dyadic[dyadic$year==years[t],c(1:3, 10, 12) ]
  data <- data[!is.na(data[,4]),]
  countrycode <- unique(c(data$ccode2, data$ccode1))
  countrycode<-countrycode[order(countrycode)]
  agreedata[[t]]$network <- matrix(NA, nrow=length(countrycode), ncol=length(countrycode))
  agreedata[[t]]$jointvote <- matrix(NA, nrow=length(countrycode), ncol=length(countrycode))
  for (i in 1:nrow(data)){
    country1 <- which(countrycode==data[i,1]); country2 <- which(countrycode==data[i,2])
    agreedata[[t]]$network[country1, country2] <- data[i,4]
    agreedata[[t]]$network[country2, country1] <- data[i,4]
    agreedata[[t]]$jointvote[country1, country2] <- data[i,5]
    agreedata[[t]]$jointvote[country2, country1] <- data[i,5]
  }
  rownames(agreedata[[t]]$network)<-countrycode2[countrycode2[,1] %in% countrycode,2]; colnames(agreedata[[t]]$network)<-countrycode2[countrycode2[,1] %in% countrycode,2]	
  rownames(agreedata[[t]]$jointvote)<-countrycode2[countrycode2[,1] %in% countrycode,2]; colnames(agreedata[[t]]$jointvote)<-countrycode2[countrycode2[,1] %in% countrycode,2]	
}

sapply(1:32, function(t) mean(agreedata[[t]]$network[upper.tri(agreedata[[t]]$network)], na.rm = TRUE))
save(agreedata, file='/Users/bomin8319/Desktop/FA16/agreedata.RData')

