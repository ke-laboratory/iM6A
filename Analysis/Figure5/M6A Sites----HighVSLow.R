rm(list=ls())
options(stringsAsFactors = F)
setwd("D:/m6AAI/Figure/Fig9/HighVSLow")


## High
TopK <- read.csv("iM6A_m3cRAC_Label.csv")
TopK <- TopK[order(TopK$Probabilty, decreasing = TRUE),]

LastExon <- read.csv("iM6A_m3cRAC_LastExon.csv")
LastExon <- LastExon[order(LastExon$Probabilty, decreasing = TRUE),]
Top100K <- TopK[c(1:100000),]

Overlap <- merge(LastExon, Top100K, by=c("chrom","Start","End","name","Probabilty","strand","Label"))
Overlap <- Overlap[order(Overlap$Probabilty, decreasing = TRUE),]
library(dplyr)
Overlap <- Overlap %>% distinct(name, .keep_all = TRUE)

Group <- read.csv("Group.csv")
Group <- Group[-c(4:10)]
StopCodon <- read.csv("StopCodon.csv")
StopCodon <- StopCodon[-c(4:9)]
Group <- merge(Group, StopCodon, by=c("chrom","name","strand"))
Group <- subset(Group, Group==2) 

Overlap <- merge(Overlap, Group, by=c("chrom","name","strand"))
Overlap <- Overlap[order(Overlap$Probabilty, decreasing = TRUE),]
High <- subset(Overlap, Probabilty>=0.7)
High <- High[-c(7:9)]
write.csv(High, "mm10_High_M6ASites.csv", row.names = F)


## Low
Remain <- subset(TopK, Probabilty<0.01)
LastExon <- read.csv("iM6A_m3cRAC_LastExon.csv")
Remain <- merge(LastExon, Remain, by=c("chrom","Start","End","name","Probabilty","strand","Label"))

x <- Overlap$name
y <- Remain$name
inter <- intersect(x,y)

for (i in inter){
  index <- which(colnames(Remain)=="name")
  index.row <- which(Remain[,index]==i)
  Remain <- Remain[-index.row,]
 }

Remain <- merge(Remain, Group, by=c("chrom","name","strand"))
Remain <- Remain[order(Remain$Probabilty, decreasing = TRUE),]
Remain <- Remain[c(1:6)]
write.csv(Remain, "mm10_Low_M6ASites.csv", row.names = F)


