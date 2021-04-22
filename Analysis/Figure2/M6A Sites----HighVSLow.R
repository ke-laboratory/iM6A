rm(list=ls())
options(stringsAsFactors = F)
library(dplyr)


## High Sites
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


## Low Sites
Remain <- subset(TopK, Probabilty<0.1)
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

Remain <- Remain %>% distinct(name, .keep_all = TRUE)
Remain <- Remain[c(1:6)]
write.csv(Remain, "mm10_Low_M6ASites.csv", row.names = F)

























#########################################################
############# High Prob sites in Last Exon ##############
#########################################################
LastExon <- read.csv("iM6A_m3cRAC_LastExon.csv")
LastExon <- subset(LastExon, Label==1)
LastExon <- LastExon[order(LastExon$Probabilty, decreasing = TRUE),]
LastExon <- LastExon %>% distinct(name, .keep_all = TRUE)
Group <- read.csv("Group.csv")
LastExon <- merge(LastExon, Group, by=c("chrom","name","strand"))
LastExon <- subset(LastExon, Group==2)

StopCodon <- read.csv("StopCodon.csv")
LastExon <- merge(LastExon, StopCodon, by=c("chrom","name","strand","txStart","txEnd",
                                            "cdsStart","cdsEnd","LastExonStart","LastExonEnd"))

LastExon <- LastExon[order(LastExon$Probabilty, decreasing = TRUE),]
LastExon <- LastExon[c(1:1000),]
LastExon <- LastExon[c(1,2,3,10,11,12)]
write.csv(LastExon, "High_M6A_SNM.csv", row.names=F)


#########################################################
############# Medium Prob sites in Last Exon ############
#########################################################
LastExon <- read.csv("iM6A_m3cRAC_LastExon.csv")
LastExon <- subset(LastExon, Label==1)
LastExon <- LastExon[order(LastExon$Probabilty, decreasing = TRUE),]
LastExon <- LastExon %>% distinct(name, .keep_all = TRUE)
Group <- read.csv("Group.csv")
LastExon <- merge(LastExon, Group, by=c("chrom","name","strand"))
LastExon <- subset(LastExon, Group==2)

StopCodon <- read.csv("StopCodon.csv")
LastExon <- merge(LastExon, StopCodon, by=c("chrom","name","strand","txStart","txEnd",
                                            "cdsStart","cdsEnd","LastExonStart","LastExonEnd"))

LastExon <- subset(LastExon, Probabilty<=0.6&Probabilty>=0.4)

LastExon <- LastExon[order(LastExon$Probabilty, decreasing = TRUE),]
LastExon <- LastExon[c(1:1000),]
LastExon <- LastExon[c(1,2,3,10,11,12)]
write.csv(LastExon, "Medium_M6A_SNM.csv", row.names=F)






























#########################################################
############# High Prob sites in Last Exon ##############
#########################################################
LastExon <- read.csv("iM6A_m3cRAC_LastExon.csv")
LastExon <- subset(LastExon, Label==1)
LastExon <- LastExon[order(LastExon$Probabilty, decreasing = TRUE),]
LastExon <- LastExon %>% distinct(name, .keep_all = TRUE)
Group <- read.csv("Group.csv")
LastExon <- merge(LastExon, Group, by=c("chrom","name","strand"))
LastExon <- subset(LastExon, Group==2)

StopCodon <- read.csv("StopCodon.csv")
LastExon <- merge(LastExon, StopCodon, by=c("chrom","name","strand","txStart","txEnd",
                                            "cdsStart","cdsEnd","LastExonStart","LastExonEnd"))

LastExon <- LastExon[order(LastExon$Probabilty, decreasing = TRUE),]
LastExon <- LastExon[c(1:1000),]
LastExon <- LastExon[c(1,2,3,10,11,12)]
write.csv(LastExon, "High_M6A_SNM.csv", row.names=F)


#########################################################
############# Medium Prob sites in Last Exon ############
#########################################################
LastExon <- read.csv("iM6A_m3cRAC_LastExon.csv")
LastExon <- subset(LastExon, Label==1)
LastExon <- LastExon[order(LastExon$Probabilty, decreasing = TRUE),]
LastExon <- LastExon %>% distinct(name, .keep_all = TRUE)
Group <- read.csv("Group.csv")
LastExon <- merge(LastExon, Group, by=c("chrom","name","strand"))
LastExon <- subset(LastExon, Group==2)

StopCodon <- read.csv("StopCodon.csv")
LastExon <- merge(LastExon, StopCodon, by=c("chrom","name","strand","txStart","txEnd",
                                            "cdsStart","cdsEnd","LastExonStart","LastExonEnd"))

LastExon <- subset(LastExon, Probabilty<=0.6&Probabilty>=0.4)

LastExon <- LastExon[order(LastExon$Probabilty, decreasing = TRUE),]
LastExon <- LastExon[c(1:1000),]
LastExon <- LastExon[c(1,2,3,10,11,12)]
write.csv(LastExon, "Medium_M6A_SNM.csv", row.names=F)


#########################################################
############### Low Prob sites in Last Exon #############
#########################################################
LastExon <- read.csv("iM6A_m3cRAC_LastExon.csv")
LastExon <- subset(LastExon, Label==1)
LastExon <- LastExon[order(LastExon$Probabilty, decreasing = TRUE),]
LastExon <- LastExon %>% distinct(name, .keep_all = TRUE)
Group <- read.csv("Group.csv")
LastExon <- merge(LastExon, Group, by=c("chrom","name","strand"))
LastExon <- subset(LastExon, Group==2)

StopCodon <- read.csv("StopCodon.csv")
LastExon <- merge(LastExon, StopCodon, by=c("chrom","name","strand","txStart","txEnd",
                                            "cdsStart","cdsEnd","LastExonStart","LastExonEnd"))

LastExon <- subset(LastExon, Probabilty<=0.2)

LastExon <- LastExon[order(LastExon$Probabilty, decreasing = TRUE),]
LastExon <- LastExon[c(1,2,3,10,11,12)]
write.csv(LastExon, "Low_M6A_SNM.csv", row.names=F)







## Position plot
LastExon <- read.csv("iM6A_m3cRAC_LastExon.csv")
LastExon <- subset(LastExon, Label==1)
LastExon <- LastExon[order(LastExon$Probabilty, decreasing = TRUE),]
LastExon <- subset(LastExon, Probabilty<=0.15)
Group <- read.csv("Group.csv")
LastExon <- merge(LastExon, Group, by=c("chrom","name","strand"))
LastExon <- subset(LastExon, Group==2)

StopCodon <- read.csv("StopCodon.csv")
LastExon <- merge(LastExon, StopCodon, by=c("chrom","name","strand","txStart","txEnd",
                                            "cdsStart","cdsEnd","LastExonStart","LastExonEnd"))

LastExon <- LastExon[order(LastExon$Probabilty, decreasing = TRUE),]
LastExon <- LastExon[c(1,2,3,10,11,12)]
write.csv(LastExon, "Low_M6A_SNM.csv", row.names=F)






## Motif phylop score
# High
LastExon <- read.csv("iM6A_m3cRAC_LastExon.csv")
LastExon <- subset(LastExon, Label==1)
LastExon <- LastExon[order(LastExon$Probabilty, decreasing = TRUE),]
LastExon <- LastExon %>% distinct(name, .keep_all = TRUE)
Group <- read.csv("Group.csv")
LastExon <- merge(LastExon, Group, by=c("chrom","name","strand"))
LastExon <- subset(LastExon, Group==2)

StopCodon <- read.csv("StopCodon.csv")
LastExon <- merge(LastExon, StopCodon, by=c("chrom","name","strand","txStart","txEnd",
                                            "cdsStart","cdsEnd","LastExonStart","LastExonEnd"))

LastExon <- LastExon[order(LastExon$Probabilty, decreasing = TRUE),]
LastExon <- LastExon[c(1:2000),]
LastExon <- LastExon[c(1,2,3,10,11,12)]
write.csv(LastExon, "High_M6A_SNM.csv", row.names=F)

# Low
LastExon <- read.csv("iM6A_m3cRAC_LastExon.csv")
LastExon <- subset(LastExon, Label==1)
LastExon <- LastExon[order(LastExon$Probabilty, decreasing = TRUE),]
LastExon <- subset(LastExon, Probabilty<=0.5)
Group <- read.csv("Group.csv")
LastExon <- merge(LastExon, Group, by=c("chrom","name","strand"))
LastExon <- subset(LastExon, Group==2)

StopCodon <- read.csv("StopCodon.csv")
LastExon <- merge(LastExon, StopCodon, by=c("chrom","name","strand","txStart","txEnd",
                                            "cdsStart","cdsEnd","LastExonStart","LastExonEnd"))

LastExon <- LastExon[order(LastExon$Probabilty, decreasing = TRUE),]
LastExon <- LastExon[c(1,2,3,10,11,12)]
write.csv(LastExon, "Low_M6A_SNM.csv", row.names=F)





