rm(list=ls())
options(stringsAsFactors = F)
setwd("D:/m6AAI/Figure/Fig6")
library(dplyr)

#########################################################
####### m6A sites for single nucleotide mutation ########
#########################################################

# Mouse
LastExon <- read.csv("iM6A_m3cRAC_LastExon.csv")
LastExon <- subset(LastExon, Label==1)
LastExon <- LastExon %>% distinct(name,Start,End, .keep_all = TRUE)
LastExon <- LastExon[order(LastExon$Probabilty, decreasing = TRUE),]
LastExon <- LastExon %>% distinct(name, .keep_all = TRUE)

Group <- read.csv("Group.csv")
LastExon <- merge(LastExon, Group, by=c("chrom","name","strand"))
LastExon <- subset(LastExon, Group==2)

StopCodon <- read.csv("StopCodon.csv")
LastExon <- merge(LastExon, StopCodon, by=c("chrom","name","strand","txStart","txEnd",
                                            "cdsStart","cdsEnd","LastExonStart","LastExonEnd"))

LastExon <- LastExon[order(LastExon$Probabilty, decreasing = TRUE),]
LastExon <- subset(LastExon, Probabilty>=0.4)
LastExon <- LastExon[c(1,2,3,10,11,12)]
write.csv(LastExon, "mm10_M6A_SNM.csv", row.names=F)

# Human
LastExon <- read.csv("iM6A_humanWhistleRAC_LastExon.csv")
LastExon <- subset(LastExon, Label==1)
LastExon <- LastExon %>% distinct(name,Start,End, .keep_all = TRUE)
LastExon <- LastExon[order(LastExon$Probabilty, decreasing = TRUE),]
LastExon <- LastExon %>% distinct(name, .keep_all = TRUE)

Group <- read.csv("HumanGroup.csv")
LastExon <- merge(LastExon, Group, by=c("chrom","name","strand"))
LastExon <- subset(LastExon, Group==2)

StopCodon <- read.csv("humanStopCodon.csv")
LastExon <- merge(LastExon, StopCodon, by=c("chrom","name","strand","txStart","txEnd",
                                            "cdsStart","cdsEnd","LastExonStart","LastExonEnd"))

LastExon <- LastExon[order(LastExon$Probabilty, decreasing = TRUE),]
LastExon <- subset(LastExon, Probabilty>=0.4)
LastExon <- LastExon[c(1,2,3,10,11,12)]
write.csv(LastExon, "hg19_M6A_SNM.csv", row.names=F)



