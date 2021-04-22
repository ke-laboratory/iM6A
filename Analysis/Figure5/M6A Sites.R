rm(list=ls())
options(stringsAsFactors = F)
library(dplyr)

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
LastExon <- LastExon[c(1,2,3,10,11,12)]
write.csv(LastExon, "M6A_CodonSwap.csv", row.names=F)
