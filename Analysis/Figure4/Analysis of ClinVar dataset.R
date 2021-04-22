library(ggpubr)
library(dplyr)
library(ggalt)
library(reshape2)
setwd("D:/m6A-QTL/ClinVar")


###############################################
############ Clinical Significance ############
###############################################
ClinVar <- read.csv("ClinVar_filtered.csv")
colnames(ClinVar)[3] <- "GENE"
ClinVar <- ClinVar[,c(1:7)]

Group <- read.csv("HumanGroup.csv")
StopCodon <- read.csv("HumanStopCodon.csv")
Group <- merge(Group, StopCodon, by=c("chrom","name","strand","txStart","txEnd","cdsStart","cdsEnd","LastExonEnd","LastExonEnd"))
Group <- subset(Group, Group==2)
colnames(Group)[2] <- "GENE"
Group <- Group[,c(1:3)]

ClinVar <- merge(ClinVar, Group, by=c("chrom","GENE","strand"))
ClinVar <- ClinVar %>% distinct(SNPID,GENE, .keep_all = TRUE)


Data <- read.csv("m6AClinVarProbChange_Total.csv")
Data <- merge(Data, ClinVar, by=c("chrom","GENE","strand","SNPID"))
Data <- Data %>% distinct(SNPID,Start,End,GENE, .keep_all = TRUE)



Data$Distance <- Data$Start - Data$POS + 1
Positive <- subset(Data, strand=="+")
Negative <- subset(Data, strand=="-")
Positive$Distance <- -(Positive$Distance)
Data <- rbind(Positive, Negative)


Data$Dvalue <- Data$ProbALT - Data$ProbREF
Positive <- subset(Data, strand=="+")
Negative <- subset(Data, strand=="-")
Negative$REF1[Negative$REF == "A"] <-"T"
Negative$REF1[Negative$REF == "G"] <-"C"
Negative$REF1[Negative$REF == "T"] <-"A"
Negative$REF1[Negative$REF == "C"] <-"G"
Negative$REF <- Negative$REF1
Negative$ALT1[Negative$ALT == "A"] <-"T"
Negative$ALT1[Negative$ALT == "G"] <-"C"
Negative$ALT1[Negative$ALT == "T"] <-"A"
Negative$ALT1[Negative$ALT == "C"] <-"G"
Negative$ALT <- Negative$ALT1

Negative <- Negative[,-c(15,16)]
Data <- rbind(Positive, Negative)

ClinicalSignificance <- read.csv("ClinVar_ClinicalSignificance.csv")
Data <- merge(Data, ClinicalSignificance, by=c("chrom","GENE","strand","SNPID","POS"))
Data <- subset(Data, Distance>= -500&Distance<= 500)
Data <- subset(Data, ClinicalSignificance=="Benign"|ClinicalSignificance=="Likely benign"|ClinicalSignificance=="Likely pathogenic"|ClinicalSignificance=="Pathogenic"|ClinicalSignificance=="Uncertain significance" )

dat <- subset(Data, Distance>= -500&Distance<= 500)
dat$Mark[dat$Dvalue >= 0.1 | dat$Dvalue <= -0.1] <-"Sig"
dat$Mark[dat$Dvalue < 0.1&dat$Dvalue > -0.1] <-"Non"

dat$Marker[dat$ClinicalSignificance=="Benign" | dat$ClinicalSignificance =="Likely benign"] <- "Benign"
dat$Marker[dat$ClinicalSignificance=="Likely pathogenic" | dat$ClinicalSignificance =="Pathogenic"] <- "Patho"
dat$Marker[dat$ClinicalSignificance=="Uncertain significance"] <- "VUS"

dat$Marker <- factor(dat$Marker, levels=c("VUS","Benign","Patho"))
ggplot(dat, aes(x = Marker, y = Dvalue)) + geom_jitter(position=position_jitter(0.2), size = 0.1) + ylim(-1,1)


Sig <- subset(dat, Mark=="Sig")
VUS <- subset(Sig, Marker=="VUS")
Patho <- subset(Sig, Marker=="Patho")
Benign <- subset(Sig, Marker=="Benign")


ggplot(Sig, aes(x = ProbREF, y = ProbALT)) + geom_point()+ xlim(0,1)  + ylim(0,1) + 
  geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
  stat_cor(method = 'spearman', aes(x =ProbREF, y =ProbALT))


Table <- as.data.frame(table(Sig$Distance))
names(Table) <- c("Distance","Count")
merge <- data.frame(Distance=c(-500:500))
merge <- merge(merge, Table, by="Distance", all.x=TRUE)
merge[is.na(merge)] <- 0
ggplot(data=merge, aes(x=Distance, y=Count, group=1)) + geom_line(linetype="blank") + geom_point()
ggplot()+geom_point(data=Sig, aes(x=Distance, y=Dvalue, group=Distance)) + ylim(-1,1) + xlim(-500,500)



###############################################
############## Significant Site ###############
###############################################
ClinVar <- read.csv("ClinVar_filtered.csv")
colnames(ClinVar)[3] <- "GENE"
ClinVar <- ClinVar[,c(1:7)]

Group <- read.csv("HumanGroup.csv")
StopCodon <- read.csv("HumanStopCodon.csv")
Group <- merge(Group, StopCodon, by=c("chrom","name","strand","txStart","txEnd","cdsStart","cdsEnd","LastExonEnd","LastExonEnd"))
Group <- subset(Group, Group==2)
colnames(Group)[2] <- "GENE"
Group <- Group[,c(1:3)]

ClinVar <- merge(ClinVar, Group, by=c("chrom","GENE","strand"))
ClinVar <- ClinVar %>% distinct(SNPID,GENE, .keep_all = TRUE)

Data <- read.csv("m6AClinVarProbChange_Total.csv")
Data <- merge(Data, ClinVar, by=c("chrom","GENE","strand","SNPID"))
Data <- Data %>% distinct(SNPID,Start,End,GENE, .keep_all = TRUE)

Data$Distance <- Data$Start - Data$POS + 1
Positive <- subset(Data, strand=="+")
Negative <- subset(Data, strand=="-")
Positive$Distance <- -(Positive$Distance)
Data <- rbind(Positive, Negative)

Data$Dvalue <- Data$ProbALT - Data$ProbREF
Positive <- subset(Data, strand=="+")
Negative <- subset(Data, strand=="-")
Negative$REF1[Negative$REF == "A"] <-"T"
Negative$REF1[Negative$REF == "G"] <-"C"
Negative$REF1[Negative$REF == "T"] <-"A"
Negative$REF1[Negative$REF == "C"] <-"G"
Negative$REF <- Negative$REF1
Negative$ALT1[Negative$ALT == "A"] <-"T"
Negative$ALT1[Negative$ALT == "G"] <-"C"
Negative$ALT1[Negative$ALT == "T"] <-"A"
Negative$ALT1[Negative$ALT == "C"] <-"G"
Negative$ALT <- Negative$ALT1

Negative <- Negative[,-c(15,16)]
Data <- rbind(Positive, Negative)

ClinicalSignificance <- read.csv("ClinVar_ClinicalSignificance.csv")
Data <- merge(Data, ClinicalSignificance, by=c("chrom","GENE","strand","SNPID","POS"))
Data <- subset(Data, Distance>= -500&Distance<= 500)
Data <- subset(Data, ClinicalSignificance=="Benign"|ClinicalSignificance=="Likely benign"|ClinicalSignificance=="Likely pathogenic"|ClinicalSignificance=="Pathogenic"|ClinicalSignificance=="Uncertain significance" )



dat <- subset(Data, Distance>= -500&Distance<= 500)
dat$Mark[dat$Dvalue >= 0.1 | dat$Dvalue <= -0.1] <-"Sig"
dat$Mark[dat$Dvalue < 0.1&dat$Dvalue > -0.1] <-"Non"

dat$Marker[dat$ClinicalSignificance=="Benign" | dat$ClinicalSignificance =="Likely benign"] <- "Benign"
dat$Marker[dat$ClinicalSignificance=="Likely pathogenic" | dat$ClinicalSignificance =="Pathogenic"] <- "Patho"
dat$Marker[dat$ClinicalSignificance=="Uncertain significance"] <- "VUS"

Sig <- subset(dat, Mark=="Sig")
Non <- subset(dat, Mark=="Non")
VUS <- subset(Sig, Marker=="VUS")
Patho <- subset(Sig, Marker=="Patho")
Benign <- subset(Sig, Marker=="Benign")
write.csv(Patho, "Pathogenic_SNVs.csv", row.names = F)


# Statistic analysis
Sig <- Sig %>% distinct(SNPID, .keep_all = TRUE)
Non <- Non %>% distinct(SNPID, .keep_all = TRUE)


High <- table(Sig$Marker)
Low <- table(Non$Marker)
Table <- cbind(Low, High)
fisher.test(Table, simulate.p.value=T)


Table <- as.data.frame(Table)
write.csv(Table, "Comparision.csv")

Table <- read.csv("Comparision.csv")
colnames(Table)[1] <- "Codon"
Table$Codon <- factor(Table$Codon, levels=c("VUS","Patho","Benign"))
ggplot(data=Table, aes(x=Codon, y=OR, color=Codon)) + geom_bar(stat="identity", fill="white") + ylim(-0.5,0.5)


# Donut plot
High <- as.data.frame(High)
names(High) <- c("Codon","Count")
High$Fraction = High$Count / sum(High$Count)
High$ymax = cumsum(High$Fraction)
High$ymin <- c(0, head(High$ymax, n=-1))
High$LabelPosition <- (High$ymax + High$ymin) / 2
High$label <- paste0(High$Codon, "\n value: ", High$Count)
ggplot(High, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Codon)) +
  geom_rect() +
  geom_label( x=3.5, aes(y=LabelPosition, label=label), size=4) +
  scale_fill_brewer(palette=4) +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none")


Low <- as.data.frame(Low)
names(Low) <- c("Codon","Count")
Low$Fraction = Low$Count / sum(Low$Count)
Low$ymax = cumsum(Low$Fraction)
Low$ymin <- c(0, head(Low$ymax, n=-1))
Low$LabelPosition <- (Low$ymax + Low$ymin) / 2
Low$label <- paste0(Low$Codon, "\n value: ", Low$Count)
ggplot(Low, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Codon)) +
  geom_rect() +
  geom_label( x=3.5, aes(y=LabelPosition, label=label), size=4) +
  scale_fill_brewer(palette=4) +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none")


