rm(list=ls())
options(stringsAsFactors = F)
library(ggpubr)
library(ggplot2)
library(reshape2)
library(dplyr)

###################################################
###################################################
NRACTRR <- read.csv("NRACTRR.csv")
TRACN <- read.csv("TRACN.csv")
TRRACN <- read.csv("TRRACN.csv")
Group <- read.csv("Group.csv")

Merge <- rbind(NRACTRR,TRACN,TRRACN)
Merge <- Merge %>% distinct(name, .keep_all = TRUE)


First <- read.csv("First.csv")
First <- First[-c(5)][-c(2,3)]
colnames(First)[4] <- "First"
Second <- read.csv("Second.csv")
Second <- Second[-c(6)][-c(2,3)]
colnames(Second)[4] <- "Second"
Third <- read.csv("Third.csv")
Third <- Third[-c(6)][-c(2,3)]
colnames(Third)[4] <- "Third"

Score <- merge(First, Second, by=c("chrom","name","strand"))
Score <- merge(Score, Third, by=c("chrom","name","strand"))
Score$Score <- (Score$First+Score$Second+Score$Third)/3


Data <- merge(Merge, Group, by=c("name","chrom","strand","txStart","txEnd",
                                 "cdsStart","cdsEnd","LastExonStart","LastExonEnd"))

Data <- subset(Data, Group==2)
Data$Mark[Data$Probabilty >= 0.05] <-"High"
Data$Mark[Data$Probabilty < 0.05] <-"Low"
table(Data$Mark)

data1 <- merge(Score, Data, by=c("chrom","name","strand"))
data1$Marker[data1$StopCodon=="TAA"] <- "No"
data1$Marker[data1$StopCodon=="TAG"] <- "No"
data1$Marker[data1$StopCodon=="TGA"] <- "Yes"

Yes <- subset(data1, Marker=="Yes")
No <- subset(data1, Marker=="No")
data1 <- rbind(Yes,No)
ggboxplot(data1, x="Marker", y="Score", fill = "Mark", add = "boxplot")


ggboxplot(data1, x="Marker", y="Score", fill = "Mark", add = "boxplot")

# TGA
Yes <- subset(data1, Marker=="Yes")
x=colnames(Yes)[23]
y=colnames(Yes)[7]
group=levels(factor(Yes$Mark))
Yes$Mark=factor(Yes$Mark, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
ggboxplot(Yes, x="Mark", y="Score", fill = "Mark", 
          xlab=x, ylab=y,
          legend.title=x,
          add = "boxplot")+ 
  stat_compare_means(comparisons = my_comparisons, method = "t.test")


high <- subset(Yes, Mark=="High")
low <- subset(Yes, Mark=="Low")
mean(high$Score)
mean(low$Score)



# Non-TGA
No <- subset(data1, Marker=="No")
x=colnames(No)[22]
y=colnames(No)[7]
group=levels(factor(No$Mark))
No$Mark=factor(No$Mark, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
ggboxplot(No, x="Mark", y="Score", fill = "Mark", 
          xlab=x, ylab=y,
          legend.title=x,
          add = "boxplot")+ 
  stat_compare_means(comparisons = my_comparisons, method = "t.test")

high <- subset(No, Mark=="High")
low <- subset(No, Mark=="Low")
mean(high$Score)
mean(low$Score)





ggboxplot(data1, x="Mark", y="Score", fill = "Marker", add = "boxplot")
# Low
Low <- subset(data1, Mark=="Low")
x=colnames(Low)[23]
y=colnames(Low)[7]
group=levels(factor(Low$Marker))
Low$Marker=factor(Low$Marker, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
ggboxplot(Low, x="Marker", y="Score", fill = "Marker", 
          xlab=x, ylab=y,
          legend.title=x,
          add = "boxplot")+ 
  stat_compare_means(comparisons = my_comparisons, method = "t.test")

# High
High <- subset(data1, Mark=="High")
x=colnames(High)[23]
y=colnames(High)[7]
group=levels(factor(High$Marker))
High$Marker=factor(High$Marker, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
ggboxplot(High, x="Marker", y="Score", fill = "Marker", 
          xlab=x, ylab=y,
          legend.title=x,
          add = "boxplot")+ 
  stat_compare_means(comparisons = my_comparisons, method = "t.test")












# Fisher Exact Test
data1 <- data1[,c(22,23)]
High <- subset(data1, Mark=="High")
High <- table(High$Marker)
Low <- subset(data1, Mark=="Low")
Low <- table(Low$Marker)
Table <- cbind(Low, High)
fisher.test(Table)

# Donut plot
High <- subset(data1, Mark=="High")
High <- table(High$Marker)
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

Low <- subset(data1, Mark=="Low")
Low <- table(Low$Marker)
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

Table <- read.csv("Comparision.csv")
colnames(Table)[1] <- "Codon"
Table$Codon <- factor(Table$Codon, levels=c("TGA","TAA","TAG"))
ggplot(data=Table, aes(x=Codon, y=OR, color=Codon)) + geom_bar(stat="identity", fill="white") + ylim(-0.3,0.3)















