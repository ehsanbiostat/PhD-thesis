
library(igraph)
library(foreach)
library(doMC)
options(width = 200)
condition.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/condition_name.txt", header = T)
NAME = names(table(condition.names))
NAME = NAME[-4]
# paralogous = read.table("~/Gene_duplication/Seed_genes/paralogous_total.txt", header = T)
ref = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Reference_SD.txt", header = T)
paralogous = read.table("~/Gene_duplication/Seed_genes/paralogous_NoFilter.txt", header = T)
GeneFamilyDup = read.delim("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/GeneFamilyZhenFeb16SSDLSD.txt")


RR = Score[Score$GeneFamilySizeRoot > 0.5,]

(RR$SSDGeneFamily - RR$LSDSSDGeneFamily)/RR$GeneFamilyNumber
sum((RR$SSDGeneFamily - RR$LSDSSDGeneFamily)/RR$GeneFamilyNumber > 0.66667)
666/1407
hist((RR$SSDGeneFamily - RR$LSDSSDGeneFamily)/RR$GeneFamilyNumber)
head(RR)
RR$Cluster.overlap
RR = RR[RR$Cluster.overlap > 50,]
RR$Cluster.overlap
hist((RR$SSDGeneFamily - RR$LSDSSDGeneFamily)/RR$GeneFamilyNumber)
sum((RR$SSDGeneFamily - RR$LSDSSDGeneFamily)/RR$GeneFamilyNumber > 0.66667)
(RR$SSDGeneFamily - RR$LSDSSDGeneFamily)/RR$GeneFamilyNumber
318/900
RR = RR[RR$Cluster.overlap > 100,]
sum((RR$SSDGeneFamily - RR$LSDSSDGeneFamily)/RR$GeneFamilyNumber > 0.66667)
(RR$SSDGeneFamily - RR$LSDSSDGeneFamily)/RR$GeneFamilyNumber
244/750
RR = RR[RR$Cluster.overlap > 200,]
(RR$SSDGeneFamily - RR$LSDSSDGeneFamily)/RR$GeneFamilyNumber
sum((RR$SSDGeneFamily - RR$LSDSSDGeneFamily)/RR$GeneFamilyNumber > 0.66667)
184/600
hist((RR$SSDGeneFamily - RR$LSDSSDGeneFamily)/RR$GeneFamilyNumber)
head(RR)
RR = RR[RR$Cluster.overlap > 500,]
hist((RR$SSDGeneFamily - RR$LSDSSDGeneFamily)/RR$GeneFamilyNumber)
((RR$SSDGeneFamily - RR$LSDSSDGeneFamily)/RR$GeneFamilyNumber)
RR
sum((RR$SSDGeneFamily - RR$LSDSSDGeneFamily)/RR$GeneFamilyNumber > 0.66667)
((RR$SSDGeneFamily - RR$LSDSSDGeneFamily)/RR$GeneFamilyNumber)
128/360
RR = RR[RR$Cluster.overlap > 1000,]
((RR$SSDGeneFamily - RR$LSDSSDGeneFamily)/RR$GeneFamilyNumber)
sum((RR$SSDGeneFamily - RR$LSDSSDGeneFamily)/RR$GeneFamilyNumber > 0.66667)
82/191
RR = RR[RR$Cluster.overlap > 1200,]
sum((RR$SSDGeneFamily - RR$LSDSSDGeneFamily)/RR$GeneFamilyNumber > 0.66667)
((RR$SSDGeneFamily - RR$LSDSSDGeneFamily)/RR$GeneFamilyNumber)
66/140
head(Score)
library(ggplot2)
data = data.frame("Jaccard" = Score$Cluster.overlap/(Score$Gene.cluster.size + Score$Paralogous))
ggplot()
(Score$SSDGeneFamily - Score$LSDSSDGeneFamily)/Score$GeneFamilyNumber)
(Score$SSDGeneFamily - Score$LSDSSDGeneFamily)/Score$GeneFamilyNumber
cor((Score$Cluster.overlap/(Score$Gene.cluster.size + Score$Paralogous)), ((Score$SSDGeneFamily - Score$LSDSSDGeneFamily)/Score$GeneFamilyNumber))

plot((Score$Cluster.overlap/(Score$Gene.cluster.size + Score$Paralogous)), ((Score$SSDGeneFamily - Score$LSDSSDGeneFamily)/Score$GeneFamilyNumber))
data = data.frame("Conserved" = (Score$Cluster.overlap/(Score$Gene.cluster.size + Score$Paralogous)), "SSD" = ((Score$SSDGeneFamily - Score$LSDSSDGeneFamily)/Score$GeneFamilyNumber),
	"Group" = 0)
data[Score$Gene.cluster.size > 1900 & Score$Paralogous.cluster.size > 1900, "Group"] = 1
ggplot(data, aes (x = Conserved, y = SSD, color = Group)) + geom_point()
savehistory("/ngsprojects/iwt/ehsab/Gene_duplication/Scripts/SA/History/SSDvsLSD.R")


































