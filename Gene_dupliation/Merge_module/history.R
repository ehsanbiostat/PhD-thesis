head(data.frame(c))
dim(data.frame(c))
c = data.frame(c)
c1 =c[!duplicated(c),]
head(c1)
dim(c1)
c1[1:15,]
c[1:15,]
c1[c1[,1] == "USV1"]
c1[c1[,2] == "USV1"]
c1[c1[,2] == "USV1",]
c1[c1[,1] == "USV1",]
c1[c1[,1] == "CIN2",]
c1[c1[,2] == "CIN2",]
head(c1)
c2 = data.frame(c1[,2],c1[,1])
head(c2)
dim(c1)
dim(c2)
colnames(c2) = c("Standard", "Systematic")
dim(c2)
head(c2)3
head(c2)
write.table(c2, file = "~/../../group/biocomp/users/ehsab/Yeast_interaction_data/Mapping.txt", sep="\t",row.names=F,quote=F)
ref = c2
a = read.table("~/../../group/biocomp/users/ehsab/Yeast_interaction_data/RegulationTwoColumnTable_Documented_201110.tsv")
a = read.table("~/../../group/biocomp/users/ehsab/Yeast_interaction_data/RegulationTwoColumnTable_Documented_201110.txt")
a = read.table("~/../../group/biocomp/users/ehsab/Yeast_interaction_data/RegulationTwoColumnTable_Documented_201110.txt")
a = read.table("~/../../group/biocomp/users/ehsab/Yeast_interaction_data/RegulationTwoColumnTable_Documented_2011109.txt")
dim(a)
head(a)
colnames(a)[1] = c("Standard")
head(a)
k = merge(ref, a, by="Standard")
head(k)
colnames(ref)
colnames(a)
k = merge(a, ref, by="Standard")
head(k)
head(ref)
head(a)
a = toupper(a)
head(a)
a = read.table("~/../../group/biocomp/users/ehsab/Yeast_interaction_data/RegulationTwoColumnTable_Documented_2011109.txt")
dim(a)
head(a)
a1 = toupper(a)
dim(a1)
a1 = toupper(a[1])
dim(a1)
head(a1)
a1 = toupper(as.matrix(a))
a1
savehistory()
getpwd()
getpwd
getwd
getwd()
a1 = data.frame(toupper(as.matrix(a)))
a1
colnames(a1)[1]=c("Standard")
head(a1)
k = merge(ref, a1, by="Standard")
head(k)
k = k[,-1]
head(k)
colnames(a1)[2]=c("Standard")
k = merge(ref, a1, by="Standard")
colnames(k)[2]=c("Standard")
k = merge(ref, a1, by="Standard")
head(k)
head(ref)
colnames(k)=c("S","Standard")
k = merge(ref, a1, by="Standard")
head(ref)
head(k)
k = merge(ref, k, by="Standard")
head(k)
k[1:30,]
k[1:100,]
k[1:1000,]
k[1:500,]
k[1:200,]
k[1:100,]
k[1:150,]
savehistory()
a1 = data.frame(toupper(as.matrix(a)))
dim(a1)
dim(ref)
k = a1
colnames(k)[1]=c("Standard")
head(k)
k = merge(ref, k, by="Standard")
dim(k)
head(k)
k[1:200,]
k[1:10000,]
k[1:40000,]
k[1:dim(k)[1],]
k[dim(k)[1],]
k = k[,-1]
k
b = read.table("~/../../group/biocomp/users/ehsab/Yeast_interaction_data/interaction_data_modified.txt")
dim(b)
b1 = b[b[,3] == "physical interactions",]
dim(b1)
b = read.table("~/../../group/biocomp/users/ehsab/Yeast_interaction_data/interaction_data_modified.txt", sep="\t")
dim(b)
b1 = b[b[,3] == "physical interactions",]
dim(b1)
c = rbind(k, b1)
colnames(b1)=colnames(k)
c = rbind(k, b1)
c = rbind(k, b1[,-3])
dim(c)
dim(c[!duplicated(c),])
dim(c)dim(k)
dim(k)
dim(b1)
dim(b)
dim(c[!duplicated(c),])
dim(c)
dim(c[!duplicated(c),])/dim(c)
d = rbind(k, b[,-3])
colnames(b)=colnames(k)
d = rbind(k, b[,-3])
dim(d)
dim(d[!duplicated(d),])/dim(d)
dim(d[!duplicated(d),])
d = rbind(k, b)
dim(d)
head(d)
b[c(2,1),-3]
head( b[,-3])
head( b[,c(-3,2,1)])
head( b[,-3&c(2,1)])
b2 = data.frame(b[,2],b[,1])
dim(b2)
colnames(b2)=colnames(k)
d = rbind(k, b2)
dim(d)
dim(d[!duplicated(d),])
length(unique(c(k[,1],k[,2])))
dim(b)
dim(b1)
length(unique(c(b1[,1],b1[,2])))
dim(k)
length(unique(c(k[,1],k[,2])))
length(unique(c(b[,1],b[,2])))
dim(b)
b2 = b[b[,3] == "genetic interactions",]
dim(b2)
197706+124924-322630
length(unique(c(b2[,1],b2[,2])))
savehistory()
load("~/Gene_duplication/Yeast/GRN/Results/R_images/20.08.13.RData")
library(igraph)
library(Matrix)
library(network)
library(randomForest)
library(intergraph)
library(lattice)
library(sm)
library(mvnormtest) # Multivariable normality test
library(nortest)# Shapiro - Francia normality test for large sample data sets
library(QuantPsyc)# Standardized regression beta coefficients
library(GO.db)
plot(overlap.pair.in, overlap.pair.out)
plot(overlap.pair.out, overlap.pair.in, xlab = 'Overlapping percentage for outgoing', ylab = 'Overlapping percentage for incoming')
abline(h = median(overlap.pair.in), v = median(overlap.pair.out))
jpeg("~/Gene_duplication/Yeast/GRN/Results/Plots/overlap.pair.out.in.jpeg")
plot(overlap.pair.out, overlap.pair.in, xlab = 'Overlapping percentage for outgoing', ylab = 'Overlapping percentage for incoming')
abline(h = median(overlap.pair.in), v = median(overlap.pair.out))
dev.off()
?jpeg()
?jpeg()
jpeg("~/Gene_duplication/Yeast/GRN/Results/Plots/overlap.pair.out.in.jpeg", width = 1200, height = 800, pointsize = 18)
plot(overlap.pair.out, overlap.pair.in, xlab = 'Overlapping percentage for outgoing', ylab = 'Overlapping percentage for incoming')
abline(h = median(overlap.pair.in), v = median(overlap.pair.out))
dev.off()
plot(overlap.pair.in, overlap.go.in.dup.TF.BP, xlab = 'Overlapping percentage for incoming', ylab = 'Overlapping percentage for incoming GO')
plot(overlap.pair.out, overlap.go.out.dup.TF.BP, xlab = 'Overlapping percentage for outgoing', ylab = 'Overlapping percentage for outgoing GO')
?factor
compare.less.more.overlap.go.TF.BP.y = c(rep("Less than median", 7), rep("More than median", 7))
compare.less.more.overlap.go.TF.CC.y = c(rep("Less than median", 7), rep("More than median", 7))
boxplot(compare.less.more.overlap.go.TF.CC ~ compare.less.more.overlap.go.TF.CC.y)
sm.density.compare(compare.less.more.overlap.go.TF.CC, compare.less.more.overlap.go.TF.CC.y)
colfill<-c(2:(2+length(levels(compare.less.more.overlap.go.TF.CC.y))))
compare.less.more.overlap.go.TF.CC
compare.less.more.overlap.go.TF.CC.y = c(rep(0, 7), rep(1, 7))
compare.less.more.overlap.go.TF.CC.factor = factor(compare.less.more.overlap.go.TF.CC.y, level = c(0, 1), labels = c("Less than median", "More than median"))
sm.density.compare(compare.less.more.overlap.go.TF.CC, compare.less.more.overlap.go.TF.CC.y)
colfill<-c(2:(2+length(levels(compare.less.more.overlap.go.TF.CC.factor))))
legend(locator(1), levels(overlap.go.in.out.TF.factor), fill=colfill)
sm.density.compare(compare.less.more.overlap.go.TF.CC, compare.less.more.overlap.go.TF.CC.y)
colfill<-c(2:(2+length(levels(compare.less.more.overlap.go.TF.CC.factor))))
legend(locator(1), levels(compare.less.more.overlap.go.TF.CC.factor), fill=colfill)
compare.less.more.overlap.go.TF.BP.y = c(rep(0, 7), rep(1, 7))
compare.less.more.overlap.go.TF.BP.factor = factor(compare.less.more.overlap.go.TF.BP.y, level = c(0, 1), labels = c("Less than median", "More than median"))
par(mfrow c= (2, 1))
par(mfrow = c(2, 1))
sm.density.compare(compare.less.more.overlap.go.TF.CC, compare.less.more.overlap.go.TF.CC.y)
colfill<-c(2:(2+length(levels(compare.less.more.overlap.go.TF.CC.factor))))
legend(locator(1), levels(compare.less.more.overlap.go.TF.CC.factor), fill=colfill)
sm.density.compare(compare.less.more.overlap.go.TF.BB, compare.less.more.overlap.go.TF.BB.y)
colfill<-c(2:(2+length(levels(compare.less.more.overlap.go.TF.BB.factor))))
legend(locator(1), levels(compare.less.more.overlap.go.TF.BB.factor), fill=colfill)
sm.density.compare(compare.less.more.overlap.go.TF.BP, compare.less.more.overlap.go.TF.BP.y)
colfill<-c(2:(2+length(levels(compare.less.more.overlap.go.TF.BP.factor))))
legend(locator(1), levels(compare.less.more.overlap.go.TF.BP.factor), fill=colfill)
par(mfrow = c(1, 1))
sm.density.compare(compare.less.more.overlap.go.TF.BP, compare.less.more.overlap.go.TF.BP.y)
colfill<-c(2:(2+length(levels(compare.less.more.overlap.go.TF.BP.factor))))
legend(locator(1), levels(compare.less.more.overlap.go.TF.BP.factor), fill=colfill)
?sm.density.compare
sm.density.compare(compare.less.more.overlap.go.TF.BP, compare.less.more.overlap.go.TF.BP.y, model = "equal")
colfill<-c(2:(2+length(levels(compare.less.more.overlap.go.TF.BP.factor))))
legend(locator(1), levels(compare.less.more.overlap.go.TF.BP.factor), fill=colfill)
?sm.density.compare
sm.density.compare(compare.less.more.overlap.go.TF.BP, compare.less.more.overlap.go.TF.BP.y, model = "equal", nboot = 100)
colfill<-c(2:(2+length(levels(compare.less.more.overlap.go.TF.BP.factor))))
legend(locator(1), levels(compare.less.more.overlap.go.TF.BP.factor), fill=colfill)
sm.density.compare(compare.less.more.overlap.go.TF.BP, compare.less.more.overlap.go.TF.BP.y, model = "equal", nboot = 1000)
colfill<-c(2:(2+length(levels(compare.less.more.overlap.go.TF.BP.factor))))
sm.density.compare(compare.less.more.overlap.go.TF.BP, compare.less.more.overlap.go.TF.BP.y, main = "A")
sm.density.compare(compare.less.more.overlap.go.TF.BP, compare.less.more.overlap.go.TF.BP.y, xlab = "A")
sm.density.compare(compare.less.more.overlap.go.TF.BP, compare.less.more.overlap.go.TF.BP.y, xlab = "Overlap BP GO terms percentage for two paralogous genes")
colfill<-c(2:(2+length(levels(compare.less.more.overlap.go.TF.BP.factor))))
legend(locator(1), levels(compare.less.more.overlap.go.TF.BP.factor), fill=colfill)
jpeg("~/Gene_duplication/Yeast/GRN/Results/Plots/overlap.GO.BP.pair.out.in.jpeg", width = 1200, height = 800, pointsize = 18)
sm.density.compare(compare.less.more.overlap.go.TF.BP, compare.less.more.overlap.go.TF.BP.y, xlab = "Overlap BP GO terms percentage for two paralogous genes")
colfill<-c(2:(2+length(levels(compare.less.more.overlap.go.TF.BP.factor))))
legend(locator(1), levels(compare.less.more.overlap.go.TF.BP.factor), fill=colfill)
dev.off()
legend(levels(compare.less.more.overlap.go.TF.BP.factor), fill=colfill)
?legend
legend('"topleft"', levels(compare.less.more.overlap.go.TF.BP.factor), fill=colfill)
legend("topleft", levels(compare.less.more.overlap.go.TF.BP.factor), fill=colfill)
sm.density.compare(compare.less.more.overlap.go.TF.BP, compare.less.more.overlap.go.TF.BP.y, xlab = "Overlap BP GO terms percentage for two paralogous genes")
colfill<-c(2:(2+length(levels(compare.less.more.overlap.go.TF.BP.factor))))
legend(levels(compare.less.more.overlap.go.TF.BP.factor), fill=colfill)
legend('"topleft"', levels(compare.less.more.overlap.go.TF.BP.factor), fill=colfill)
legend("topleft", levels(compare.less.more.overlap.go.TF.BP.factor), fill=colfill)
jpeg("~/Gene_duplication/Yeast/GRN/Results/Plots/overlap.GO.BP.pair.out.in.jpeg", width = 1200, height = 800, pointsize = 18)
sm.density.compare(compare.less.more.overlap.go.TF.BP, compare.less.more.overlap.go.TF.BP.y, xlab = "Overlap BP GO terms percentage for two paralogous genes")
colfill<-c(2:(2+length(levels(compare.less.more.overlap.go.TF.BP.factor))))
legend("topleft", levels(compare.less.more.overlap.go.TF.BP.factor), fill=colfill)
dev.off()
jpeg("~/Gene_duplication/Yeast/GRN/Results/Plots/overlap.GO.CC.pair.out.in.jpeg", width = 1200, height = 800, pointsize = 18)
sm.density.compare(compare.less.more.overlap.go.TF.CC, compare.less.more.overlap.go.TF.CC.y, xlab = "Overlap CC GO terms percentage for two paralogous genes")
colfill<-c(2:(2+length(levels(compare.less.more.overlap.go.TF.CC.factor))))
legend("topleft", levels(compare.less.more.overlap.go.TF.CC.factor), fill=colfill)
dev.off()
jpeg("~/Gene_duplication/Yeast/GRN/Results/Plots/overlap.GO.BP.pair.out.in.jpeg", width = 1200, height = 800, pointsize = 18)
sm.density.compare(compare.less.more.overlap.go.TF.BP, compare.less.more.overlap.go.TF.BP.y, xlab = "Overlap BP GO terms percentage target genes of two paralogous")
colfill<-c(2:(2+length(levels(compare.less.more.overlap.go.TF.BP.factor))))
legend("topleft", levels(compare.less.more.overlap.go.TF.BP.factor), fill=colfill)
dev.off()
jpeg("~/Gene_duplication/Yeast/GRN/Results/Plots/overlap.GO.CC.pair.out.in.jpeg", width = 1200, height = 800, pointsize = 18)
sm.density.compare(compare.less.more.overlap.go.TF.CC, compare.less.more.overlap.go.TF.CC.y, xlab = "Overlap CC GO terms percentage for target genes of two paralogous")
colfill<-c(2:(2+length(levels(compare.less.more.overlap.go.TF.CC.factor))))
legend("topleft", levels(compare.less.more.overlap.go.TF.CC.factor), fill=colfill)
dev.off()
jpeg("~/Gene_duplication/Yeast/GRN/Results/Plots/overlap.GO.BP.pair.out.in.jpeg", width = 1200, height = 800, pointsize = 18)
sm.density.compare(compare.less.more.overlap.go.TF.BP, compare.less.more.overlap.go.TF.BP.y, xlab = "Overlap BP GO terms percentage for target genes of two paralogous")
colfill<-c(2:(2+length(levels(compare.less.more.overlap.go.TF.BP.factor))))
legend("topleft", levels(compare.less.more.overlap.go.TF.BP.factor), fill=colfill)
dev.off()
head(GO_resultgenes.overlap.pair.less.median.go.BP.out.all, 20)
head(GO_resultgenes.overlap.pair.less.median.go.BP.out.all, 25)
head(GO_resultgenes.overlap.pair.more.median.go.BP.out.all, 25)
head(GO_result_against_more.genes.overlap.pair.less.median.go.BP.out.all, 20)
head(GO_result_against_less.genes.overlap.pair.more.median.go.BP.out.all, 20)
ls()
sequnce.identity.yeastract.raw = read.delim("~/Gene_duplication/Yeast_interaction_data/Sequence/alignscore.txt", header = F)
colnames(sequnce.identity.yeastract.raw)[1] = c("Systematic")
sequnce.identity.yeastract.raw = merge(ref.yeastract.SID, sequnce.identity.yeastract.raw, by = "Systematic")
head(sequnce.identity.yeastract.raw)
sequnce.identity.yeastract.raw = read.delim("~/Gene_duplication/Yeast_interaction_data/Sequence/alignscore.txt", header = F)
dim(sequnce.identity.yeastract.raw)
colnames(sequnce.identity.yeastract.raw)[1] = c("Systematic")
sequnce.identity.yeastract.raw = merge(ref.yeastract.SID, sequnce.identity.yeastract.raw, by = "Systematic")
dim(sequnce.identity.yeastract.raw)
head(sequnce.identity.yeastract.raw)
sequnce.identity.yeastract.raw = sequnce.identity.yeastract.raw[, -1]
head(sequnce.identity.yeastract.raw)
colnames(sequnce.identity.yeastract.raw)[4] = c("Systematic")
head(sequnce.identity.yeastract.raw)
sequnce.identity.yeastract.raw = merge(ref.yeastract.SID, sequnce.identity.yeastract.raw, by = "Systematic")
dim(sequnce.identity.yeastract.raw)
434*2
yeastract.ref
ref.yeastract
yeastract.conv
yeastract.ref
ref.yeastract.SID
ref.yeastract
dim(ref.yeastract.SID)
ref.yeastract.SID
head(ref.yeastract.SID, 50)
sequnce.identity.yeastract.raw = read.delim("~/Gene_duplication/Yeast_interaction_data/Sequence/alignscore.txt", header = F)
head(sequnce.identity.yeastract.raw)
which(ref.yeastract.SID[, 1] %in% sequnce.identity.yeastract.raw[, 1])
which!(ref.yeastract.SID[, 1] %in% sequnce.identity.yeastract.raw[, 1])
!(ref.yeastract.SID[, 1] %in% sequnce.identity.yeastract.raw[, 1])
(ref.yeastract.SID[, 1] %in% sequnce.identity.yeastract.raw[, 1])
which(!(ref.yeastract.SID[, 1] %in% sequnce.identity.yeastract.raw[, 1]))
which(!(sequnce.identity.yeastract.raw[, 1] %in% sequnce.identity.yeastract.raw[, 1]))
which((sequnce.identity.yeastract.raw[, 1] %in% sequnce.identity.yeastract.raw[, 1]))
which((sequnce.identity.yeastract.raw[, 1] %in% ref.yeastract.SID[, 1]))
which(!(sequnce.identity.yeastract.raw[, 1] %in% ref.yeastract.SID[, 1]))
sequnce.identity.yeastract.raw[which(!(sequnce.identity.yeastract.raw[, 1] %in% ref.yeastract.SID[, 1])), ]
sequnce.identity.yeastract.raw = read.delim("~/Gene_duplication/Yeast_interaction_data/Sequence/alignscore.txt", header = F)
colnames(sequnce.identity.yeastract.raw)[1] = c("Systematic")
sequnce.identity.yeastract.raw = merge(ref.yeastract.SID, sequnce.identity.yeastract.raw, by = "Systematic")
rm(sequnce.identity.yeastract.raw)
sequence.identity.yeastract.raw = read.delim("~/Gene_duplication/Yeast_interaction_data/Sequence/alignscore.txt", header = F)
colnames(sequence.identity.yeastract.raw)[1] = c("Systematic")
## After converting first column to the reference codes we have
sequence.identity.yeastract.raw[which(!(sequence.identity.yeastract.raw[, 1] %in% ref.yeastract.SID[, 1])), ]
sequence.identity.yeastract.raw = merge(ref.yeastract.SID, sequence.identity.yeastract.raw, by = "Systematic")
dim(sequence.identity.yeastract.raw)
head(sequence.identity.yeastract.raw)
sequence.identity.yeastract.raw = sequence.identity.yeastract.raw[, -1]
colnames(sequence.identity.yeastract.raw)[4] = c("Systematic")
sequence.identity.yeastract.raw = merge(ref.yeastract.SID, sequence.identity.yeastract.raw, by = "Systematic")
head(sequence.identity.yeastract.raw)
sequence.identity.yeastract.raw = sequence.identity.yeastract.raw[, -c(1, 4, 7)]
head(sequence.identity.yeastract.raw)
colnames(sequence.identity.yeastract.raw)[5] = c("Seq.identity")
head(sequence.identity.yeastract.raw)
sequence.identity.yeastract.raw = read.delim("~/Gene_duplication/Yeast_interaction_data/Sequence/alignscore.txt", header = F)
colnames(sequence.identity.yeastract.raw)[1] = c("Systematic")
sequence.identity.yeastract.raw[which(!(sequence.identity.yeastract.raw[, 1] %in% ref.yeastract.SID[, 1])), ]
sequence.identity.yeastract.raw = sequence.identity.yeastract.raw[, -1]
head(sequence.identity.yeastract.raw)
sequence.identity.yeastract.raw = read.delim("~/Gene_duplication/Yeast_interaction_data/Sequence/alignscore.txt", header = F)
colnames(sequence.identity.yeastract.raw)[1] = c("Systematic")
sequence.identity.yeastract.raw[which(!(sequence.identity.yeastract.raw[, 1] %in% ref.yeastract.SID[, 1])), ]
##V1      V2   V3
##25  YOR069W YKR078W 0.18
##123 YKR028W YJL098W 0.43
##174 YKL166C YJL164C 0.74
##225 YPL105C YBR172C 0.33
##237 YHR158C YGR238C 0.44
##276 YGL082W YPL191C 0.46
##289 YBR270C YJL058C 0.45
##389 YKR050W YJL129C 0.55
## After converting first column to the reference codes we have 439 paralogous genes
sequence.identity.yeastract.raw = merge(ref.yeastract.SID, sequence.identity.yeastract.raw, by = "Systematic")
head(sequence.identity.yeastract.raw)
sequence.identity.yeastract.raw = sequence.identity.yeastract.raw[, -1]
colnames(sequence.identity.yeastract.raw)[4] = c("Systematic")
head(sequence.identity.yeastract.raw)
colnames(sequence.identity.yeastract.raw)[4] = c("Systematic")
head(sequence.identity.yeastract.raw)
sequence.identity.yeastract.raw[which(!(sequence.identity.yeastract.raw[, 4] %in% ref.yeastract.SID[, 1])), ]
sequence.identity.yeastract.raw = merge(ref.yeastract.SID, sequence.identity.yeastract.raw, by = "Systematic")
dim(sequence.identity.yeastract.raw)
sequence.identity.yeastract.raw = sequence.identity.yeastract.raw[, -c(1, 4, 7)]
colnames(sequence.identity.yeastract.raw)[5] = c("Seq.identity")
head(sequence.identity.yeastract.raw)
dup.pair.ref
dup.pair
dup.pair.ref.pur
dup.pair.ref.TF
dup.pair.ref.TF.pur
dim(dup.pair.ref.TF.pur)
dup.pair.ref.pur
dim(dup.pair.ref.pur)
sequence.identity.yeastract.raw.temp = sequence.identity.yeastract.raw
head(sequence.identity.yeastract.raw.temp)
head(dup.pair.ref.pur)
sequence.identity.yeastract.raw.temp = merge(dup.pair.ref.pur, sequence.identity.yeastract.raw.temp, by = "ref")
colnames(sequence.identity.yeastract.raw.temp)[1] = c("ref")
sequence.identity.yeastract.raw.temp = merge(dup.pair.ref.pur, sequence.identity.yeastract.raw.temp, by = "ref")
head(sequence.identity.yeastract.raw.temp)
dim(sequence.identity.yeastract.raw.temp)
sequence.identity.yeastract.raw.temp = sequence.identity.yeastract.raw.temp[, -1]
dim(sequence.identity.yeastract.raw.temp)
head(sequence.identity.yeastract.raw.temp)
sequence.identity.yeastract.raw.temp = sequence.identity.yeastract.raw
colnames(sequence.identity.yeastract.raw.temp)[1] = c("ref")
sequence.identity.yeastract.raw.temp = merge(dup.pair.ref.pur, sequence.identity.yeastract.raw.temp, by = "ref")
head(sequence.identity.yeastract.raw.temp)
colnames(sequence.identity.yeastract.raw.temp)[c(1, 5)] = c("ref.x", "ref")
head(sequence.identity.yeastract.raw.temp)
sequence.identity.yeastract.raw.temp = merge(dup.pair.ref.pur, sequence.identity.yeastract.raw.temp, by = "ref")
head(sequence.identity.yeastract.raw.temp)
dim(sequence.identity.yeastract.raw.temp)
(sequence.identity.yeastract.raw.temp[, 3]
)
(sequence.identity.yeastract.raw.temp[, 6])
sum(sequence.identity.yeastract.raw.temp[, 6])
sum(sequence.identity.yeastract.raw.temp[, 3])
quntile(sequence.identity.yeastract.raw.temp[, 3])
quantile(sequence.identity.yeastract.raw.temp[, 3])
quantile(sequence.identity.yeastract.raw.temp[, 6])
head(sequence.identity.yeastract.raw.temp)
colnames(sequence.identity.yeastract.raw.temp)[1] = c("ref.y")
head(sequence.identity.yeastract.raw.temp)
sequence.identity.yeastract.raw.temp = data.frame(sequence.identity.yeastract.raw.temp[, c(4, 1, 7, 8, 2, 5, 3, 6, 7)])
head(sequence.identity.yeastract.raw.temp)
sequence.identity.yeastract.raw.temp = sequence.identity.yeastract.raw
colnames(sequence.identity.yeastract.raw.temp)[1] = c("ref")
sequence.identity.yeastract.raw.temp = merge(dup.pair.ref.pur, sequence.identity.yeastract.raw.temp, by = "ref")
colnames(sequence.identity.yeastract.raw.temp)[c(1, 5)] = c("ref.x", "ref")
sequence.identity.yeastract.raw.temp = merge(dup.pair.ref.pur, sequence.identity.yeastract.raw.temp, by = "ref")
colnames(sequence.identity.yeastract.raw.temp)[1] = c("ref.y")
head(sequence.identity.yeastract.raw.temp)
sequence.identity.yeastract.raw.temp = data.frame(sequence.identity.yeastract.raw.temp[, c(4, 1, 7, 8, 2, 5, 3, 6, 9)])
head(sequence.identity.yeastract.raw.temp)
head(sequence.identity.yeastract.raw.temp[order(sequence.identity.yeastract.raw.temp[, 7]), ])
sequence.identity.yeastract.raw = sequence.identity.yeastract.raw.temp[order(sequence.identity.yeastract.raw.temp[, 7]), ]
head(sequence.identity.yeastract.raw)
sequence.identity.yeastract.raw$pairs.x
dup.pair.ref.TF.pur
dup.pair.ref.TF.pur$pairs %in% seq.iden$pairs.x
seq.iden = sequence.identity.yeastract.raw
dup.pair.ref.TF.pur$pairs %in% seq.iden$pairs.x
seq.iden$pairs.x %in% dup.pair.ref.TF.pur$pairs
which(seq.iden$pairs.x %in% dup.pair.ref.TF.pur$pairs)
seq.iden[which(seq.iden$pairs.x %in% dup.pair.ref.TF.pur$pairs), ]
dup.pair.ref.TF.pur
seq.iden.dup.pair.ref.TF.pur = seq.iden[which(seq.iden$pairs.x %in% dup.pair.ref.TF.pur$pairs), 9]
seq.iden.dup.pair.ref.TF.pur
plot(overlap.pair.out, seq.iden.dup.pair.ref.TF.pur)
plot(overlap.pair.in, seq.iden.dup.pair.ref.TF.pur)
plot(seq.iden.dup.pair.ref.TF.pur, overlap.pair.in)
plot(seq.iden.dup.pair.ref.TF.pur, overlap.pair.out)
cor.test(seq.iden.dup.pair.ref.TF.pur, overlap.pair.out)
cor.test(seq.iden.dup.pair.ref.TF.pur, overlap.pair.in)
dup.pair.ref.TF.pur
seq.iden.dup.pair.ref.TF.pur
seq.iden[which(seq.iden$pairs.x %in% dup.pair.ref.TF.pur$pairs), ]
seq.iden.dup.pair.ref.TF.pur[index.overlap.pair.more.median]
seq.iden.dup.pair.ref.TF.pur[index.overlap.pair.less.median]
compare.less.more.seq.iden.dup.pair.ref.TF.pur = c(seq.iden.dup.pair.ref.TF.pur[index.overlap.pair.less.median], seq.iden.dup.pair.ref.TF.pur[index.overlap.pair.more.median])
compare.less.more.seq.iden.dup.pair.ref.TF.pur
boxplot(compare.less.more.seq.iden.dup.pair.ref.TF.pur ~ compare.less.more.seq.iden.dup.pair.ref.TF.pur.y)
compare.less.more.seq.iden.dup.pair.ref.TF.pur.y = c(rep(0, 7), rep(1, 7))
compare.less.more.seq.iden.dup.pair.ref.TF.pur.factor = factor(compare.less.more.seq.iden.dup.pair.ref.TF.pur.y, level = c(0, 1), labels = c("Less than median", "More than median"))
boxplot(compare.less.more.seq.iden.dup.pair.ref.TF.pur ~ compare.less.more.seq.iden.dup.pair.ref.TF.pur.y)
boxplot(compare.less.more.seq.iden.dup.pair.ref.TF.pur ~ compare.less.more.seq.iden.dup.pair.ref.TF.pur.factor)
sm.density.compare(compare.less.more.seq.iden.dup.pair.ref.TF.pur, compare.less.more.seq.iden.dup.pair.ref.TF.pur.y, xlab = "Overlap BP GO terms percentage for target genes of two paralogous")
colfill<-c(2:(2+length(levels(compare.less.more.seq.iden.dup.pair.ref.TF.pur.factor))))
legend("topleft", levels(compare.less.more.seq.iden.dup.pair.ref.TF.pur.factor), fill=colfill)
sm.density.compare(compare.less.more.seq.iden.dup.pair.ref.TF.pur, compare.less.more.seq.iden.dup.pair.ref.TF.pur.y, xlab = "Overlap BP GO terms percentage for target genes of two paralogous", method = "equal", nboot = 1000)
sm.density.compare(compare.less.more.seq.iden.dup.pair.ref.TF.pur, compare.less.more.seq.iden.dup.pair.ref.TF.pur.y, xlab = "Overlap BP GO terms percentage for target genes of two paralogous", model = "equal", nboot = 1000)
overlap.go.out.dup.TF.BP
plot(seq.iden.dup.pair.ref.TF.pur, overlap.go.out.dup.TF.BP)
plot(seq.iden.dup.pair.ref.TF.pur, overlap.go.in.dup.TF.BP)
cor.test(seq.iden.dup.pair.ref.TF.pur, overlap.go.out.dup.TF.BP)
plot(seq.iden.dup.pair.ref.TF.pur, overlap.go.out.dup.TF.BP)
plot(seq.iden.dup.pair.ref.TF.pur[index.overlap.pair.more.median], overlap.go.out.dup.TF.BP[index.overlap.pair.more.median])
plot(seq.iden.dup.pair.ref.TF.pur[index.overlap.pair.less.median], overlap.go.out.dup.TF.BP[index.overlap.pair.less.median])
cor.test(seq.iden.dup.pair.ref.TF.pur[index.overlap.pair.less.median], overlap.go.out.dup.TF.BP[index.overlap.pair.less.median])
cor.test(seq.iden.dup.pair.ref.TF.pur[index.overlap.pair.less.median], overlap.go.out.dup.TF.BP[index.overlap.pair.less.median])
cor.test(seq.iden.dup.pair.ref.TF.pur[index.overlap.pair.more.median], overlap.go.out.dup.TF.BP[index.overlap.pair.more.median])
plot(seq.iden.dup.pair.ref.TF.pur[index.overlap.pair.more.median], overlap.go.out.dup.TF.BP[index.overlap.pair.more.median])
seq.iden.dup.pair.ref.TF.pur[index.overlap.pair.more.median]
overlap.go.out.dup.TF.BP[index.overlap.pair.more.median]
index.overlap.pair.more.median
index.overlap.pair.less.median
library(igraph)
library("isa2")
library("foreach")
library("doMC")
library(igarph)
library(igraph)
data.exp = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated CORNET2.0/Integrated_SD.txt", colClasses = "numeric", nrows = 19285, comment.char = "", header = T)
T = matrix(NA, 19285, 13)
for (i in 1:13){
load(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Merge/Gene_score/Compendium_wise/abiotic/Permutation/V2/G_200_5_100",i,".RData", sep = ""))
T[,i] = gene_score_module
}
}
for (i in 1:13){{
for (i in 1:13){
load(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Merge/Gene_score/Compendium_wise/abiotic/Permutation/V2/G_200_5_100",i,".RData", sep = ""))
T[,i] = gene_score_module
}
head(T)
sum(T)
which(T[, 1] == 1)
intersect(which(T[, 1] == 1), which(T[, 2] == 1))
table(apply(T, 1, sum))
intersect(which(T[, 1] == 1), which(T[, 2] == 1))/union(which(T[, 1] == 1), which(T[, 2] == 1))
length(intersect(which(T[, 1] == 1), which(T[, 2] == 1)))/length(union(which(T[, 1] == 1), which(T[, 2] == 1)))
Jacc = matrix(NA, 13,13)
for(i in 1:13){
for(j in 1:13){
Jacc[i,j] = length(intersect(which(T[, i] == 1), which(T[, j] == 1)))/length(union(which(T[, i] == 1), which(T[, j] == 1)))
}
}
Jacc
options(width = 200)
Jacc
image(Jacc)
hclust(as.dist(Jacc))
plot(hclust(as.dist(Jacc)))
?hclust
1 - Jacc
plot(hclust(as.dist(1 - Jacc)))
plot(hclust(as.dist(1 - Jacc)), method = "ward")
plot(hclust(as.dist(1 - Jacc), method = "ward")))
plot(hclust(as.dist(1 - Jacc), method = "ward"))
Jacc
plot(hclust(as.dist(1 - Jacc), method = "median"))
plot(hclust(as.dist(1 - Jacc), method = "average"))
plot(hclust(as.dist(1 - Jacc), method = "single"))
savehistory("/ngsprojects/iwt/ehsab/Gene_duplication/Scripts/SA/Merge_module/history.R")
