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
library(igraph)
library(foreach)
library(doMC)
registerDoMC(2)
k = 6
condition.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated CORNET2.0/condition_name.txt", header = T)
NAME = names(table(condition.names))
PO = read.delim("/group/biocomp/projects/CORNET/MA_annotation/cornet_allMA_desc_120509.txt")
NAME.spe = NAME[k]
source("/ngsprojects/iwt/ehsab/Gene_duplication/Scripts/General/myImagePlot.R")
source("/ngsprojects/iwt/ehsab/Gene_duplication/Scripts/General/myImagePlot_sym.R")
## Condition scores
list = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Scripts/General/order.condition_score.txt")
list = as.vector(unlist(list))
## Read text or RData file with read.table or load function
setwd(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Condition_score/Compendium_wise/",NAME.spe,"/Raw", sep = ""))
# C = lapply(list, read.table, sep= "\t" , header=T)
C = lapply(list, function(x) mget(load(x)))
total_cond_scores = matrix(NA, length(unlist(C[1])), 19285)
for (i in 1:19285){
print(i)
total_cond_scores[, i] = as.matrix(unlist(C[i]))
}
data.exp = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated CORNET2.0/Integrated_SD.txt", colClasses = "numeric", nrows = 19285, comment.char = "", header = T)
data.exp = data.exp[, condition.names == NAME.spe]
data.exp.scale = apply(data.exp,1,scale)
data.exp.scale = t(data.exp.scale)
ref = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated CORNET2.0/Reference_SD.txt", header = T)
## Load generated graph from SA result
load(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/",NAME.spe,"/Graph_V2.RData", sep = ""))
j = 8732
jpar = 3977
memb.comun = neighborhood(g,1,j)[[1]]
memb.comun.para = neighborhood(g,1,jpar)[[1]]
memb.comun
memb.comun.para
CondRank1 = foreach(i = 1:length(memb.comun), .combine = rbind) %do%{
rank(total_cond_scores[, memb.comun[i]])
}
CondRank2 = foreach(i = 1:length(memb.comun.para), .combine = rbind) %do%{
rank(total_cond_scores[, memb.comun.para[i]])
}
rownames(CondRank2) = c(1:length(memb.comun.para))
rownames(CondRank1) = c(1:length(memb.comun))
CondRank = rbind(c(1:dim(total_cond_scores)[1]), CondRank1, CondRank2)
colnames(CondRank) = colnames(data.exp)
CondRank0 = CondRank[,order(CondRank[2,], decreasing = T)]
CondRank = rbind(c(1:dim(total_cond_scores)[1]), CondRank2, CondRank1)
colnames(CondRank) = colnames(data.exp)
CondRank01 = CondRank[,order(CondRank[2,], decreasing = T)]
myImagePlot_sym(data.exp.scale[c(memb.comun, memb.comun.para),CondRank0[1,]])
myImagePlot_sym(data.exp.scale[c(memb.comun, memb.comun.para),CondRank01[1,]])
paralogous = read.table("~/Gene_duplication/Seed_genes/paralogous_total.txt", header = T)
Inter1 = intersect(memb.comun.para, paralogous[paralogous$Ref %in% memb.comun, "Anch"])
Inter2 = intersect(memb.comun, paralogous[paralogous$Ref %in% memb.comun.para, "Anch"])
aa = paralogous[paralogous$Ref %in% Inter2, ]
aa
bb = aa[aa$Anch %in% Inter1,]
colnames(bb)[1] = "References"
Inter_name = merge(ref, bb, by = "References")
colnames(Inter_name)[c(1,3)] = c("Ref", "References")
Inter_name = merge(ref, Inter_name, by = "References")
Inter_name
i = 87323977
write.table(Inter_name, file = paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Local/Compendium_wise/",NAME.spe,"/Score99/Module_Anchorpoints_V3/Module",i,".txt", sep = ""), row.names = F, quote = F, sep = "\t")
write.table(Inter_name, file = paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Local/Compendium_wise/",NAME.spe,"/Score99/Module_Anchorpoints_V2/Module",i,".txt", sep = ""), row.names = F, quote = F, sep = "\t")
k = i
Result = read.table(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Local/Compendium_wise/",NAME.spe,"/Score99/Module_Anchorpoints_V2/Module",k,".txt", sep = ""), header = T)
j = Result$Ref #Gene.y
jpar = Result$References #Gene.x
CondRank1 = foreach(i = 1:length(j), .combine = rbind) %do%{
rank(total_cond_scores[, j[i]])
}
CondRank2 = foreach(i = 1:length(jpar), .combine = rbind) %do%{
rank(total_cond_scores[, jpar[i]])
}
CondRank1_score = t(total_cond_scores[, j])
CondRank2_score = t(total_cond_scores[, jpar])
colnames(CondRank1) = colnames(CondRank2) = colnames(data.exp)
CondRank111 = apply(CondRank1, 2, median)
CondRank11 = apply(CondRank1, 2, mean)
CondRank211 = apply(CondRank2, 2, median)
CondRank21 = apply(CondRank2, 2, mean)
TOP1 = names(sort(CondRank111, decreasing=T)[1:20])
TOP2 = names(sort(CondRank211, decreasing=T)[1:20])
TOP1.C = foreach (i = 1:length(TOP1), .combine = c) %do% {
round(as.numeric(substring(TOP1[i],2)))
}
TOP2.C = foreach (i = 1:length(TOP2), .combine = c) %do% {
round(as.numeric(substring(TOP2[i],2)))
}
CondRank = rbind(c(1:dim(total_cond_scores)[1]), CondRank1, CondRank2)
colnames(CondRank) = colnames(data.exp)
ORDER = which(Result$Gene.y %in% R$Result.Gene.x[k])[1]
CondRank0 = CondRank[,order(CondRank[ORDER+1, ], decreasing = T)]
CondRank = rbind(c(1:dim(total_cond_scores)[1]), CondRank2, CondRank1)
colnames(CondRank) = colnames(data.exp)
ORDER = which(Result$Gene.y %in% R$Result.Gene.x[k])[1]
CondRank01 = CondRank[,order(CondRank[ORDER+1, ], decreasing = T)]
R = read.table(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Local/Compendium_wise/",NAME.spe,"/Score99/Module_center_node95_V3.txt", sep = ""), header = T)
R = read.table(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Local/Compendium_wise/",NAME.spe,"/Score99/Module_center_node95_V2.txt", sep = ""), header = T)
head(R)
ORDER = 1
head(CondRank)
ORDER = which(Result$Gene.y %in% R$Result.Gene.x[k])[1]
CondRank0 = CondRank[,order(CondRank[ORDER+1, ], decreasing = T)]
CondRank = rbind(c(1:dim(total_cond_scores)[1]), CondRank2, CondRank1)
colnames(CondRank) = colnames(data.exp)
ORDER = which(Result$Gene.y %in% R$Result.Gene.x[k])[1]
CondRank01 = CondRank[,order(CondRank[ORDER+1, ], decreasing = T)]
myImagePlot_sym(data.exp.scale[c(j, jpar),CondRank0[1,]])
myImagePlot_sym(data.exp.scale[c(j, jpar),CondRank0[1,1:20]])
myImagePlot_sym(data.exp.scale[c(j, jpar),CondRank01[1,1:20]])
myImagePlot_sym(data.exp.scale[c(j, jpar),CondRank0[1,1:20]])
myImagePlot_sym(data.exp.scale[c(j, jpar),CondRank01[1,1:20]])
myImagePlot_sym(data.exp.scale[c(j, jpar),c(CondRank0[1,1:20],CondRank01[1,1:20])])
CondRank1 = foreach(i = 1:length(j), .combine = rbind) %do%{
rank(total_cond_scores[, j[i]])
}
CondRank2 = foreach(i = 1:length(jpar), .combine = rbind) %do%{
rank(total_cond_scores[, jpar[i]])
}
CondRank1_score = t(total_cond_scores[, j])
CondRank2_score = t(total_cond_scores[, jpar])
CondRank1
CondRank2
CondRank = rbind(c(1:dim(total_cond_scores)[1]), CondRank1, CondRank2)
colnames(CondRank) = colnames(data.exp)
ORDER
ORDER = 1
CondRank0 = CondRank[,order(CondRank[ORDER+1, ], decreasing = T)]
CondRank = rbind(c(1:dim(total_cond_scores)[1]), CondRank2, CondRank1)
colnames(CondRank) = colnames(data.exp)
ORDER
CondRank01 = CondRank[,order(CondRank[ORDER+1, ], decreasing = T)]
myImagePlot_sym(data.exp.scale[c(j, jpar),c(CondRank0[1,1:20],CondRank01[1,1:20])])
j
jpar
savehistory("/ngsprojects/iwt/ehsab/Gene_duplication/Scripts/SA/History/20.10.14.11.40.r")
