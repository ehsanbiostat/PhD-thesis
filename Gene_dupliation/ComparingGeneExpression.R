

library(igraph)
library(foreach)
library(doMC)
library(ggplot2)
# registerDoMC(2)
k = 6
condition.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/condition_name.txt", header = T)
NAME = names(table(condition.names))
NAME.spe = NAME[k]
source("/ngsprojects/iwt/ehsab/Gene_duplication/Scripts/General/myImagePlot.R")
source("/ngsprojects/iwt/ehsab/Gene_duplication/Scripts/General/myImagePlot_sym.R")
## Condition scores

data.exp = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Integrated_SD.txt", colClasses = "numeric", nrows = 19285, comment.char = "", header = T)
data.exp.scale = apply(data.exp,1,scale)
data.exp.scale = t(data.exp.scale)
data.exp.scale = apply(data.exp.scale, 2, scale)
data.exp.scale = data.exp.scale[, condition.names == NAME.spe]
data.exp.compen = data.exp[, condition.names == NAME.spe]


data.exp.scale.med = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Integrated_SD_Norm_Med.txt", colClasses = "numeric", nrows = 19285, comment.char = "", header = T)

data.exp.scale.med = data.exp.scale.med[, condition.names == NAME.spe]

Sig = ifelse(data.exp.scale > 0, 1, ifelse(data.exp.scale < 0, -1, 0))
Sig.med = ifelse(data.exp.scale.med > 0, 1, ifelse(data.exp.scale.med < 0, -1, 0))

Sig.med.zero = ifelse(data.exp.scale.med > 0, 1, ifelse(data.exp.scale.med < 0, 0, 0))


ref = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Reference_SD.txt", header = T)
rownames(ref) = ref[,1]
# paralogous = read.table("~/Gene_duplication/Seed_genes/paralogous_total.txt", header = T)
paralogous = read.table("~/Gene_duplication/Seed_genes/paralogous_NoFilter.txt", header = T)


memb.comun = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/AT1G50090_95.txt") # Root
memb.comun.para = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/AT3G19710_95.txt") # Leaf

memb.comun = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/Geneset1_95.txt") # Leaf
memb.comun.para = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/Geneset2_95.txt") # Pollen

memb.comun = ref[as.character(memb.comun[,1]),2]
memb.comun.para = ref[as.character(memb.comun.para[,1]),2] 

Inter1 = intersect(memb.comun.para, paralogous[paralogous$Ref %in% memb.comun, "Anch"])
Inter2 = intersect(memb.comun, paralogous[paralogous$Ref %in% memb.comun.para, "Anch"])
aa = paralogous[paralogous$Ref %in% Inter2, ]
bb = aa[aa$Anch %in% Inter1,]
colnames(bb)[1] = "References"
Inter_name = merge(ref, bb, by = "References")
colnames(Inter_name)[c(1,3)] = c("Ref", "References")
Inter_name = merge(ref, Inter_name, by = "References")
	
j = Inter_name$Ref # Leaf ## Root
jpar = Inter_name$References # Pollen ## Leaf


# Condition scores for both modules
load("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Condition_score/Compendium_wise/development/Raw/Geneset1_pollen.RData")
CondRank1 = cond_scores
load("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Condition_score/Compendium_wise/development/Raw/Geneset1_leaf.RData")
CondRank2 = cond_scores
CondRank = rbind(c(1:length(CondRank1)), data.frame(t(CondRank1)), data.frame(t(CondRank2)))
CondRankPollen = CondRank[,order(CondRank[2, ], decreasing = T)]
CondRankLeaf = CondRank[,order(CondRank[3, ], decreasing = T)]
myImagePlot_sym(data.exp.scale[c(j, jpar), unlist(c(CondRankLeaf[1,1:10], CondRankPollen[1,1:10]))])


# For entire module
data.exp.Leaf = data.exp.compen[memb.comun,]
data.exp.Pollen = data.exp.compen[memb.comun.para,]
denPlot = data.frame("Exp" = as.numeric(unlist(data.exp.Pollen)), "Module" = rep("Pollen", length(as.numeric(unlist(data.exp.Pollen)))))
denPlot = rbind(denPlot, data.frame("Exp" = as.numeric(unlist(data.exp.Leaf)), "Module" = rep("Leaf", length(as.numeric(unlist(data.exp.Leaf))))))

jpeg ("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/GenomeDominance/DensityAll_90.jpeg", width = 11, height = 7, units = "in", res = 200, pointsize=8)
ggplot(denPlot, aes(x = Exp, fill = Module)) + geom_density(alpha = 0.3)
dev.off()

# KS test 
ks.test(as.numeric(unlist(data.exp.Pollen)), as.numeric(unlist(data.exp.Leaf)), alternative = c("greater"))

# Comparison per condition
percondition = data.frame("median" = apply(data.exp.Leaf, 2, median), "Module" = rep("Leaf", 134))
percondition = rbind(percondition, data.frame("median" = apply(data.exp.Pollen, 2, median), "Module" = rep("Pollen", 134)))
percondition["cond"] = rep(1:134,2)

jpeg ("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/GenomeDominance/PerConditionAll.jpeg", width = 11, height = 7, units = "in", res = 200, pointsize=8)
ggplot(percondition, aes(y = median, x = cond, color = Module)) + geom_line()
dev.off()


# T-test for comparing gene expression of each compendium
TTest.pvalue = foreach(i = 1:ncol(data.exp.Leaf), .combine = c) %do% {
	t.test(data.exp.Leaf[,i], data.exp.Pollen[,i], alternative = "greater")$p.value
}
TTest.pvalue = p.adjust(TTest.pvalue , method = "fdr")

jpeg ("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/GenomeDominance/PerConditionAllTtest.jpeg", width = 11, height = 7, units = "in", res = 200, pointsize=8)
ggplot(percondition, aes(y = median, x = cond, color = Module)) + geom_line() + geom_vline(xintercept = which(TTest.pvalue > 0.001), colour= "black", linetype = "dotted")
dev.off()

# Remove conditions with high score conditions in the leaf module
perconditionLeaf = percondition[!(percondition$cond %in% as.numeric(CondRankLeaf[1,1:10])), ]
ggplot(perconditionLeaf, aes(y = median, x = cond, color = Module)) + geom_line()

# Remove conditions with high score conditions in the pollen module
perconditionPollen = percondition[!(percondition$cond %in% as.numeric(CondRankPollen[1,1:10])), ]
ggplot(perconditionPollen, aes(y = median, x = cond, color = Module)) + geom_line()




# For anchorpoints
data.exp.j.Leaf = data.exp.compen[j,]
data.exp.jpar.Pollen = data.exp.compen[jpar,]

Exp.diff = data.exp.j.Leaf - data.exp.jpar.Pollen

jpeg("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/GenomeDominance/PerConditionAll.jpeg", width = 11, height = 7, units = "in", res = 200, pointsize=8)
myImagePlot(Exp.diff)
dev.off()

myImagePlot(Exp.diff[, -c(unlist(c(CondRankLeaf[1,1:10], CondRankPollen[1,1:10])))])
Exp.logRatio = log(data.exp.j.Leaf/data.exp.jpar.Pollen)

jpeg("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/GenomeDominance/AnchExpLogRatioAll.jpeg", width = 11, height = 7, units = "in", res = 200, pointsize=8)
myImagePlot(Exp.logRatio)
dev.off()

# Remove top 10 conditions for both modules, resulting 114 conditions
Exp.logRatio.general = Exp.logRatio[, -c(c(unlist(c(CondRankLeaf[1,1:10], CondRankPollen[1,1:10]))))]
jpeg("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/GenomeDominance/AnchExpLogRatioGeneral.jpeg", width = 11, height = 7, units = "in", res = 200, pointsize=8)
myImagePlot(Exp.logRatio.general)
dev.off()

# Fraction of gene expression matrix in which genes in the leaf module have higher absolute expression values when removing top 10 (20) conditions
sum(Exp.logRatio.general > 0)/(nrow(Exp.logRatio.general) * ncol(Exp.logRatio.general))

# Fraction of gene expression matrix in which genes in the leaf module have higher absolute expression values
sum(Exp.logRatio > 0)/(nrow(Exp.logRatio) * ncol(Exp.logRatio))




# -------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Root VS Leaf
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------
memb.comun = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/AT3G19710_95.txt") # Leaf
memb.comun.para = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/AT1G50090_95.txt") # Root


memb.comun = ref[as.character(memb.comun[,1]),2]
memb.comun.para = ref[as.character(memb.comun.para[,1]),2] 

Inter1 = intersect(memb.comun.para, paralogous[paralogous$Ref %in% memb.comun, "Anch"])
Inter2 = intersect(memb.comun, paralogous[paralogous$Ref %in% memb.comun.para, "Anch"])
aa = paralogous[paralogous$Ref %in% Inter2, ]
bb = aa[aa$Anch %in% Inter1,]
colnames(bb)[1] = "References"
Inter_name = merge(ref, bb, by = "References")
colnames(Inter_name)[c(1,3)] = c("Ref", "References")
Inter_name = merge(ref, Inter_name, by = "References")
	
j = Inter_name$Ref # Leaf
jpar = Inter_name$References # Root


# Condition scores for both modules
load("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Condition_score/Compendium_wise/development/Raw/C2858.RData")
CondRank1 = cond_scores
load("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Condition_score/Compendium_wise/development/Raw/C9603.RData")
CondRank2 = cond_scores
CondRank = rbind(c(1:length(CondRank1)), data.frame(t(CondRank1)), data.frame(t(CondRank2)))
CondRankRoot = CondRank[,order(CondRank[2, ], decreasing = T)]
CondRankLeaf = CondRank[,order(CondRank[3, ], decreasing = T)]
myImagePlot_sym(data.exp.scale[c(j, jpar), unlist(c(CondRankLeaf[1,1:10], CondRankRoot[1,1:10]))])


# For entire module
data.exp.Leaf = data.exp.compen[memb.comun,]
data.exp.Root = data.exp.compen[memb.comun.para,]
denPlot = data.frame("Exp" = as.numeric(unlist(data.exp.Root)), "Module" = rep("Root", length(as.numeric(unlist(data.exp.Root)))))
denPlot = rbind(denPlot, data.frame("Exp" = as.numeric(unlist(data.exp.Leaf)), "Module" = rep("Leaf", length(as.numeric(unlist(data.exp.Leaf))))))

jpeg ("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/GenomeDominance/DensityAll_90.jpeg", width = 11, height = 7, units = "in", res = 200, pointsize=8)
ggplot(denPlot, aes(x = Exp, fill = Module)) + geom_density(alpha = 0.3)
dev.off()

# KS test 
ks.test(as.numeric(unlist(data.exp.Root)), as.numeric(unlist(data.exp.Leaf)), alternative = c("greater"))

# Comparison per condition
percondition = data.frame("median" = apply(data.exp.Leaf, 2, median), "Module" = rep("Leaf", 134))
percondition = rbind(percondition, data.frame("median" = apply(data.exp.Root, 2, median), "Module" = rep("Root", 134)))
percondition["cond"] = rep(1:134,2)

jpeg ("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/GenomeDominance/PerConditionAll.jpeg", width = 11, height = 7, units = "in", res = 200, pointsize=8)
ggplot(percondition, aes(y = median, x = cond, color = Module)) + geom_line()
dev.off()


# T-test for comparing gene expression of each compendium
TTest.pvalue = foreach(i = 1:ncol(data.exp.Leaf), .combine = c) %do% {
	t.test(data.exp.Leaf[,i], data.exp.Root[,i], alternative = "greater")$p.value
}
TTest.pvalue = p.adjust(TTest.pvalue , method = "fdr")

jpeg ("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/GenomeDominance/PerConditionAllTtest.jpeg", width = 11, height = 7, units = "in", res = 200, pointsize=8)
ggplot(percondition, aes(y = median, x = cond, color = Module)) + geom_line() + geom_vline(xintercept = which(TTest.pvalue > 0.001), colour= "black", linetype = "dotted")
dev.off()

# Remove conditions with high score conditions in the leaf module
perconditionLeaf = percondition[!(percondition$cond %in% as.numeric(CondRankLeaf[1,1:10])), ]
ggplot(perconditionLeaf, aes(y = median, x = cond, color = Module)) + geom_line()

# Remove conditions with high score conditions in the Root module
perconditionRoot = percondition[!(percondition$cond %in% as.numeric(CondRankRoot[1,1:10])), ]
ggplot(perconditionRoot, aes(y = median, x = cond, color = Module)) + geom_line()




# For anchorpoints
data.exp.j.Leaf = data.exp.compen[j,]
data.exp.jpar.Root = data.exp.compen[jpar,]

Exp.diff = data.exp.j.Leaf - data.exp.jpar.Root

jpeg("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/GenomeDominance/PerConditionAll.jpeg", width = 11, height = 7, units = "in", res = 200, pointsize=8)
myImagePlot(Exp.diff)
dev.off()

myImagePlot(Exp.diff[, -c(unlist(c(CondRankLeaf[1,1:10], CondRankRoot[1,1:10])))])
Exp.logRatio = log(data.exp.j.Leaf/data.exp.jpar.Root)

jpeg("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/GenomeDominance/AnchExpLogRatioAll.jpeg", width = 11, height = 7, units = "in", res = 200, pointsize=8)
myImagePlot(Exp.logRatio)
dev.off()

# Remove top 10 conditions for both modules, resulting 114 conditions
Exp.logRatio.general = Exp.logRatio[, -c(c(unlist(c(CondRankLeaf[1,1:10], CondRankRoot[1,1:10]))))]
jpeg("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/GenomeDominance/AnchExpLogRatioGeneral.jpeg", width = 11, height = 7, units = "in", res = 200, pointsize=8)
myImagePlot(Exp.logRatio.general)
dev.off()

# Fraction of gene expression matrix in which genes in the leaf module have higher absolute expression values when removing top 10 (20) conditions
sum(Exp.logRatio.general > 0)/(nrow(Exp.logRatio.general) * ncol(Exp.logRatio.general))

# Fraction of gene expression matrix in which genes in the leaf module have higher absolute expression values
sum(Exp.logRatio > 0)/(nrow(Exp.logRatio) * ncol(Exp.logRatio))




