
library(igraph)
library(foreach)
library(doMC)
registerDoMC(1)

## Read condition names
condition.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/condition_name.txt", header = T)
NAME = names(table(condition.names))
k = 1
NAME.spe = NAME[k]

paralogous = read.table("~/Gene_duplication/Seed_genes/paralogous_NoFilter.txt", header = T)
PARA = unique(paralogous$Ref)

Score = read.delim(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Local/Compendium_wise/",NAME[k],"/Score99/Scores_V2/Anch_Score.txt", sep = ""))
Score = Score[Score$Jaccard != 0, ]
PARA = unique(c(Score$Gene, Score$Paralogous))

# --------------------------------------------------------------------------------------------------------------------------------------------------
## Condition scores
# --------------------------------------------------------------------------------------------------------------------------------------------------
list = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Scripts/General/order.condition_score.txt")
list = as.vector(unlist(list))
## Read text or RData file with read.table or load function
setwd(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Condition_score/Compendium_wise/",NAME.spe,"/Permutation/V2", sep = ""))
# C = lapply(list, read.table, sep= "\t" , header=T)
C = lapply(list, function(x) mget(load(x)))
total_cond_scores = matrix(NA, length(unlist(C[1])), 19285)

for (i in 1:19285){
	print(i)
	total_cond_scores[, i] = as.matrix(unlist(C[i]))
}
# --------------------------------------------------------------------------------------------------------------------------------------------------




##------------------------------------------------------------------------------------------------------------------------------
# Reading Gene scores from single file and integrate them to a matrix
##------------------------------------------------------------------------------------------------------------------------------
# change directory to location which Gene scores files are there
setwd(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/",NAME.spe,"/Permutation/V2", sep = ""))

## Generate name of Gene score files by a Perl script and read it to an object
list = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Scripts/General/order.gene.score.txt")
list = as.vector(unlist(list)[1:19285])

## Read text file with read.table function
g = lapply(list, function(x) mget(load(x)))

total_gene_scores = matrix(NA, length(unlist(g[1])), length(g))

for (i in 1:length(g)){
	print(i)
	total_gene_scores[, i] =  as.matrix(unlist(g[i]))
}
##------------------------------------------------------------------------------------------------------------------------------





# --------------------------------------------------------------------------------------------------------------------------------------------------
# Merging gene and condition score for all genes
# We need to read condition and gene score matrices, then merge them and apply a distance function
# --------------------------------------------------------------------------------------------------------------------------------------------------

total_cond_scores = t(total_cond_scores) 
SS = apply(total_cond_scores, 1, sum)
TT = matrix(0, 19285, 255)


for (i in 1:sum(SS > 0)) {
	print(i)
	GG = matrix(0, 19285, 255)
	GG[which(total_gene_scores[which(SS > 0)[i],] > 0), which(total_cond_scores[which(SS > 0)[i],] > 0)] = 1
	TT = TT + GG
}

save(TT, file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Merge_module/abiotic/TT.RData")
library(gplots)
library(lattice)
jpeg("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Merge_module/abiotic/Matrix.jpg", width=8, height=5, units="in", res=500, pointsize = 3)
heatmap.2(TT1, Rowv=FALSE, Colv=FALSE, dendrogram="none", tracecol="#303030", trace="none", key = FALSE)
dev.off()

jpeg("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Merge_module/abiotic/Matrix2.jpg", width=8, height=5, units="in", res=500, pointsize = 3)
levelplot(TT)
dev.off()
