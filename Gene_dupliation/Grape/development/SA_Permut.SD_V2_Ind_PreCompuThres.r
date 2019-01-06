

library("isa2")
library(foreach)
library(doMC)
registerDoMC(50)
# expressionMat is the gene expression matrix

foreach(j = 1:23629) %dopar% {
	print(j)
	# load thresholds
	thresh_gene_score = readRDS(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Permutation/Gene_score/V2/Threshold/T",j,".rds", sep = ""))
	# thresh_cond_score = readRDS(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Permutation/Condition_score/V2/Threshold/T",j,".rds", sep = ""))
	# Reading the actual values of gene and condition scores 

	gene_scores = read.table(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Raw/Gene_score/G_Total.SD",j,".txt", sep = ""), header = T)
	# cond_scores = read.table(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Raw/Condition_score/C_Total.SD",j,".txt", sep = ""), header = T)
	
	# Comparing actual values with threshold and assigning 0 and 1 codes
	gene_score_module = ifelse(gene_scores > thresh_gene_score[, "99%"], 1, 0)
	# cond_score_module = ifelse(cond_scores > thresh_cond_score[, "99%"], 1, 0)

	# colnames(cond_score_module) = colnames(gene_score_module) = gene.names[j]
	# colnames(gene_score_module) = gene.names[j]
	
	save(gene_score_module, file = paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Permutation/Gene_score/99/G",j,".RData", sep = ""))
	# save(cond_score_module, file = paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Condition_score/Compendium_wise/",NAME[k],"/Permutation/V2/90/C",j,".RData", sep = ""))
}
	
