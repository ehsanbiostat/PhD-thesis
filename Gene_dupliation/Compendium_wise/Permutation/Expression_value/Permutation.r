
library(isa2)
load("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Input_data/Input_data.RData")

smartSeed <- matrix(0, nrow = dim(data.norm$Ec), ncol = 1) # every column is a seedvector
smartSeed[1, 1] <- 1

#
# Perform the ISA steps
#

cond_scores_general <- data.norm$Er %*% smartSeed
cond_scores_matrix = matrix(NA, dim(data.norm$Er), 1000)
gene_scores_matrix = matrix(NA, dim(data.norm$Ec), 1000)

for (i in 1:1000){
	print(i)
# 1) Calculate the condition scores for the seed gene ==> in case of only one seed gene this boils down to a vector of its gene expression values
	cond_scores <- as.matrix(sample(cond_scores_general))
	cond_scores <- apply(cond_scores, 2, function(x) x/sqrt(sum(x^2))) # normalize by the length of the vector ==> vector of length 1

# 2) Get the gene score vector
	
	gene_scores <- apply(gene_scores, 2, function(x) x/sqrt(sum(x^2))) # normalize by the length of the vector ==> vector of length 1

# 3) Recalculate the condition scores
	cond_scores <- data.norm$Er %*% gene_scores
	cond_scores <- apply(cond_scores, 2, function(x) x/sqrt(sum(x^2))) # normalize by the length of the vector ==> vector of length 1

	cond_scores_matrix[, i] = cond_scores
	gene_scores_matrix[, i] = gene_scores
}


thresh_gene_score = apply(gene_scores_matrix, 1, quantile,probs = c(0.95))
thresh_cond_score = apply(cond_scores_matrix, 1, quantile,probs = c(0.95))
SA_gene_score = read.table("/group/biocomp/users/ehsab/Gene_duplication/Results/SA/Gene_score/Total/G_Total.SD1.txt", header = T)
SA_cond_score = read.table("/group/biocomp/users/ehsab/Gene_duplication/Results/SA/Condition_score/Total/C_Total.SD1.txt", header = T)

gene_score_module = ifelse(SA_gene_score > thresh_gene_score, 1, 0)
cond_score_module = ifelse(SA_cond_score > thresh_cond_score, 1, 0)

colnames(cond_score_module) = colnames(gene_score_module) = gene.names[1]
write.table(gene_score_module, file = "/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Permutation/Scores/G_Total.SD1.txt", row.names = F, col.names = T, sep = "\t")
write.table(cond_score_module, file = "/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Condition_score/Permutation/Scores/C_Total.SD1.txt", row.names = F, col.names = T, sep = "\t")

write.table(thresh_gene_score, file = "/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Permutation/Threshold/G_Thresh.SD1.txt", row.names = F, col.names = F, sep = "\t")
write.table(thresh_cond_score, file = "/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Condition_score/Permutation/Threshold/C_Thresh.SD1.txt", row.names = F, col.names = F, sep = "\t")






















