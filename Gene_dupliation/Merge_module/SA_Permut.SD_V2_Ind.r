

library("isa2")
library("foreach")
library("doMC")
registerDoMC(6)

# expressionMat is the gene expression matrix
condition.names = read.delim("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/condition_name.txt")
NAME = names(table(condition.names))
k = 1
gene.names <- read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/gene.names.low.mod.SD.txt")
gene.names = as.matrix(gene.names)

load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Local/Compendium_wise/abiotic/data.norm.RData")
Seed = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Merge_module/abiotic/Seed_100_20_80.txt")

foreach (j = 1:length(table(Seed$Cluster))) %dopar% {
	# Select one compendium
	
	#
	# Prepare seed set 
	#
	## Reading 



	smartSeed <- matrix(0, nrow = dim(data.norm$Ec), ncol = 1) # every column is a seedvector
	smartSeed[Seed[Seed$Cluster == j,"PARA"], 1] <- 1

	#
	# Perform the ISA steps
	#

	# Extract the gene expression value for seed genes and keep it intact
	cond_scores_general <- data.norm$Er %*% smartSeed
	
	# Generate two matrix in order to save the gene and condition score resulting from each permutation step
	cond_scores_matrix = matrix(NA, dim(data.norm$Er), 1000)
	gene_scores_matrix = matrix(NA, dim(data.norm$Ec), 1000)
	
	# Forloop for running 1000 permutaions
	for (i in 1:1000){
		print(i)
		
		load(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Local/Compendium_wise/",NAME[k],"/Randomization/95/",i,".RData", sep = ""))
	# 1) Calculate the condition scores for the seed gene ==> in case of only one seed gene this boils down to a vector of its gene expression values
		# Random sampling from original gene expression values of gene of interest
		cond_scores <- data.norm$Er %*% smartSeed
		cond_scores <- apply(cond_scores, 2, function(x) x/sqrt(sum(x^2))) # normalize by the length of the vector ==> vector of length 1
	
	
	# 2) Get the gene score vector
		gene_scores <- sampled.Ec %*% cond_scores
		gene_scores <- apply(gene_scores, 2, function(x) x/sqrt(sum(x^2))) # normalize by the length of the vector ==> vector of length 1
	
	
	# 3) Recalculate the condition scores
		cond_scores <- sampled.Er %*% gene_scores
		cond_scores <- apply(cond_scores, 2, function(x) x/sqrt(sum(x^2))) # normalize by the length of the vector ==> vector of length 1

		cond_scores_matrix[, i] = cond_scores
		gene_scores_matrix[, i] = gene_scores
	}

	# Compute threshold values
	thresh_gene_score = apply(gene_scores_matrix, 1, quantile,probs = c(0.95))
	thresh_cond_score = apply(cond_scores_matrix, 1, quantile,probs = c(0.95))
	# Reading the actual values of gene and condition scores 
	load(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Merge/Gene_score/Compendium_wise/",NAME[k],"/Raw/G_100_20_80",j,".RData", sep = ""))
	load(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Merge/Condition_score/Compendium_wise/",NAME[k],"/Raw/C_100_20_80",j,".RData", sep = ""))
	
	
	# Comparing actual values with threshold and assigning 0 and 1 codes
	gene_score_module = ifelse(gene_scores > thresh_gene_score, 1, 0)
	cond_score_module = ifelse(cond_scores > thresh_cond_score, 1, 0)

	colnames(cond_score_module) = colnames(gene_score_module) = gene.names[j]
	
	save(gene_score_module, file = paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Merge/Gene_score/Compendium_wise/",NAME[k],"/Permutation/V2/G_100_20_80",j,".RData", sep = ""))
	save(cond_score_module, file = paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Merge/Condition_score/Compendium_wise/",NAME[k],"/Permutation/V2/C_100_20_80",j,".RData", sep = ""))
}	



