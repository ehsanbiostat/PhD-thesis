

library("isa2")
library("foreach")
library("doMC")
registerDoMC(8)
# expressionMat is the gene expression matrix
data.exp = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated CORNET2.0/Integrated_SD.txt", colClasses = "numeric", nrows = 19285, comment.char = "", header = T)
condition.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated CORNET2.0/condition_name.txt", header = T)
NAME = names(table(condition.names))

gene.names <- read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated CORNET2.0/gene.names.low.mod.SD.txt")
gene.names = as.matrix(gene.names)

data = as.matrix(data.exp[, condition.names == NAME[1]])
# Normalize gene expression in gene-wise and condition-wise
data.norm <- isa.normalize(data)


foreach (j = 1:nrow(data)) %dopar% {
	print(j)
	# Select one compendium
	
	#
	# Prepare seed set 
	#
	## Reading 



	smartSeed <- matrix(0, nrow = dim(data.norm$Ec), ncol = 1) # every column is a seedvector
	smartSeed[j, 1] <- 1

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
	# 1) Calculate the condition scores for the seed gene ==> in case of only one seed gene this boils down to a vector of its gene expression values
		# Random sampling from original gene expression values of gene of interest
		cond_scores <- as.matrix(sample(cond_scores_general))
		cond_scores <- apply(cond_scores, 2, function(x) x/sqrt(sum(x^2))) # normalize by the length of the vector ==> vector of length 1

	# 2) Get the gene score vector
		gene_scores <- data.norm$Ec %*% cond_scores
		gene_scores <- apply(gene_scores, 2, function(x) x/sqrt(sum(x^2))) # normalize by the length of the vector ==> vector of length 1

	# 3) Recalculate the condition scores
		cond_scores <- data.norm$Er %*% gene_scores
		cond_scores <- apply(cond_scores, 2, function(x) x/sqrt(sum(x^2))) # normalize by the length of the vector ==> vector of length 1

		cond_scores_matrix[, i] = cond_scores
		gene_scores_matrix[, i] = gene_scores
	}

	# Compute threshold values
	thresh_gene_score = apply(gene_scores_matrix, 1, quantile,probs = c(0.99))
	thresh_cond_score = apply(cond_scores_matrix, 1, quantile,probs = c(0.99))
	# Reading the actual values of gene and condition scores 
	load(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/",as.name(NAME[1]),"/Raw/G",j,".RData", sep = ""))
	load(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Condition_score/Compendium_wise/",as.name(NAME[1]),"/Raw/C",j,".RData", sep = ""))
	
	
	# Comparing actual values with threshold and assigning 0 and 1 codes
	gene_score_module = ifelse(gene_scores > thresh_gene_score, 1, 0)
	cond_score_module = ifelse(cond_scores > thresh_cond_score, 1, 0)

	colnames(cond_score_module) = colnames(gene_score_module) = gene.names[j]
	
	save(gene_score_module, file = paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/",as.name(NAME[1]),"/Permutation/G",j,".RData", sep = ""))
	save(cond_score_module, file = paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Condition_score/Compendium_wise/",as.name(NAME[1]),"/Permutation/C",j,".RData", sep = ""))
	
}

