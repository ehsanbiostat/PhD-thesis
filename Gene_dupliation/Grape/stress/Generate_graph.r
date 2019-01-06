

library(igraph)
library("foreach")
library("doMC")
# registerDoMC(12)


	##------------------------------------------------------------------------------------------------------------------------------
	# Reading Gene scores from single file and integrate them to a matrix
	##------------------------------------------------------------------------------------------------------------------------------
	# change directory to location which Gene scores files are there
	setwd("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Permutation/stress/Gene_score/V2/") # For Grapevine 
	# Use ~/Gene_duplication/Scripts/General/order.pl to make a list for G_TotalSDi.txt files

	## Generate name of Gene score files by a Perl script and read it to an object
	list = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Scripts/General/order.gene.score.txt")
	list = as.vector(unlist(list[1:7179,]))

	## Read text file with read.table function
	g = lapply(list, function(x) mget(load(x)))

	total_gene_scores = matrix(NA, length(unlist(g[1])), length(g))

	for (i in 1:length(g)){
		print(i)
		total_gene_scores[, i] =  as.matrix(unlist(g[i]))
	}

	gene.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated CORNET2.0/Vitis_vinifera/Stress/gene.names.txt", header = F)
	
	rownames(total_gene_scores) = colnames(total_gene_scores) = gene.names[,1]
	
	# Build an undirected graph from adjacency matrix
	
	g = graph.adjacency(total_gene_scores, mode = c("undirected"), diag = F)


	save(g, file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/stress/Graph_V2.R") # For Grapevine


#}


