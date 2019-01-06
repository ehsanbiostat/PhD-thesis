

library(igraph)
library("foreach")
library("doMC")
# registerDoMC(12)


	##------------------------------------------------------------------------------------------------------------------------------
	# Reading Gene scores from single file and integrate them to a matrix
	##------------------------------------------------------------------------------------------------------------------------------
	# change directory to location which Gene scores files are there
	# setwd(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/",as.name(NAME[j]),"/Permutation/V3", sep = ""))
	setwd("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Permutation/Gene_score/99/") # For Grapevine 
	# Use ~/Gene_duplication/Scripts/General/order.pl to make a list for G_TotalSDi.txt files

	## Generate name of Gene score files by a Perl script and read it to an object
	list = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Scripts/General/order.gene_score.Grape.txt")
	list = as.vector(unlist(list))

	## Read text file with read.table function
	g = lapply(list, function(x) mget(load(x)))

	total_gene_scores = matrix(NA, length(unlist(g[1])), length(g))

	for (i in 1:length(g)){
		print(i)
		total_gene_scores[, i] =  as.matrix(unlist(g[i]))
	}

	# gene.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated CORNET2.0/Vitis_vinifera/gene.names.txt", header = F)
	# rownames(total_gene_scores) = colnames(total_gene_scores) = gene.names[,1]
	
	# Build an undirected graph from adjacency matrix
	
	g = graph.adjacency(t(total_gene_scores), mode = c("directed"), diag = F)

	# save(g, file = paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/",NAME[j],"/Graph_V3.RData", sep = ""))
	
	
	save(g, file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Graph_V2.RData") # For Grapevine


#}


