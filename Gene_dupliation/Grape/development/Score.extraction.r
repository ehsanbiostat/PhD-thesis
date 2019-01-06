library(foreach)
library(doMC)
library(preprocessCore)
registerDoMC(10)

load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/development/Score99/Result/Total1.RData")

# ICC results
load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/development/ICC.RData")

# Select paralogous modules and their ancesteral modules from grapvine
# A pair of paralogous modules should have a same grapevine modules which has a significant enrichment score
A = sapply(g, unlist)
AA = sapply(A, length)

# Load paralogous modules from Arabidopsis data
load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/development/Score99/Result.RData")

# Load significant grape modules modules with arabidopsis modules
load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/development/Score99/G.score.RData")

load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/development/Score99/Result/Intersect_Orth.RData")

# Repeat row function
rep.row <- function(x,n){
   matrix(rep(x,each=n),nrow=n)
}


## Grape
	library(igraph)
	condition.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/condition_name.txt", header = T)
	NAME = names(table(condition.names))

	k = 6
	
	## Grape modules
	load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Graph_V2.R")
	g.grape = g
	
	## Arabidopsis modules
	load(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/",NAME[k],"/Graph_V2.RData", sep = ""))
	
	orthologous = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Orthology/Grape-Arabidopsis/All.txt", header = T)
	gene.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Vitis_vinifera/gene.names.txt")
	map = read.delim("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Orthology/Grape-Arabidopsis/Map.12x.plaza.txt")
	Ref = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Reference_SD.txt", header = T)
	paralogous = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Seed_genes/paralogous_total.txt", header = T)


# Quantile normalization of all orthology scores
TT.scale = normalize.quantiles(G.score[, -c(1:9)])
colnames(TT.scale) = colnames(G.score)[-c(1:9)]
G.score.copy = G.score
G.score = cbind(TT.scale, "orthologous" = G.score[, "orthologous"], 
	"Gene" = G.score[, "Gene"], "orthologous.cluster" = G.score[, "orthologous.cluster"])


# Quantile normalization of all paralogy scores
TT.scale.Para = normalize.quantiles(Result[, -c(1:12)])
colnames(TT.scale.Para) = colnames(Result)[-c(1:12)]
Result.copy = Result
Result = cbind(TT.scale.Para, "Paralogous" = Result[, "Paralogous"], 
	"Gene" = Result[, "Gene"], "Paralogous.cluster" = Result[, "Paralogous.cluster"])




TT = foreach (m = 1:length(which(AA > 1)), .combine = rbind) %dopar% { 
	print(m)
	
	# Extract all potential pre-duplicated modules for each paralogous module
	i = unlist(A[which(AA > 1)[m]])
	# Extract Arabidopsis modules which shared potential pre-duplicated modules
	j = Result[which(AA > 1)[m], "Gene"]
	k = Result[which(AA > 1)[m], "Paralogous"]


	# Select orthologous scores
	graph.line = data.frame(G.score[G.score[, "Gene"] %in% i & (G.score[, "orthologous"] == j | G.score[, "orthologous"] == k), c("Normalized.orthologous.cluster"
		, "Normalized.orthologous.cluster.grape" ,"Normalized.orthologous.cluster.Arab" , "Normalized.orthologous.whitin.grape" 
		, "Normalized.orthologous.whitin.Arab" , "orthologous", "Gene", "orthologous.cluster")])

	graph.line = graph.line[order(graph.line[,"Gene"]),]
	# Select paralogous scores
	Arab.scores = data.frame(rep.row(Result[which(AA > 1)[m], c(1:6, 9)], nrow(graph.line)/2))
	colnames(Arab.scores) = colnames(Result)[c(1:6, 9)]

	# ICC score
	CORR = foreach (p = 1:length(Correlation_total[[m]]), .combine = rbind) %do% {
		sapply(Correlation_total[[m]][[p]], median)
	}
	colnames(CORR) = c("Cor1", "Cor2", "Cor1.mutual", "Cor2.mutual")

	# Generate row for each Arabidopsis-Grapevine module combination and all
	# related scores
	cbind(graph.line[seq.int(1, nrow(graph.line), 2), ], graph.line[seq.int(2, nrow(graph.line), 2), ], data.frame(INTER.Tot1[m]), data.frame(CORR), Arab.scores)
}


save(TT, file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/development/Combined_Score.RData")

