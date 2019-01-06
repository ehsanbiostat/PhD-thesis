
load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/abiotic/Score99/Result/Total1.RData")
A = sapply(g, unlist)
AA = sapply(A, length)
load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/abiotic/Score99/G.score.RData")
load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/abiotic/Score99/Result.RData")


## Grape
	library(igraph)
	data.grape = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated CORNET2.0/Vitis_vinifera/Stress/Exp.matrix.txt", nrows = 7179, comment.char = "", header = T)
	gene.names <- rownames(data.grape)
	gene.names = as.matrix(gene.names)
	
	condition.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated CORNET2.0/condition_name.txt", header = T)
	NAME = names(table(condition.names))

	## Arabidopsis
	data.arab = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated CORNET2.0/Integrated_SD.txt", colClasses = "numeric", nrows = 19285, comment.char = "", header = T)
	data.arab = as.matrix(data.arab[, condition.names == "abiotic"])
	
	
	
	k = 1
	
	## Grape modules
	load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/stress/Graph_V2.R")
	g.grape = g
	
	## Arabidopsis modules
	load(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/",NAME[k],"/Graph_V2.RData", sep = ""))
	
	orthologous = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Orthology/Grape-Arabidopsis/All.txt", header = T)
	map = read.delim("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Orthology/Grape-Arabidopsis/Map.12x.plaza.txt")
	Ref = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated CORNET2.0/Reference_SD.txt", header = T)
	paralogous = read.table("~/Gene_duplication/Seed_genes/paralogous_total.txt", header = T)
	rownames(data.arab) = Ref[,1]


	

	
library(foreach)
library(doMC)
registerDoMC(16)	

Correlation_total = foreach (m = 1:length(which(AA > 1))) %dopar% {
	print(m)
	LEN = length(unlist(A[which(AA > 1)[m]]))
	Correlation = list()	
	for (n in 1:LEN)	{
		
		i = unlist(A[which(AA > 1)[m]])[n]
		
		memb.comun.grape.para = data.frame("PLAZA.ID" = gene.names[neighborhood(g.grape, 1, i)[[1]],])
		## Remove duplicates 
		memb.comun.grape.para = data.frame("PLAZA.ID" = memb.comun.grape.para[!(duplicated(memb.comun.grape.para$PLAZA.ID) 
							| duplicated(memb.comun.grape.para$PLAZA.ID, fromLast = TRUE)), ])
		memb.comun.grape.para = merge(memb.comun.grape.para, orthologous, by = "PLAZA.ID")
		
		# First pair arabidopsis module to vitis module
		j = Result[which(AA > 1)[m], "Gene"]
		memb.comun.Arab = data.frame("gene_id.y" = Ref[neighborhood(g, 1, j)[[1]],1])
		# Orthologous genes between each module
		Orth1 = merge(memb.comun.Arab, memb.comun.grape.para, by = "gene_id.y")$PLAZA.ID
		Orth1.arab = merge(memb.comun.Arab, memb.comun.grape.para, by = "gene_id.y")$gene_id.y
		# Gene expression profile for orthologous genes
		gene.exp1.grape = data.grape[Orth1, ]
		gene.exp1.arab = data.arab[Orth1.arab, ]
		
		# Second pair arabidopsis module to vitis module
		j = Result[which(AA > 1)[m], "Paralogous"]
		memb.comun.Arab = data.frame("gene_id.y" = Ref[neighborhood(g, 1, j)[[1]],1])
		# Orthologous genes between each module
		Orth2 = merge(memb.comun.Arab, memb.comun.grape.para, by = "gene_id.y")$PLAZA.ID
		Orth2.arab = merge(memb.comun.Arab, memb.comun.grape.para, by = "gene_id.y")$gene_id.y
		# Gene expression profile for orthologous genes
		gene.exp2.grape = data.grape[Orth2, ]
		gene.exp2.arab = data.arab[Orth2.arab, ]
		
		# Whitin correlation between paralogous genes in each module
		cor1.arab = cor(t(gene.exp1.arab))
		cor2.arab = cor(t(gene.exp2.arab))
		cor1.grape = cor(t(gene.exp1.grape))
		cor2.grape = cor(t(gene.exp2.grape))
		
		COrrelation1 = COrrelation2 = c()
		
		# Between orthologous gene correlation 
		for (i in 1:nrow(cor1.arab)) {
			COrrelation1[i] = cor(cbind(cor1.arab[i,], cor1.grape[i,]))[1,2]
		}
		
		for (i in 1:nrow(cor2.arab)) {
			COrrelation2[i] = cor(cbind(cor2.arab[i,], cor2.grape[i,]))[1,2]
		}
		
		Correlation[n] = list(list(COrrelation1, COrrelation2))
		
	}	
	Correlation

}

save(Correlation_total, file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/abiotic/ICC.RData")

# for (i in 1:LEN){
# 	print(ks.test(unlist(Correlation[i][[1]][1]), unlist(Correlation[i][[1]][2]), 
# 		exact = T)$p)
# }

# CORR = foreach (i = 1:length(which(AA > 1)), .combine = rbind) %do% {
# 	print(i)
# 	sapply(Correlation_total[[i]][[1]], median)
# }
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	

	
	
