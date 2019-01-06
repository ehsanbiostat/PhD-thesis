
load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/development/Score99/Result/Total.RData")
A = sapply(g, unlist)
AA = sapply(A, length)
load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/development/Score99/Result.RData")
load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/development/Score99/G.score.RData")


## Grape
	library(igraph)
	data.grape = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated CORNET2.0/Vitis_vinifera/GE.mean.SD.txt", nrows = 19589, comment.char = "", header = T)
	gene.names <- read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated CORNET2.0/Vitis_vinifera/gene.names.txt")
	gene.names = as.matrix(gene.names)
	
	condition.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated CORNET2.0/condition_name.txt", header = T)
	NAME = names(table(condition.names))

	
	data.arab = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated CORNET2.0/Integrated_SD.txt", colClasses = "numeric", nrows = 19285, comment.char = "", header = T)
	data.arab = as.matrix(data.arab[, condition.names == "development"])
	
	
	
	k = 6
	
	## Grape modules
	load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Graph_V2.R")
	g.grape = g
	
	## Arabidopsis modules
	load(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/",NAME[k],"/Graph_V2.RData", sep = ""))
	
	orthologous = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Orthology/Grape-Arabidopsis/All.txt", header = T)
	gene.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated CORNET2.0/Vitis_vinifera/gene.names.txt")
	map = read.delim("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Orthology/Grape-Arabidopsis/Map.12x.plaza.txt")
	Ref = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated CORNET2.0/Reference_SD.txt", header = T)
	paralogous = read.table("~/Gene_duplication/Seed_genes/paralogous_total.txt", header = T)
	rownames(data.arab) = Ref[,1]


	

	
library(foreach)
library(doMC)
registerDoMC(50)	

Mod1 = G.score[G.score[, "Gene"] %in% unlist(A[which(AA > 1)[m]]) & (G.score[, "orthologous"] == Result[which(AA > 1)[m], "Gene"]),]
Mod2 = G.score[G.score[, "Gene"] %in% unlist(A[which(AA > 1)[m]]) & (G.score[, "orthologous"] == Result[which(AA > 1)[m], "Paralogous"]),]

Mod1 = Mod1[order(-Mod1[, "Normalized.orthologous.cluster"]), ]
Mod2 = Mod2[order(-Mod2[, "Normalized.orthologous.cluster"]), ]

Mod1 = Mod1[order(-Mod1[, "Gene"]), ]
Mod2 = Mod2[order(-Mod2[, "Gene"]), ]

QQ = (Mod1[, "orthologous.cluster"] + Mod2[, "orthologous.cluster"])/(Mod1[, "Gene.cluster.size"] + Mod2[, "Gene.cluster.size"]
	+ Mod1[, "orthologous.cluster.size"] + Mod2[, "orthologous.cluster.size"])

Mod1 = cbind(Mod1, "Score" = QQ)
Mod2 = cbind(Mod2, "Score" = QQ)
Mod1 = Mod1[order(-Mod1[, "Score"]), ]
Mod2 = Mod2[order(-Mod2[, "Score"]), ]

Total = list()

INTER.Tot1 = foreach (m = 1:length(which(AA > 1)), .combine = rbind) %dopar% {
	print(m)
	LEN = length(unlist(A[which(AA > 1)[m]]))
	Correlation = list()	
	for (n in 1:LEN)	{
		print(n)
		i = unlist(A[which(AA > 1)[m]])[n]
		
		memb.comun.grape = data.frame("name" = gene.names[neighborhood(g.grape, 1, i)[[1]],])
		memb.comun.grape.para = merge(memb.comun.grape, map, by = "name")
		## Remove duplicates 
		memb.comun.grape.para = memb.comun.grape.para[!(duplicated(memb.comun.grape.para$name) 
							| duplicated(memb.comun.grape.para$name, fromLast = TRUE)), ]
		memb.comun.grape.para = merge(memb.comun.grape.para, orthologous, by = "PLAZA.ID")
		
		j = Result[which(AA > 1)[m], "Gene"]
		memb.comun.Arab = data.frame("gene_id.y" = Ref[neighborhood(g, 1, j)[[1]],1])
		Orth1 = merge(memb.comun.Arab, memb.comun.grape.para, by = "gene_id.y")$name
		Orth1.arab = merge(memb.comun.Arab, memb.comun.grape.para, by = "gene_id.y")$gene_id.y
		gene.exp1.grape = data.grape[Orth1, ]
		gene.exp1.arab = data.arab[Orth1.arab, ]
		
		j = Result[which(AA > 1)[m], "Paralogous"]
		memb.comun.Arab = data.frame("gene_id.y" = Ref[neighborhood(g, 1, j)[[1]],1])
		Orth2 = merge(memb.comun.Arab, memb.comun.grape.para, by = "gene_id.y")$name
		Orth2.arab = merge(memb.comun.Arab, memb.comun.grape.para, by = "gene_id.y")$gene_id.y
		gene.exp2.grape = data.grape[Orth2, ]
		gene.exp2.arab = data.arab[Orth2.arab, ]
		
		cor1.arab = cor(t(gene.exp1.arab))
		cor2.arab = cor(t(gene.exp2.arab))
		cor1.grape = cor(t(gene.exp1.grape))
		cor2.grape = cor(t(gene.exp2.grape))
		
		COrrelation1 = COrrelation2 = c()
		
		for (i in 1:nrow(cor1.arab)) {
			COrrelation1[i] = cor(cbind(cor1.arab[i,], cor1.grape[i,]))[1,2]
		}
		
		for (i in 1:nrow(cor2.arab)) {
			COrrelation2[i] = cor(cbind(cor2.arab[i,], cor2.grape[i,]))[1,2]
		}
		
		Correlation[n] = list(list(COrrelation1, COrrelation2))
		
	}	
	
		
	}
	max(INTER)
}	

for (i in 1:LEN){
	print(ks.test(unlist(Correlation[i][[1]][1]), unlist(Correlation[i][[1]][2]), 
		exact = T)$p)
}

save(INTER.Tot1, file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/development/Score99/Result/Intersect_Orth.RData")
	

	
	go.terms = merge(go.names, memb.comun.grape, by = "name")[, "go"]
	
	freq = table(matrix(go.terms)) # Frequency of GO among
	 pvalue = odd = global = loc = c()

	 for (i in 1:length(freq)){
		 z = matrix(c(length(go.terms) - freq[i], freq[i]), 2)
		 x = matrix(table(go.names[, "go"] == names(freq[i])))
		 x = cbind(x ,z) 
		 pvalue[i] = fisher.test(x, alternative = "greater")$p
		 global[i] = x[2,1]
		 loc[i] = x[2,2]
		 odd[i] = (x[1,1] * x[2,2])/(x[1,2] * x[2,1])
	 }
	
	 cat("First modules finished\n")
	 names(pvalue) = names(odd) = names(global) = names(loc) = names(freq)
	 f = data.frame(pvalue, odd, global, loc)
	 f["FDR"] = p.adjust(f$pvalue, method = c("fdr"))
	 f = f[order(f$FDR), ]
	 top = c(rownames(f[f$FDR < 0.05,]))
	 f = f[f$FDR < 0.05,]
	 f["GO"] = rownames(f)
	 RRf = merge(RR, f, by = "GO")
	 RRf = RRf[order(RRf$FDR, decreasing = F), ]
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	

	
	