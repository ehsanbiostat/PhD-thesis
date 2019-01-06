
# List of significant shared ancesteral ortholog
load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/development/Score99/Result/Total.RData")

# Select paralogous modules and their ancesteral modules from grapvine
# A pair of paralogous modules should have a same grapevine modules which has a significant enrichment score
A = sapply(g, unlist)
AA = sapply(A, length)

# Load paralogous modules from Arabidopsis data
load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/development/Score99/Result.RData")

# Load significant grape modules modules with arabidopsis modules
load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/development/Score99/G.score.RData")


## Grape
	library(igraph)
	condition.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated CORNET2.0/condition_name.txt", header = T)
	NAME = names(table(condition.names))

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


	
library(foreach)
library(doMC)
registerDoMC(50)	
# INTER.Tot = c()
INTER.Tot1 = foreach (m = 1:length(which(AA > 1)), .combine = rbind) %dopar% {
	print(m)
	LEN = length(unlist(A[which(AA > 1)[m]]))
	INTER = c()
	for (n in 1:LEN)	{
		
		i = unlist(A[which(AA > 1)[m]])[n]

		memb.comun.grape = data.frame("name" = gene.names[neighborhood(g.grape, 1, i, mode = c("out"))[[1]],])
		memb.comun.grape.para = merge(memb.comun.grape, map, by = "name")
		## Remove duplicates 
		memb.comun.grape.para = memb.comun.grape.para[!(duplicated(memb.comun.grape.para$name) 
							| duplicated(memb.comun.grape.para$name, fromLast = TRUE)), ]
		memb.comun.grape.para = merge(memb.comun.grape.para, orthologous, by = "PLAZA.ID")
		
		j = Result[which(AA > 1)[m], "Gene"]
		memb.comun.Arab = data.frame("gene_id.y" = Ref[neighborhood(g, 1, j)[[1]],1])
		Orth1 = merge(memb.comun.Arab, memb.comun.grape.para, by = "gene_id.y")$name
		
		j = Result[which(AA > 1)[m], "Paralogous"]
		memb.comun.Arab = data.frame("gene_id.y" = Ref[neighborhood(g, 1, j)[[1]],1])
		Orth2 = merge(memb.comun.Arab, memb.comun.grape.para, by = "gene_id.y")$name
		# INTER[n] = length(intersect(Orth2, Orth1))
		
	# Arabidopsis
		
		i = Result[which(AA > 1)[m], "Gene"]

		memb.comun = data.frame("Ref" = neighborhood(g, 1, i, mode = c("out"))[[1]])
		memb.comun.para = merge(memb.comun, paralogous, by = "Ref")

		j = Result[which(AA > 1)[m], "Paralogous"]
		Para = Ref$Gene[intersect(neighborhood(g, 1, j, mode = c("out"))[[1]], memb.comun.para$Anch)]
		
		Anch = data.frame("Anch" = neighborhood(g, 1, j, mode = c("out"))[[1]])
		Anch = merge(Anch, memb.comun.para, by = "Anch")
		
		INTER[n] = length(intersect(Orth2, Para))
		# INTER.Tot[m] = 
	}
	max(INTER)
}	
	
save(INTER.Tot1, file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/development/Score99/Result/Intersect_Orth.RData")
	
	
## Line plot
foreach (m = 1:20) %dopar% { 
	print(m)
	i = unlist(A[which(AA > 1)[m]])
	j = Result[which(AA > 1)[m], "Gene"]
	k = Result[which(AA > 1)[m], "Paralogous"]
	graph.line = data.frame(G.score[G.score[, "Gene"] %in% i & (G.score[, "orthologous"] == j | G.score[, "orthologous"] == k),
		c("Normalized.orthologous.cluster", "orthologous", "Gene")])

	graph.line["Group"] = Ref[graph.line$orthologous, "Gene"]
	graph.line["X"] = gene.names[graph.line$Gene, ]
	jpeg(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Grape/Multiple_othologous_modules/",m,".jpg", sep = "")
		, width=8, height=5, units="in", res=500, pointsize = 1)
	ggplot(data=graph.line,
		aes(x=X, y=Normalized.orthologous.cluster, group=Group, colour=Group)) + geom_line() + geom_point() + geom_text(size = 9)
	dev.off()
}
	
	
	# Arabidopsis
	CondRank1 = foreach(i = 1:nrow(Anch), .combine = rbind) %do%{
		rank(total_cond_scores[, Anch$Anch[i]])
	}


	CondRank2 = foreach(i = 1:nrow(Anch), .combine = rbind) %do%{
		rank(total_cond_scores[, Anch$Ref[i]])
	}
	
	
	CondRank = rbind(c(1:dim(total_cond_scores)[1]), CondRank1, CondRank2)
	CondRank0 = CondRank[,order(CondRank[2, ], decreasing = T)]
	
	CondRank = rbind(c(1:dim(total_cond_scores)[1]), CondRank2, CondRank1)
	CondRank01 = CondRank[,order(CondRank[2, ], decreasing = T)]
	
	myImagePlot_sym(data.exp.scale[c(Anch$Anch, Anch$Ref), CondRank01[1,]])
	myImagePlot_sym(data.exp.scale[c(Anch$Anch, Anch$Ref), CondRank0[1,]])
	myImagePlot_sym(data.exp.scale[c(Anch$Anch, Anch$Ref),c(CondRank01[1,1:20],CondRank0[1,1:20])])
	
	
	# Grape
require(gridExtra)
library(ggplot2)
library(igraph)
library(foreach)
library(doMC)
# registerDoMC(4)

## Read condition names
condition.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated CORNET2.0/condition_name.txt", header = T)
NAME = names(table(condition.names))
NAME.spe = NAME[6]

## Condition scores
list = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Scripts/General/order.condition_score.Grape.txt")
list = as.vector(unlist(list))
## Read text or RData file with read.table or load function
setwd("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Raw/Condition_score/")
# C = lapply(list, read.table, sep= "\t" , header=T)
C = lapply(list, read.table, header = T)
total_cond_scores = matrix(NA, length(unlist(C[1])), 19589)
data.grape.scale = apply(data.grape,1,scale)
data.grape.scale = t(data.grape.scale)
	
for (i in 1:19589){
	# print(i)
	total_cond_scores[, i] = as.matrix(unlist(C[i]))
}

	colnames(total_cond_scores) = gene.names[,1]
	
	CondRank.grape = foreach(i = 1:nrow(memb.comun.grape), .combine = rbind) %do%{
		rank(total_cond_scores[, memb.comun.grape[i,1]])
	}
	
	CondRank.grape = rbind(c(1:dim(total_cond_scores)[1]), CondRank.grape)
	CondRank.grape = CondRank.grape[,order(CondRank.grape[2, ], decreasing = T)]
	myImagePlot_sym(data.grape.scale[memb.comun.grape[ ,1], CondRank.grape[1,]])
	
	quantile(cor(t(data.grape.scale[memb.comun.grape[ ,1], CondRank.grape[1,]])))
	
	
	
	

	venn.diagram(X, category = c("A", "B", "C", "D", "E"), filename = "/home/ehsab/Gene_duplication/test1.png", imagetype = "png")
	
	
	
	
	
	
