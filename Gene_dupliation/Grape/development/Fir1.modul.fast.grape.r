

	 
library(fields)
library(igraph)
library(foreach)
library(doMC)
# registerDoMC(14)

	condition.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/condition_name.txt", header = T)
	NAME = names(table(condition.names))

	k = 6
	
	## Grape modules
	load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Graph_V2.RData")
	g.grape = g
	
	## Arabidopsis modules
	load(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/",NAME[k],"/Graph_V2.RData", sep = ""))
	
	orthologous = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Orthology/Grape-Arabidopsis/All.txt", header = T)
	gene.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Vitis_vinifera/gene.names.txt")
	map = read.delim("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Orthology/Grape-Arabidopsis/Map.12x.plaza.txt")
	Ref = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Reference_SD.txt", header = T)
	
	
	Score = matrix(NA, length(V(g)), 6)
	colnames(Score) = c("Gene", "orthologous", "orthologous.cluster", "Gene.cluster.size",
	 "orthologous.cluster.size", "Arab.orthol")
	Score = as.data.frame(Score)

	i = 1

	memb.comun.grape = data.frame("name" = gene.names[neighborhood(g.grape, 1, i, mode = c("out"))[[1]],])
        memb.comun.grape.para = merge(memb.comun.grape, map, by = "name")
        ## Remove duplicates
        memb.comun.grape.para = memb.comun.grape.para[!(duplicated(memb.comun.grape.para$name)
                                                | duplicated(memb.comun.grape.para$name, fromLast = TRUE)), ]
        memb.comun.grape.para = merge(memb.comun.grape.para, orthologous, by = "PLAZA.ID")
        memb.comun.grape.para.size = length(unique(memb.comun.grape.para$name))




	for (j in 1:length(V(g))){
			print(j)
			
			memb.comun.Arab = data.frame("gene_id.y" = Ref[neighborhood(g, 1, j, mode = c("out"))[[1]],1])
			Score[j, "Gene"] = i
			Score[j, "orthologous"] = j
			Score[j, "orthologous.cluster"] = length(unique(merge(memb.comun.Arab, memb.comun.grape.para, by = "gene_id.y")$name))
			Score[j, "Arab.orthol"] = length(unique(merge(orthologous, memb.comun.Arab, by = "gene_id.y")$PLAZA.ID))
			Score[j, "Gene.cluster.size"] = dim(memb.comun.grape)[1]
			Score[j, "orthologous.cluster.size"] = dim(memb.comun.Arab)[1]
	}


	Score["Normalized.orthologous.cluster"] = Score$orthologous.cluster/(Score$Gene.cluster.size + Score$orthologous.cluster.size)
	Score["Normalized.orthologous.cluster.grape"] = Score$orthologous.cluster/memb.comun.grape.para.size
	Score["Normalized.orthologous.cluster.Arab"] = Score$orthologous.cluster/Score$Arab.orthol
	Score["Normalized.orthologous.whitin.grape"] = memb.comun.grape.para.size/Score$Gene.cluster.size
	Score["Normalized.orthologous.whitin.Arab"] = Score$Arab.orthol/Score$orthologous.cluster.size
	Score[is.na(Score)] = 0

	

	write.table(Score, file = paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/",NAME[k],"/Score99/Scores_V2/Score1.txt", sep = ""), row.names = F, quote = F, sep = "\t")

# }


