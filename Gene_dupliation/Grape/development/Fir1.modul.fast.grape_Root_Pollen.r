

	 
library(fields)
library(igraph)
library(foreach)
library(doMC)
# registerDoMC(14)

condition.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/condition_name.txt", header = T)
NAME = names(table(condition.names))

k = 6
	
## Grape modules
load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Graph_V2_90.RData")
		

PollenRootLeafclean = read.delim("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Orthology/Grape-Arabidopsis/PollenRootLeafclean.txt")
gene.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Vitis_vinifera/gene.names.txt")
	
	
Score = matrix(NA, length(V(g)), 4)
colnames(Score) = c("Grape", "Grape.size", "PollenLeaf", "RootLeaf")
Score = as.data.frame(Score)

for (j in 1:length(V(g))){
	print(j)
	memb.comun.grape = data.frame("name" = gene.names[neighborhood(g, 1, j, mode = c("out"))[[1]],])
    sum(memb.comun.grape[, "name"] %in% PollenRootLeafclean[PollenRootLeafclean$Module == "Pollen", "Name"])

    Score[j, "Grape"] = j
	Score[j, "Grape.size"] = nrow(memb.comun.grape)
	Score[j, "PollenLeaf"] = sum(memb.comun.grape[, "name"] %in% PollenRootLeafclean[PollenRootLeafclean$Module == "Pollen", "Name"])
	Score[j, "RootLeaf"] = sum(memb.comun.grape[, "name"] %in% PollenRootLeafclean[PollenRootLeafclean$Module == "Root", "Name"])
}


Score["Normalized.orthologous.cluster"] = Score$orthologous.cluster/(Score$Gene.cluster.size + Score$orthologous.cluster.size)
Score["Normalized.orthologous.cluster.grape"] = Score$orthologous.cluster/memb.comun.grape.para.size
Score["Normalized.orthologous.cluster.Arab"] = Score$orthologous.cluster/Score$Arab.orthol
Score["Normalized.orthologous.whitin.grape"] = memb.comun.grape.para.size/Score$Gene.cluster.size
Score["Normalized.orthologous.whitin.Arab"] = Score$Arab.orthol/Score$orthologous.cluster.size
Score[is.na(Score)] = 0



	write.table(Score, file = paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/",NAME[k],"/Score99/Scores_V2/Score1.txt", sep = ""), row.names = F, quote = F, sep = "\t")

# }


