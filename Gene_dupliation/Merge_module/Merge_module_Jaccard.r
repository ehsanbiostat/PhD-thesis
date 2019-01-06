
library(igraph)
library(ggplot2)
library(foreach)
library(data.table)

condition.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/condition_name.txt", header = T)
NAME = names(table(condition.names))
k = 1
paralogous = read.table("~/Gene_duplication/Seed_genes/paralogous_NoFilter.txt", header = T)
PARA = sort(unique(paralogous$Ref))

setwd(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Local/Compendium_wise/",NAME[k],"/Score99/Scores_V2/", sep = ""))


# Extract Jaccard index only for anchorpoints modules
Jacc = matrix(NA, length(PARA), length(PARA))
for (i in 1:length(PARA)) {
	print(i)
	Score = read.delim(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Local/Compendium_wise/",NAME[k],"/Score99/Scores_V2/Score",PARA[i],".txt", sep = ""))
	Jacc[,i] = Score[PARA, "Jaccard"]
}

write.table(Jacc, file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Local/Compendium_wise/abiotic/Score99/Scores_V2/Jaccard_All_Anch.txt", row.names = F, col.names=F, sep = "\t")

Jacc = fread("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Local/Compendium_wise/abiotic/Score99/Scores_V2/Jaccard_All_Anch.txt")
Score = fread(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Local/Compendium_wise/",NAME[k],"/Score99/Scores_V2/Score1.txt", sep = ""), header = T)
Jacc.3 = Jacc[Score[PARA,"Paralogous.cluster.size"] > 2, Score[PARA,"Paralogous.cluster.size"] > 2]
PP2 = data.frame("x" = Jacc.3[upper.tri(Jacc.3, diag = F)] * 2)
ggplot(PP2, aes(x = x)) + geom_density(colour="darkgreen", size=0.5, fill = "blue") + scale_x_continuous(limits=c(0, 1))

# Creat a data frame of Jaccard index and all combination of anchorpoints
Jacc.C = data.frame("Jacc" = Jacc.C)
A = data.frame("Gene" = rep(PARA, length(PARA)), "Paralogous" = c(sapply(PARA, rep , length(PARA))))
Jaccard = cbind(A, Jacc.C)
Jaccard = Jaccard[Jaccard$Gene != Jaccard$Paralogous,]
Score = Jaccard[Jaccard$Jacc >= 0.4, ]


g = graph.data.frame(Score[, 1:2])
E(g)$weight = Score$Jacc
Mat = get.adjacency(g, attr = "weight")
Mat = as.matrix(Mat)
Mat[Mat != 0] = 1
Mat = 0.5 - Mat
Mat[Mat != 0.5] = 0
diag(Mat) = 0
MAt.dist = as.dist(Mat)
clust <- hclust(MAt.dist, method = "single")
ord <- order(cutree(clust, k = 3))
coph <- cophenetic(clust)

image(as.matrix(MAt.dist)[ord, ord], main = "Original distances")
image(as.matrix(coph)[ord, ord], main = "Cophenetic distances")
Clusters = data.frame("PARA" = names(cutree(clust, k = 3)), "cluster" = cutree(clust, k = 3))
Clusters[Clusters$cluster == 3, "PARA"]

load(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/",NAME[k],"/Graph_V2.RData", sep = ""))
sort(table(unlist(neighborhood(g, 1, Clusters[Clusters$cluster == 3, "PARA"], mode = c("out")))))

# ---------------------------------------------------------------------------------------------------------------------------------------------------------
# Randomization test
# ---------------------------------------------------------------------------------------------------------------------------------------------------------
Score = read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Local/Compendium_wise/development/Score99/Scores_V2/Score1.txt", header = T)
cluster.size = sort(unique(Score[PARA, "Paralogous.cluster.size"]))
cluster.size = cluster.size[cluster.size > 2]
# Multiple plots function
source("/ngsprojects/iwt/ehsab/Gene_duplication/Scripts/General/Multiple_plot.r")

# Sample from nodes with the cluster's size


Int = function(x,y){length(intersect(x, y))/length(union(x, y))}
plots = list()
Quantile = foreach (i = 1:9, .combine = rbind) %do% { 
	cluster = lapply(cluster.size, function(x){sample(1:nrow(Score), x)})
	Anc = sapply(1:length(cluster), function(x){mapply(Int, cluster, cluster[x])})
	PP = data.frame("x" = Anc[upper.tri(Anc, diag = F)] * 2)
	p1 = ggplot(PP, aes(x = x)) + geom_density(colour="darkgreen", size=0.5, fill = "blue") + scale_x_continuous(limits=c(0, 1))
	plots[[i]] = p1
	quantile(PP$x)
}
multiplot(plotlist = plots, cols = 3)
# ---------------------------------------------------------------------------------------------------------------------------------------------------------



