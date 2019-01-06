


library(igraph)
options(width = 200)
library(vegan)
library(corrplot)
library(MASS)
library(ggplot2)


source("/ngsprojects/iwt/ehsab/Gene_duplication/Scripts/Functions/Scaling.R")



###############################################################################################################################################################
## Grape
###############################################################################################################################################################
data.grape = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Vitis_vinifera/GE.mean.SD.txt", nrows = 23629, comment.char = "", header = T)
cond.devel = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Vitis_vinifera/development.annotaion.txt")
colnames(data.grape) = cond.devel[,2]
data.grape.scale = Scaling(data.grape)
gene.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Vitis_vinifera/gene.names.txt")
gene.names = as.matrix(gene.names)
## Grape modules
load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Graph_V2.R")
g.grape = g

###############################################################################################################################################################




###############################################################################################################################################################
## Arabidopsis
###############################################################################################################################################################
data.arab = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Integrated_SD.txt", colClasses = "numeric", nrows = 19285, comment.char = "", header = T)
condition.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/condition_name.txt", header = T)
NAME = names(table(condition.names))
data.arab = as.matrix(data.arab[, condition.names == "development"])
Ref = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Reference_SD.txt", header = T)
rownames(data.arab) = Ref[,1]
k = 6
## Arabidopsis modules
load(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/",NAME[k],"/Graph_V2.RData", sep = ""))
###############################################################################################################################################################



orthologous = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Orthology/Grape-Arabidopsis/All.txt", header = T)
onetone = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Orthology/Grape-Arabidopsis/ZhenProcessed.txt", header = T)
gene.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Vitis_vinifera/gene.names.txt")
map = read.delim("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Orthology/Grape-Arabidopsis/Map.12x.plaza.txt")
Ref = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Reference_SD.txt", header = T)
# paralogous = read.table("~/Gene_duplication/Seed_genes/paralogous_total.txt", header = T)
paralogous = read.table("~/Gene_duplication/Seed_genes/paralogous_NoFilter.txt", header = T)
rownames(data.arab) = Ref[,1]


##########################################################################
# Extracting APs with corresponding ORTH
##########################################################################
OrthZhen = read.delim("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Orthology/Grape-Arabidopsis/PollenRootLeafclean.txt")
Pollen = unique(as.character(OrthZhen[OrthZhen$Module == "Pollen", "name"]))
Root = unique(as.character(OrthZhen[OrthZhen$Module == "Root", "name"]))
RootPollen = intersect(Pollen, Root)
Pollen = Pollen[!(Pollen %in% RootPollen)]
Root = Root[!(Root %in% RootPollen)]
##########################################################################


###############################################################################################################################################################
# MDS plot for the all conditions in the grapevone development compendium 
###############################################################################################################################################################


CorPlot = function(Data, Genes, order, address){
	require("corrplot")
	CC = cor(Data[Genes,])
	pdf(address, width = 10, height = 10)
	corrplot(CC, tl.cex = 0.4, order = order, method = "color", , hclust.method="ward")
	dev.off()
}

##########################################################################

Dist = dist(t(data.grape.scale))

Dist.Cor = as.dist(1 - CC)
cmd = as.data.frame(cmdscale(Dist.Cor))
cmd$type = rownames(cmd)
kclus <- kmeans(t(data.grape.scale), centers= 10, iter.max=1000, nstart=10000)
kclus <- kmeans(cmd[,1:2], centers= 6, iter.max=1000, nstart=10000)
cmd$group = kclus$cluster




CorPlot(data.grape.scale, Pollen, order = "FPC", address = "/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Grape/Plots/ConditionCorrelationMatrixFPCPollen.pdf")
CorPlot(data.grape.scale, Root, order = "FPC", address = "/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Grape/Plots/ConditionCorrelationMatrixFPCRoot.pdf")

CorPlot(data.grape.scale, Pollen, order = "hclust", address = "/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Grape/Plots/ConditionCorrelationMatrixHclustPollen.pdf")
CorPlot(data.grape.scale, Root, order = "hclust", address = "/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Grape/Plots/ConditionCorrelationMatrixHclustRoot.pdf")

CorPlot(data.grape.scale, gene.names, order = "FPC", address = "/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Grape/Plots/ConditionCorrelationMatrixFPC.pdf")
CorPlot(data.grape.scale, gene.names, order = "hclust", address = "/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Grape/Plots/ConditionCorrelationMatrixHclust.pdf")


pdf("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Grape/Plots/ConditionCorrelationMDSKmeans.pdf", width = 15, height = 10)
ggplot(cmd, aes(x = V1, y = V2, color = factor(group))) + geom_point() + geom_text(aes(label = type), cex = 3) + scale_colour_discrete(guide = FALSE)
dev.off()

# pdf("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Grape/Plots/ConditionCorrelationMDS.pdf", width = 15, height = 10)
# autoplot(sammon(Dist.Cor), shape = FALSE, label.colour = 'blue', label.size = 2.5)
# dev.off()

###############################################################################################################################################################
