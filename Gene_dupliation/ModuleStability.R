
library(igraph)
library(foreach)
library(doMC)
library(corrplot)
library(NMF)
library(pheatmap)
options(width = 200)


k = 6
condition.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/condition_name.txt", header = T)
NAME = names(table(condition.names))
NAME.spe = NAME[k]


ref = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Reference_SD.txt", header = T)
paralogous = read.table("~/Gene_duplication/Seed_genes/paralogous_NoFilter.txt", header = T)
PARA = unique(paralogous$Ref)


load(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/",NAME.spe,"/Graph_V2_90.RData", sep = ""))
Anchorpoint = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/2candidates/Genes/Leaf_Pollen_90_Anch.txt")
Anchorpoint1 = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/2candidates/Genes/Root_Leaf_90_Anch.txt")

AnchPollen = as.character(Anchorpoint$Pollen)
AnchRoot = as.character(Anchorpoint1$Root)

AnchPollenLeaf = as.character(Anchorpoint$Leaf)
AnchRootLeaf = as.character(Anchorpoint1$Leaf)

Anch = unique(c(AnchPollen, AnchPollenLeaf, AnchRoot, AnchRootLeaf))

AnchPollenRoot = intersect(AnchRoot, AnchPollen)
AnchLeafLeaf = intersect(AnchRootLeaf, AnchPollenLeaf)

AnchPollenS = AnchPollen[!(AnchPollen %in% AnchPollenRoot)]
AnchPollenLeafS = AnchPollenLeaf[!(AnchPollenLeaf %in% AnchLeafLeaf)]
AnchRootS = AnchRoot[!(AnchRoot %in% AnchPollenRoot)]
AnchRootLeafS = AnchRootLeaf[!(AnchRootLeaf %in% AnchLeafLeaf)]

AnchS = unique(c(AnchPollenS, AnchPollenLeafS, AnchRootS, AnchRootLeafS, AnchPollenRoot, AnchLeafLeaf))

# Membership matrix
annotation = data.frame("Module" = ifelse(Anch %in% AnchPollen, "Pollen", ifelse(Anch %in% AnchPollenLeaf, "PollenLeaf", ifelse(Anch %in% AnchRoot, "Root", "RootLeaf"))))
annotationS = data.frame("Module" = ifelse(AnchS %in% AnchPollenS, "Pollen", ifelse(AnchS %in% AnchPollenLeafS, "PollenLeaf",
	ifelse(AnchS %in% AnchRootS, "Root", ifelse(AnchS %in% AnchRootLeafS, "RootLeaf", ifelse(AnchS %in% AnchPollenRoot, "PollenRoot", "LeafLeaf"))))))



stabilityAnchMatrix = function(Anch) {
	AA = matrix(0, length(Anch), length(Anch))
	for (i in 1:length(Anch)) {
		#memb.comun = data.frame("Ref" = neighborhood(g, 1, Anch[i], mode = c("out"))[[1]])
		#AA[i, which(Anch %in% memb.comun$Ref)] = 1
		memb.comun = V(g)$name[neighborhood(g, 1, Anch[i], mode = "out")[[1]]]
		AA[i, which(Anch %in% memb.comun)] = 1
	}
	colnames(AA) = rownames(AA) = Anch
	return(AA)
}

Matrix = stabilityAnchMatrix(Anch)
MatrixS = stabilityAnchMatrix(AnchS)


ReorderMatrix = function(mydata, row, col) { 
	# set the custom distance and clustering functions, per your example
	hclustfunc <- function(x) hclust(x, method="ward.D2")
	distfunc <- function(x) dist(x, method="binary")

	# perform clustering on rows and columns
	cl.col <- hclustfunc(distfunc(t(mydata)))
	cl.row <- hclustfunc(distfunc(mydata))
	if(row & col) mydata = mydata[cl.row$order, cl.col$order]
	if(!row & col) mydata = mydata[, cl.col$order]
	if(row & !col) mydata = mydata[cl.row$order,]
	return(mydata)
}

RMatrix = ReorderMatrix(Matrix, row = T, col = T)

png("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Local/Validation/Stability/ModuleStabilityMatrix.png", width = 35, height = 20, units = "in", res = 300, pointsize=8)
pdf("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Local/Validation/Stability/ModuleStabilityMatrix.pdf", width = 30, height = 20, pointsize=12)
par(mfrow = c(2,3))
aheatmap(Matrix, color = "black", Rowv = NA, Colv = NA, main = "Module stability",
	annCol = annotation, annRow=annotation, annColors = list(c("blue", "red", "green", "yellow")), cellwidth = 3, cellheight = 3, legend = F)
aheatmap(Matrix, color = "black", Rowv = T, Colv = NA, distfun = "binary", hclustfun = "ward", main = "Module stability with clustering",
	annCol = annotation, annRow=annotation, annColors = list(c("blue", "red", "green", "yellow")), cellwidth = 3, cellheight = 3, legend = F)
aheatmap(Matrix, color = "black", Rowv = T, Colv = T, distfun = "binary", hclustfun = "ward", main = "Module stability with clustering",
	annCol = annotation, annRow=annotation, annColors = list(c("blue", "red", "green", "yellow")), cellwidth = 3, cellheight = 3, legend = F)
aheatmap(MatrixS, color = "black", Rowv = NA, Colv = NA, main = "Module stability No overlap",
	annCol = annotationS, annRow=annotationS, annColors = list(c("blue", "red", "green", "yellow", "black", "brown")), cellwidth = 3, cellheight = 3, legend = F)
aheatmap(MatrixS, color = "black", Rowv = T, Colv = NA, distfun = "binary", hclustfun = "ward", main = "Module stability with clustering No overlap",
	annCol = annotationS, annRow=annotationS, annColors = list(c("blue", "red", "green", "yellow", "black", "brown")), cellwidth = 3, cellheight = 3, legend = F)
aheatmap(MatrixS, color = "black", Rowv = T, Colv = T, distfun = "binary", hclustfun = "ward", main = "Module stability with clustering No overlap",
	annCol = annotationS, annRow=annotationS, annColors = list(c("blue", "red", "green", "yellow", "black", "brown")), cellwidth = 3, cellheight = 3, legend = F)
dev.off()




rownames(annotation) = rownames(Matrix)
rownames(annotationS) = rownames(MatrixS)
pdf("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Local/Validation/Stability/ModuleStabilityMatrix.pdf", width = 30, height = 20, pointsize=12)
par(mfrow = c(2,3))
pheatmap(Matrix, color = colorRampPalette(c("white", "black"))(10), clustering_distance_rows = "binary", clustering_distance_cols = "euclidean", clustering_method = "ward.D2",
	labCol = NULL, annotation_row = annotation, annotation_col = annotation, cellwidth = 5, cellheight = 5, legend = FALSE, fontsize_col = 5, fontsize_row = 5, cluster_col = F,
	border_color = NA, annotation_colors = list(Module = c(Pollen = "red", PollenLeaf = "green", Root = "black", RootLeaf = "brown")), cluster_row = F, main = "Module stability"
)
pheatmap(MatrixS, color = colorRampPalette(c("white", "black"))(10), clustering_distance_rows = "binary", clustering_distance_cols = "euclidean", clustering_method = "ward.D2",
	labCol = NULL, annotation_row = annotationS, annotation_col = annotationS, cellwidth = 5, cellheight = 5, legend = FALSE, fontsize_col = 5, fontsize_row = 5, cluster_col = F,
	border_color = NA, annotation_colors = list(Module = c(LeafLeaf = "darkblue", Pollen = "red", PollenLeaf = "green", PollenRoot = "yellow", Root = "black", RootLeaf = "brown")),
	cluster_row = F, main = "Module stability No overlap"
)
pheatmap(Matrix, color = colorRampPalette(c("white", "black"))(10), clustering_distance_rows = "binary", clustering_distance_cols = "euclidean", clustering_method = "ward.D2",
	labCol = NULL, annotation_row = annotation, annotation_col = annotation, cellwidth = 5, cellheight = 5, legend = FALSE, fontsize_col = 5, fontsize_row = 5, cluster_col = F,
	border_color = NA, annotation_colors = list(Module = c(Pollen = "red", PollenLeaf = "green", Root = "black", RootLeaf = "brown")), main = "Module stability with clustering"
)
pheatmap(MatrixS, color = colorRampPalette(c("white", "black"))(10), clustering_distance_rows = "binary", clustering_distance_cols = "euclidean", clustering_method = "ward.D2",
	labCol = NULL, annotation_row = annotationS, annotation_col = annotationS, cellwidth = 5, cellheight = 5, legend = FALSE, fontsize_col = 5, fontsize_row = 5, cluster_col = F,
	border_color = NA, annotation_colors = list(Module = c(LeafLeaf = "darkblue", Pollen = "red", PollenLeaf = "green", PollenRoot = "yellow", Root = "black", RootLeaf = "brown")),
	main = "Module stability with clustering No overlap"
)
pheatmap(Matrix, color = colorRampPalette(c("white", "black"))(10), clustering_distance_rows = "binary", clustering_distance_cols = "binary", clustering_method = "ward.D2",
	labCol = NULL, annotation_row = annotation, annotation_col = annotation, cellwidth = 5, cellheight = 5, legend = FALSE, fontsize_col = 5, fontsize_row = 5, cluster_col = T,
	border_color = NA, annotation_colors = list(Module = c(Pollen = "red", PollenLeaf = "green", Root = "black", RootLeaf = "brown")), main = "Module stability with clustering"
)
pheatmap(MatrixS, color = colorRampPalette(c("white", "black"))(10), clustering_distance_rows = "binary", clustering_distance_cols = "binary", clustering_method = "ward.D2",
	labCol = NULL, annotation_row = annotationS, annotation_col = annotationS, cellwidth = 5, cellheight = 5, legend = FALSE, fontsize_col = 5, fontsize_row = 5, cluster_col = T,
	border_color = NA, annotation_colors = list(Module = c(LeafLeaf = "darkblue", Pollen = "red", PollenLeaf = "green", PollenRoot = "yellow", Root = "black", RootLeaf = "brown")),
	main = "Module stability with clustering No overlap"
)
dev.off()



write.table(RMatrix, file = "/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Local/Validation/Stability/ModuleStabilityMatrix.txt", quote = F, sep = "\t")








stabilityAnch = function(x) {
		foreach(i = 1:length(x), .combine = c) %do% {
		memb.comun = data.frame("Ref" = V(g)$name[neighborhood(g, 1, x[i], mode = c("out"))[[1]]])
		sum(memb.comun$Ref %in% x)/length(x)
	}
}

stability.AnchorRoot = stabilityAnch(AnchRoot)
stability.AnchorPollen = stabilityAnch(AnchPollen)
stability.AnchorPollenLeaf = stabilityAnch(AnchPollenLeaf)
stability.AnchorRootLeaf = stabilityAnch(AnchRootLeaf)
stability.AnchorS = c(stability.AnchorRootLeaf, stability.AnchorPollenLeaf, stability.AnchorPollen, stability.AnchorRoot)

Anchorpoint1 = data.frame(Anchorpoint1, "stabilityLeaf" = stability.AnchorRootLeaf, "stabilityRoot" = stability.AnchorRoot)
Anchorpoint = data.frame(Anchorpoint, "stabilityPollen" = stability.AnchorPollen, "stabilityLeaf" = stability.AnchorPollenLeaf)

Anchorpoint.All = rbind.data.frame(data.frame("x" = Anchorpoint[, c("Pollen")], "y" = Anchorpoint[, c("stabilityPollen")]),
 data.frame("x" = Anchorpoint[, c("Leaf")], "y" = Anchorpoint[, c("stabilityLeaf")]),
 data.frame("x" = Anchorpoint1[, c("Root")], "y" = Anchorpoint1[, c("stabilityRoot")]),
 data.frame("x" = Anchorpoint1[, c("Leaf")], "y" = Anchorpoint1[, c("stabilityLeaf")]))

mean(Anchorpoint.All[!(duplicated(Anchorpoint.All$x)), 2])
median(Anchorpoint.All[!(duplicated(Anchorpoint.All$x)), 2])

stabilityPermutation = matrix(NA, 1000, length(stability.AnchorPollen))
for (j in 1:1000) {
	print(j)
	Gene = as.character(sample(PARA, length(stability.AnchorPollen)))
	stabilityPermutation[j,] = stabilityAnch(Gene)
}


AA = rbind.data.frame(data.frame("x" = c(stabilityPermutation), "Group" = "Randomized"),
		data.frame("x" = stability.AnchorPollen, "Group" = "Pollen"),
		data.frame("x" = stability.AnchorPollenLeaf, "Group" = "Leaf_P"),
		data.frame("x" = stability.AnchorRoot, "Group" = "Root"),
		data.frame("x" = stability.AnchorRootLeaf, "Group" = "Leaf_R"))




jpeg("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Local/Validation/Stability/ModuleStability.jpeg", width = 19, height = 10, units = "in", res = 800, pointsize=8)
ggplot(AA, aes(x, fill = Group)) + geom_density(alpha = 0.4) + scale_fill_manual(values = c("green", "red", "yellow", "royalblue1", "darkorchid1")) +
 geom_vline(xintercept = c(0.5471698), color = "red", linetype = "longdash") +
 geom_vline(xintercept = c(0.5849057), color = "yellow", linetype = "longdash") +
 geom_vline(xintercept = c(0.1296296), color = "black", linetype = "longdash") +
 geom_vline(xintercept = c(0.5740741), color = "royalblue1", linetype = "longdash") +
 geom_vline(xintercept = c(0.6018519), color = "darkorchid1", linetype = "longdash") +
 annotate("text", x = 0.54, y = 9, label = "Pollen", size = 3.3, angle = 90) +
 annotate("text", x = 0.59, y = 9, label = "Leaf_P", size = 3.3, angle = 90) +
 annotate("text", x = 0.122, y = 3, label = "Randomized", size = 3.3, angle = 90) +
 annotate("text", x = 0.567, y = 9, label = "Root", size = 3.3, angle = 90) +
 annotate("text", x = 0.608, y = 9, label = "Leaf_R", size = 3.3, angle = 90) +
 xlab("Stability")
dev.off()



stability.Anchor = foreach(i = 1:length(PARA)) %dopar% {
	print(i)
	memb.comun = data.frame("Ref" = V(g)$name[neighborhood(g, 1, PARA[i], mode = c("out"))[[1]]])
	memb.comun.para = merge(memb.comun, paralogous, by = "Ref")
	stabilityAnch(as.character(memb.comun.para$Ref))
}



BB = data.frame("x" = sapply(stability.Anchor, median))
BB["Group"] = "All anchorpoint modules"
AAA = rbind(AA, BB)
jpeg("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Local/Validation/Stability/ModuleStability.jpeg", width = 19, height = 10, units = "in", res = 800, pointsize=8)
ggplot(AAA, aes(x, fill = Group)) + geom_density(alpha = 0.4) + scale_fill_manual(values = c("black", "red", "yellow", "green", "royalblue1", "darkorchid1")) +
 geom_vline(xintercept = c(0.4111748), color = "black", linetype = "longdash") +
 geom_vline(xintercept = c(0.5471698), color = "red", linetype = "longdash") +
 geom_vline(xintercept = c(0.5849057), color = "yellow", linetype = "longdash") +
 geom_vline(xintercept = c(0.1296296), color = "black", linetype = "longdash") +
 geom_vline(xintercept = c(0.5740741), color = "royalblue1", linetype = "longdash") +
 geom_vline(xintercept = c(0.6018519), color = "darkorchid1", linetype = "longdash") +
 annotate("text", x = 0.40, y = 9, label = "All anchorpoint modules", size = 3.3, angle = 90) +
 annotate("text", x = 0.54, y = 9, label = "Pollen", size = 3.3, angle = 90) +
 annotate("text", x = 0.59, y = 9, label = "PollenLeaf", size = 3.3, angle = 90) +
 annotate("text", x = 0.122, y = 3, label = "Randomized", size = 3.3, angle = 90) +
 annotate("text", x = 0.567, y = 9, label = "Root", size = 3.3, angle = 90) +
 annotate("text", x = 0.608, y = 9, label = "RootLeaf", size = 3.3, angle = 90) +
 xlab("Stability")
dev.off()








##################################################################################################################################################################################
# SSD stability
##################################################################################################################################################################################

library(igraph)
library(foreach)
library(doMC)
library(corrplot)
library(NMF)
library(pheatmap)
options(width = 200)


k = 6
condition.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/condition_name.txt", header = T)
NAME = names(table(condition.names))
NAME.spe = NAME[k]


ref = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Reference_SD.txt", header = T)
paralogous = read.table("~/Gene_duplication/Seed_genes/paralogous_NoFilter.txt", header = T)
PARA = unique(paralogous$Ref)


load(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/",NAME.spe,"/Graph_V2_90.RData", sep = ""))
Anchorpoint = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/2candidates/Genes/Leaf_Pollen_90_Anch.txt")
Anchorpoint1 = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/2candidates/Genes/Root_Leaf_90_Anch.txt")

AnchPollen = as.character(Anchorpoint$Pollen)
AnchRoot = as.character(Anchorpoint1$Root)

AnchPollenLeaf = as.character(Anchorpoint$Leaf)
AnchRootLeaf = as.character(Anchorpoint1$Leaf)


SSDPollenLeaf = as.character(read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/SSDPollenLeafGenes.txt")[,2])
SSDPollen = as.character(read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/SSDPollenGenes.txt")[,2])
SSDRootLeaf = as.character(read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/SSDRootLeafGenes.txt")[,2])
SSDRoot = as.character(read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/SSDRootGenes.txt")[,2])


Anch = unique(c(AnchPollen, AnchPollenLeaf, AnchRoot, AnchRootLeaf))


stabilityAnchMatrix = function(Anch, SSD) {
	AA = matrix(0, length(Anch), length(SSD))
	for (i in 1:length(Anch)) {
		#memb.comun = data.frame("Ref" = neighborhood(g, 1, Anch[i], mode = c("out"))[[1]])
		#AA[i, which(Anch %in% memb.comun$Ref)] = 1
		memb.comun = V(g)$name[neighborhood(g, 1, Anch[i], mode = "out")[[1]]]
		AA[i, which(SSD %in% memb.comun)] = 1
	}
	colnames(AA) = SSD
	rownames(AA) = Anch
	return(AA)
}


PollenStability = stabilityAnchMatrix(AnchPollen, SSDPollen)
PollenLeafStability = stabilityAnchMatrix(AnchPollenLeaf, SSDPollenLeaf)
RootStability = stabilityAnchMatrix(AnchRoot, SSDRoot)
RootLeafStability = stabilityAnchMatrix(AnchRootLeaf, SSDRootLeaf)



pdf("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Local/Validation/Stability/ModuleStabilityMatrixSSDPollen.pdf", width = 15, height = 10, pointsize=12)
pheatmap(PollenStability, color = colorRampPalette(c("white", "black"))(10), clustering_distance_rows = "binary", clustering_distance_cols = "binary", clustering_method = "ward.D2",
	labCol = NULL,  cellwidth = 3, cellheight = 3, legend = FALSE, fontsize_col = 5, fontsize_row = 5, cluster_col = T,
	border_color = NA,  cluster_row = T, main = "Module stability"
)
dev.off()

pdf("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Local/Validation/Stability/ModuleStabilityMatrixSSDPollenLeaf.pdf", width = 15, height = 10, pointsize=12)
pheatmap(PollenLeafStability, color = colorRampPalette(c("white", "black"))(10), clustering_distance_rows = "binary", clustering_distance_cols = "binary", clustering_method = "ward.D2",
	labCol = NULL,  cellwidth = 3, cellheight = 3, legend = FALSE, fontsize_col = 5, fontsize_row = 5, cluster_col = T,
	border_color = NA,  cluster_row = T, main = "Module stability"
)
dev.off()

pdf("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Local/Validation/Stability/ModuleStabilityMatrixSSDRoot.pdf", width = 15, height = 10, pointsize=12)
pheatmap(RootStability, color = colorRampPalette(c("white", "black"))(10), clustering_distance_rows = "binary", clustering_distance_cols = "binary", clustering_method = "ward.D2",
	labCol = NULL,  cellwidth = 3, cellheight = 3, legend = FALSE, fontsize_col = 5, fontsize_row = 5, cluster_col = T,
	border_color = NA,  cluster_row = T, main = "Module stability"
)
dev.off()

pdf("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Local/Validation/Stability/ModuleStabilityMatrixSSDRootLeaf.pdf", width = 15, height = 10, pointsize=12)
pheatmap(RootLeafStability, color = colorRampPalette(c("white", "black"))(10), clustering_distance_rows = "binary", clustering_distance_cols = "binary", clustering_method = "ward.D2",
	labCol = NULL,  cellwidth = 3, cellheight = 3, legend = FALSE, fontsize_col = 5, fontsize_row = 5, cluster_col = T,
	border_color = NA,  cluster_row = T, main = "Module stability"
)
dev.off()


annotation_row = annotation, annotation_col = annotation, annotation_colors = list(Module = c(Pollen = "red", PollenLeaf = "green", Root = "black", RootLeaf = "brown")),

##################################################################################################################################################################################









##################################################################################################################################################################################
# SSD & LSD
##################################################################################################################################################################################

DupRootLeaf = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/DupRootLeaf.txt")
DupPollenLeaf = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/DupPollenLeaf.txt")
DupRootLeafGene = unique(c(as.character(DupRootLeaf$Leaf), as.character(DupRootLeaf$Root)))
DupRootGene = unique(as.character(DupRootLeaf$Root))
DupLeafRootGene = unique(as.character(DupRootLeaf$Leaf))
DupPollenLeafGene = unique(c(as.character(DupPollenLeaf$Leaf), as.character(DupPollenLeaf$Root)))
DupPollenGene = unique(as.character(DupPollenLeaf$Pollen))
DupLeafPollenGene = unique(as.character(DupPollenLeaf$Leaf))

AnchSSD = unique(c(Anch, DupRootLeafGene, DupPollenLeafGene))
Matrix = stabilityAnchMatrix(AnchSSD)
RMatrix = ReorderMatrix(Matrix, row = T, col = T)

annotation = data.frame("Module" = ifelse(AnchSSD %in% c(AnchPollen, DupPollenGene),
	"Pollen", ifelse(AnchSSD %in% c(AnchPollenLeaf, DupLeafPollenGene), "PollenLeaf", ifelse(AnchSSD %in% c(AnchRoot, DupRootGene), "Root", "RootLeaf"))),
	"Duplication" = ifelse(AnchSSD %in% Anch, "LSD", "SSD"))

rownames(annotation) = rownames(Matrix)
pdf("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Local/Validation/Stability/ModuleStabilityMatrix_SSD_LSD.pdf", width = 30, height = 20, pointsize=12)
pheatmap(Matrix, color = colorRampPalette(c("white", "black"))(10), clustering_distance_rows = "binary", clustering_distance_cols = "euclidean", clustering_method = "ward.D2",
	labCol = NULL, annotation_row = annotation, annotation_col = annotation, cellwidth = 1, cellheight = 1, legend = FALSE, fontsize_col = 1, fontsize_row = 1, cluster_col = F,
	border_color = NA, annotation_colors = list(Module = c(Pollen = "red", PollenLeaf = "green", Root = "black", RootLeaf = "brown")), cluster_row = F, main = "Module stability")
pheatmap(Matrix, color = colorRampPalette(c("white", "black"))(10), clustering_distance_rows = "binary", clustering_distance_cols = "euclidean", clustering_method = "ward.D2",
	labCol = NULL, annotation_row = annotation, annotation_col = annotation, cellwidth = 1, cellheight = 1, legend = FALSE, fontsize_col = 1, fontsize_row = 1, cluster_col = F,
	border_color = NA, annotation_colors = list(Module = c(Pollen = "red", PollenLeaf = "green", Root = "black", RootLeaf = "brown")), cluster_row = T, main = "Module stability"
)
pheatmap(Matrix, color = colorRampPalette(c("white", "black"))(10), clustering_distance_rows = "binary", clustering_distance_cols = "euclidean", clustering_method = "ward.D2",
	labCol = NULL, annotation_row = annotation, annotation_col = annotation, cellwidth = 1, cellheight = 1, legend = FALSE, fontsize_col = 1, fontsize_row = 1, cluster_col = T,
	border_color = NA, annotation_colors = list(Module = c(Pollen = "red", PollenLeaf = "green", Root = "black", RootLeaf = "brown")), cluster_row = T, main = "Module stability"
)
dev.off()

write.table(RMatrix, file = "/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Local/Validation/Stability/ModuleStabilityMatrix_SSD_LSD_Order.txt", quote = F, sep = "\t")
write.table(Matrix, file = "/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Local/Validation/Stability/ModuleStabilityMatrix_SSD_LSD.txt", quote = F, sep = "\t")

##################################################################################################################################################################################











##################################################################################################################################################################################
# Module stability
##################################################################################################################################################################################
ModulePollen = read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/Pollen_90.txt")
ref = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Reference_SD.txt", header = T)
Gene = ref[ref[,1] %in% ModulePollen$V1,2]

stability.all = foreach(i = 1:length(Gene), .combine = c) %do% {
	print(i)
	memb.comun = data.frame("Ref" = neighborhood(g, 1, Gene[i], mode = c("out"))[[1]])
	sum(memb.comun$Ref %in% Gene)/length(Gene)
}

load("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Local/Validation/StabilityAnchorRandomization.RData")

top20[top20$Threshold == 90,"GeneStabilityScore"] = data.frame(matrix(stability.Anchor, length(Gene)/2, 2))[,1]
top20[top20$Threshold == 90,"ParaStabilityScore"] = data.frame(matrix(stability.Anchor, length(Gene)/2, 2))[,2]


# Image scores integration
top20[top20$Threshold == 90,"GeneImageScore"] = Image_scoring[,3]
top20[top20$Threshold == 90,"ParaImageScore"] = Image_scoring[,4]


# Stability for all anchorpoint modules
stability.Anchor = foreach(i = 1:length(PARA), .combine = rbind) %dopar% {
	print(i)
	memb.comun = data.frame("Ref" = neighborhood(g, 1, PARA[i], mode = c("out"))[[1]])
	memb.comun.para = merge(memb.comun, paralogous, by = "Ref")
	Anchors = memb.comun.para$Ref
	A = foreach(k = 1:length(Anchors), .combine = c) %do% {
		sum(neighborhood(g, 1, Anchors[k], mode = c("out"))[[1]] == PARA[i])
		# sum(sum(neighborhood(g, 1, Anchors[k], mode = c("out"))[[1]] %in% Gene[c(2:4)]) == 3)
	}
	data.frame("Seed" = PARA[i], sum(A), length(A))
}


stability.Anchor = foreach(i = 1:length(PARA)) %dopar% {
	print(i)
	memb.comun = data.frame("Ref" = neighborhood(g, 1, PARA[i], mode = c("out"))[[1]])
	memb.comun.para = merge(memb.comun, paralogous, by = "Ref")
	Anchors = memb.comun.para$Ref
	stabilityAnch(Anchors)
}



# Permutation test for stability test

stabilityPermutation90 = foreach(m = 1:length(Gene), .combine = rbind) %do% {
	memb.comun = data.frame("Ref" = neighborhood(g, 1, Gene[m], mode = c("out"))[[1]])
	# memb.comun = data.frame("Ref" = as.numeric(unlist(gene_score_module13393[3])))
	memb.comun.para = merge(memb.comun, paralogous, by = "Ref")
	Anchors = memb.comun.para$Ref
	B = foreach(i = 1:1000, .combine = c) %dopar% {
		print(i)
		A = foreach(k = 1:length(Anchors), .combine = c) %do% {
			sum(sample(1:19285, length(neighborhood(g, 1, Anchors[k], mode = c("out"))[[1]])) %in% Gene[m])
			# sum(sum(sample(1:19285, length(neighborhood(g, 1, Anchors[k], mode = c("out"))[[1]])) %in% Gene[c(2:4)]) == 3)
		}
		sum(A)/length(A)
	}
	B
}




# Jaccard matrix
Jaccard.Matrix = foreach (i = 1:length(Gene), .combine = rbind) %do% {
	JJ = foreach(j = 1:length(Gene), .combine = c) %do% {
		length(intersect(neighborhood(g, 1, Gene[i], mode = c("out"))[[1]], neighborhood(g, 1, Gene[j], mode = c("out"))[[1]]))/length(union(neighborhood(g, 1, Gene[i], mode = c("out"))[[1]], neighborhood(g, 1, Gene[j], mode = c("out"))[[1]]))
	}
	JJ
}


# Venn diagram
GeneSet1 = Gene[c(1,6,5)]
memb.comun = data.frame("Ref" = unlist(neighborhood(g,1,GeneSet1, mode = c("out"))))
memb.comun.para = merge(memb.comun, paralogous, by = "Ref")
universe = memb.comun.para$Ref

Counts <- matrix(NA, nrow=length(universe), ncol=length(GeneSet1))
dim(Counts)
for (i in 1:length(universe)) {
	print(i)
	for (j in 1:length(GeneSet1)) {
		Counts[i,j] <- universe[i] %in% merge(data.frame("Ref" = neighborhood(g,1,GeneSet1, mode = c("out"))[[j]]), paralogous, by = "Ref")$Ref
	}
}

GeneSet1.name = Gene.name[c(1,6,5)]
colnames(Counts) <- GeneSet1.name
library(limma)
cols<-c("Red", "Green", "Blue")
vennDiagram(vennCounts(Counts), circle.col=cols)


Geneset1 = c(10182,1352,2663,6258,1932,13393,7682,13899,8399,15525,13906,6874,9603,2228,39,1694,10546,6713)
Geneset2 = c(5792,15713,399,16681,17367,17229,2858,16546,7903,14025,440,18332,12364,10298,12450)

# Overlapping in genes, looking at gene score and figure out genes with high gene scores are overlapped
# across candidate seed genes


GeneScore1 = foreach (i = Geneset1, .combine = cbind) %do% {
	load(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/development/Raw/G",i,".RData", sep = ""))
	# gene_scores = cbind(gene_scores, 1:19285)
	# gene_scores = gene_scores[order(-gene_scores[,1]),]
	gene_scores
}


Jaccard.Matrix.GeneScore = foreach (i = 1:length(Geneset1), .combine = rbind) %do% {
	JJ = foreach(j = 1:length(Geneset1), .combine = c) %do% {
		length(intersect(GeneScore[1:300, i], GeneScore[1:300, j]))/length(union(GeneScore[1:300, i], GeneScore[1:300, j]))
	}
	JJ
}


# Correlationplot for gene scores - two set of anchorpoints which show a clear gene expression (anti)correlation
# Comparing the nice pattern observed in development compendium across other GeneSetCorrplotOtherComepndia
jpeg(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/GeneSetCorrplotOtherComepndia.jpeg", sep = ""), width = 11, height = 7, units = "in", res = 500, pointsize=8)
par(mfrow = c(3,5))
for(k in 1:length(NAME)) {
	GeneScore = foreach (i = Geneset, .combine = cbind) %do% {
		load(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/",NAME[k],"/Raw/G",i,".RData", sep = ""))
		gene_scores
	}
	corrplot(cor(GeneScore), method = "color", order = "AOE", title = NAME[k], tl.cex=1, tl.pos = "n", cl.pos = "n")
}
dev.off()




############################################################################################################################################################
#Old staffs
############################################################################################################################################################
# Read module scores by different thresholds
a99 = read.delim(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Local/Compendium_wise/",NAME[k],"/Score99/Scores_V2/Anch_Score_5000.txt", sep = ""))
a95 = read.delim(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Local/Compendium_wise/",NAME[k],"/Score95/Scores_V2/Anch_Score_5000.txt", sep = ""))
a90 = read.delim(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Local/Compendium_wise/",NAME[k],"/Score90/Scores_V2/Anch_Score_5000.txt", sep = ""))
a = rbind(a99, a95, a90)
a["Threshold"] = c(sapply(list(99,95,90), function(x){rep(x, nrow(a99))}))

# Significant number of shared anchorpoints, permutation test
load("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Local/Compendium_wise/biotic/Randomization/95/Threshold.RData")

aMerge = merge(a, Threshold, by.x = c("Gene.cluster.size", "Paralogous.cluster.size"), by.y = c("cluster.one", "cluster.two"))
aMergeSig = aMerge[aMerge$Paralogous.cluster > aMerge$Threshold95 & aMerge$Cluster.overlap < aMerge$Threshold_Jac95, ]
aMergeSig = aMergeSig[order(aMergeSig$Gene, aMergeSig$Paralogous, aMergeSig$Threshold), ]

# Screen for Module's pair by #membership overlap

QQ = foreach(k = 1:20, .combine = rbind) %do% {

	ZeroJac = aMergeSig[aMergeSig$Cluster.overlap == k , ]

	ZZ = foreach(i = 1:nrow(ZeroJac), .combine = rbind) %do% {
		AA = aMergeSig[aMergeSig$Gene == ZeroJac$Gene[i] & aMergeSig$Paralogous == ZeroJac$Paralogous[i], ]
		if(nrow(AA) > 1) data.frame("Gene" = ZeroJac$Gene[i], "Paralogous" = ZeroJac$Paralogous[i], "Score" = abs(diff(AA[, "Cluster.overlap"])), 
			"Pvalue" = prop.test(AA[, 6], AA[, 1] + AA[, 2] - AA[, 6])$p.value, "Cluster.overlap" = k)
	}
	ZZ
}
QQ = QQ[!duplicated(QQ[, 1:2]), ]
QQ = QQ[order(-QQ$Pvalue, -QQ$Cluster.overlap), ]
QQSig = QQ[QQ$Pvalue > 0.05,]
top20 = QQSig[!duplicated(QQSig[, 3:4]), ]

# top10 = ZeroJac[duplicated(ZeroJac[, c("Gene", "Paralogous")]), c("Gene", "Paralogous")]
top20 = merge(aMergeSig, top20, by = c("Gene", "Paralogous"))[, c("Gene", "Paralogous", "Gene.cluster.size", "Paralogous.cluster.size", "Paralogous.cluster", 
 	"Cluster.overlap.x", "Cluster.overlap.Anch", "Threshold")]


write.table(top20, file = paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/",NAME.spe,"/3threshold_20Jaccard.txt", sep = ""), sep = "\t", quote = F, row.names = F)


top = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/3threshold_0Jaccard.txt")
Gene = unique(c(top$Gene, top$Paralogous))

rownames(ref) = ref[,1]
Gene.name = ref[Gene, 1]
paralogous = read.table("~/Gene_duplication/Seed_genes/paralogous_NoFilter.txt", header = T)
PARA = unique(paralogous$Ref)

load(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/",NAME.spe,"/Graph_V2.RData", sep = ""))


stability.All = foreach(i = 1:length(Gene), .combine = c) %dopar% {
	print(i)
	memb.comun = data.frame("Ref" = neighborhood(g, 1, Gene[i], mode = c("out"))[[1]])
	memb.comun.para = merge(memb.comun, paralogous, by = "Ref")
	Anchors = memb.comun.para$Ref
	A = foreach(k = 1:length(Anchors), .combine = rbind) %do% {
			B = foreach(m = 1:length(Anchors), .combine = c) %do% {
				sum(neighborhood(g, 1, Anchors[k], mode = c("out"))[[1]] == Anchors[m])
			}
		B
	}
	sum(A)/length(A)
}


stability.Anchor = foreach(i = 1:length(Gene), .combine = c) %do% {
	print(i)
	memb.comun = data.frame("Ref" = neighborhood(g, 1, Gene[i], mode = c("out"))[[1]])
	# memb.comun = data.frame("Ref" = as.numeric(unlist(gene_score_module13393[3])))
	memb.comun.para = merge(memb.comun, paralogous, by = "Ref")
	Anchors = memb.comun.para$Ref
	A = foreach(k = 1:length(Anchors), .combine = c) %do% {
		sum(neighborhood(g, 1, Anchors[k], mode = c("out"))[[1]] == Gene[i])
		# sum(sum(neighborhood(g, 1, Anchors[k], mode = c("out"))[[1]] %in% Gene[c(2:4)]) == 3)
	}
	sum(A)/length(A)
}

