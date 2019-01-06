
library(igraph)
library(foreach)
library(doMC)
library(corrplot)
library(NMF)
library(pheatmap)
options(width = 200)


##########################################################################
# Functions
##########################################################################


stabilityAnchMatrix = function(Anch) {
	AA = matrix(0, length(Anch), length(Anch))
	for (i in 1:length(Anch)) {
		memb.comun = data.frame("Ref" = neighborhood(g, 1, Anch[i], mode = c("out"))[[1]])
		AA[i, which(Anch %in% memb.comun$Ref)] = 1
	}
	colnames(AA) = rownames(AA) = as.character(gene.names[Anch, 1])
	return(AA)
}



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



##########################################################################
# Removing ORTH with multiple gene names to one PLAZA ID
##########################################################################
OrthZhen = read.delim("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Orthology/Grape-Arabidopsis/PollenRootLeafclean.txt")
PLAZA.Dup = as.character(OrthZhen$PLAZA.ID[duplicated(OrthZhen$PLAZA.ID)])
PLAZA.Dup = PLAZA.Dup[sapply(1:length(PLAZA.Dup), function(i){length(unique(OrthZhen[OrthZhen$PLAZA.ID %in% PLAZA.Dup[i], "name"]))}) > 1]
OrthZhen = OrthZhen[!(OrthZhen$PLAZA.ID %in% PLAZA.Dup),]
##########################################################################


##########################################################################
# Extracting APs with corresponding ORTH
##########################################################################
Pollen = unique(as.character(OrthZhen[OrthZhen$Module == "Pollen", "name"]))
Root = unique(as.character(OrthZhen[OrthZhen$Module == "Root", "name"]))
RootPollen = intersect(Pollen, Root)
Pollen = Pollen[!(Pollen %in% RootPollen)]
Root = Root[!(Root %in% RootPollen)]
##########################################################################



##########################################################################
# graph data and gene names for grapevine
##########################################################################
load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Graph_V2_90.RData")
gene.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Vitis_vinifera/gene.names.txt")
rownames(gene.names) = gene.names[,1]
gene.names$Row = 1:nrow(gene.names)
##########################################################################

##########################################################################
# Extract genes for the graph
##########################################################################
PollenGrape = gene.names[Pollen,"Row"]
RootGrape = gene.names[Root,"Row"]
Anch = (unique(c(Pollen, Root)))
PollenRootGrape = gene.names[Anch,"Row"]
##########################################################################


MatrixPollen = stabilityAnchMatrix(PollenGrape)
MatrixRoot = stabilityAnchMatrix(RootGrape)
MatrixPollenRoot = stabilityAnchMatrix(PollenRootGrape)


RMatrix = ReorderMatrix(Matrix, row = T, col = T)



annotation = data.frame("Module" = c(rep("Pollen", length(Pollen)), rep("Root", length(Root))))
rownames(annotation) = rownames(MatrixPollenRoot)

pheatmap(MatrixPollenRoot, color = colorRampPalette(c("white", "black"))(10), clustering_distance_rows = "binary", clustering_distance_cols = "binary", clustering_method = "ward.D2",
	labCol = NULL, annotation_row = annotation, annotation_col = annotation, cellwidth = 10, cellheight = 10, legend = FALSE, fontsize_col = 8, fontsize_row = 8, cluster_col = F,
	border_color = NA, annotation_colors = list(Module = c(Pollen = "red", Root = "darkblue")), main = "Module stability for Root & Pollen orthologues", cluster_row = F,
	gaps_row = length(Pollen), gaps_col = length(Pollen),
	filename = "/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Local/Validation/Stability/ModuleStabilityMatrixGrape.pdf"
)




pdf("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Local/Validation/Stability/ModuleStabilityMatrixGrape.pdf", width = 30, height = 20, pointsize=12)
par(mfrow = c(3,3))
pheatmap(MatrixPollen, color = colorRampPalette(c("white", "black"))(10), clustering_distance_rows = "binary", clustering_distance_cols = "binary", clustering_method = "ward.D2",
	labCol = NULL, cellwidth = 10, cellheight = 10, legend = FALSE, fontsize_col = 8, fontsize_row = 8, cluster_col = F,
	border_color = NA,  cluster_row = F, main = "Module stability for Pollen orthologues"
)
pheatmap(MatrixPollen, color = colorRampPalette(c("white", "black"))(10), clustering_distance_rows = "binary", clustering_distance_cols = "binary", clustering_method = "ward.D2",
	labCol = NULL, cellwidth = 10, cellheight = 10, legend = FALSE, fontsize_col = 8, fontsize_row = 8, cluster_col = F,
	border_color = NA,  cluster_row = T, main = "Module stability for Pollen orthologues"
)
pheatmap(MatrixPollen, color = colorRampPalette(c("white", "black"))(10), clustering_distance_rows = "binary", clustering_distance_cols = "binary", clustering_method = "ward.D2",
	labCol = NULL, cellwidth = 10, cellheight = 10, legend = FALSE, fontsize_col = 8, fontsize_row = 8, cluster_col = T,
	border_color = NA,  cluster_row = T, main = "Module stability for Pollen orthologues"
)
pheatmap(MatrixRoot, color = colorRampPalette(c("white", "black"))(10), clustering_distance_rows = "binary", clustering_distance_cols = "binary", clustering_method = "ward.D2",
	labCol = NULL, cellwidth = 10, cellheight = 10, legend = FALSE, fontsize_col = 8, fontsize_row = 8, cluster_col = F,
	border_color = NA,  cluster_row = F, main = "Module stability for Root orthologues"
)
pheatmap(MatrixRoot, color = colorRampPalette(c("white", "black"))(10), clustering_distance_rows = "binary", clustering_distance_cols = "binary", clustering_method = "ward.D2",
	labCol = NULL, cellwidth = 10, cellheight = 10, legend = FALSE, fontsize_col = 8, fontsize_row = 8, cluster_col = F,
	border_color = NA,  cluster_row = T, main = "Module stability for Root orthologues"
)
pheatmap(MatrixRoot, color = colorRampPalette(c("white", "black"))(10), clustering_distance_rows = "binary", clustering_distance_cols = "binary", clustering_method = "ward.D2",
	labCol = NULL, cellwidth = 10, cellheight = 10, legend = FALSE, fontsize_col = 8, fontsize_row = 8, cluster_col = T,
	border_color = NA,  cluster_row = T, main = "Module stability for Root orthologues"
)
pheatmap(MatrixPollenRoot, color = colorRampPalette(c("white", "black"))(10), clustering_distance_rows = "binary", clustering_distance_cols = "binary", clustering_method = "ward.D2",
	labCol = NULL, annotation_row = annotation, annotation_col = annotation, cellwidth = 10, cellheight = 10, legend = FALSE, fontsize_col = 8, fontsize_row = 8, cluster_col = F,
	border_color = NA, annotation_colors = list(Module = c(Pollen = "red", Root = "darkblue")), main = "Module stability for Root & Pollen orthologues", cluster_row = T
)
pheatmap(MatrixPollenRoot, color = colorRampPalette(c("white", "black"))(10), clustering_distance_rows = "binary", clustering_distance_cols = "binary", clustering_method = "ward.D2",
	labCol = NULL, annotation_row = annotation, annotation_col = annotation, cellwidth = 10, cellheight = 10, legend = FALSE, fontsize_col = 8, fontsize_row = 8, cluster_col = T,
	border_color = NA, annotation_colors = list(Module = c(Pollen = "red", Root = "darkblue")), main = "Module stability for Root & Pollen orthologues", cluster_row = T
)
dev.off()



write.table(RMatrix, file = "/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Local/Validation/Stability/ModuleStabilityMatrix.txt", quote = F, sep = "\t")




