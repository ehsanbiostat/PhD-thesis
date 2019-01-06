library(igraph)
library("isa2")
library("foreach")
library("doMC")
library(igarph)
library(igraph)
data.exp = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated CORNET2.0/Integrated_SD.txt", colClasses = "numeric", nrows = 19285, comment.char = "", header = T)
library(igraph)
library(foreach)
library(doMC)
library(seriation)
options(width = 200)
paralogous = read.table("~/Gene_duplication/Seed_genes/paralogous_NoFilter.txt", header = T)
PARA = unique(paralogous$Ref)
GeneScore = foreach (i = PARA, .combine = cbind) %do% {
load(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/development/Raw/G",i,".RData", sep = ""))
# gene_scores = cbind(gene_scores, 1:19285)
# gene_scores = gene_scores[order(-gene_scores[,1]),]
gene_scores
}
COR = cor(GeneScore[, 1:500])
library(corrplot)
?corrplot
corrMatOrder(COR, order = "AOE")
ORDER = corrMatOrder(COR, order = "AOE")
source("/ngsprojects/iwt/ehsab/Gene_duplication/Scripts/General/myImagePlot.R")
myImagePlot(COR[ORDER,ORDER])
ORDER = corrMatOrder(COR, order = "FPE")
ORDER = corrMatOrder(COR, order = "FPC")
myImagePlot(COR[ORDER,ORDER])
ORDER = corrMatOrder(COR, order = "hclust")
myImagePlot(COR[ORDER,ORDER])
COR = cor(GeneScore)
Order.AOE = corrMatOrder(COR, order = "AOE")
Order.FPC = corrMatOrder(COR, order = "FPC")
Order.hclust = corrMatOrder(COR, order = "hclust")
myImagePlot(COR[Order.AOE,Order.AOE])
jpeg(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Validation/Compendium_wise/development/Heatmap/AnchHclust", sep = ""), width = 11, height = 7, units = "in", res = 300, pointsize=8)
myImagePlot(COR[Order.hclust,Order.hclust])
dev.off()
Order.hclust[1:100]
PARA[Order.hclust[1:100]]
jpeg(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Validation/Compendium_wise/development/Heatmap/AnchFPC", sep = ""), width = 11, height = 7, units = "in", res = 300, pointsize=8)
myImagePlot(COR[Order.FPC,Order.FPC])
dev.off()
Order.FPC[1:10]
head(Order.FPC)
tail(Order.FPC)
paralogous[head(Order.FPC)]
paralogous[head(Order.FPC),]
paralogous[PARA[head(Order.FPC)],]
PARA[head(Order.FPC)]
PARA[tail(Order.FPC)]
paralogous[paralogous$Ref %in% PARA[head(Order.FPC)],]
hclust(COR)
dim(COR)
table(is.na(COR)))
table(is.na(COR))
?hclust
hclust(as.dist(COR))
plot(hclust(as.dist(COR)))
?kmeans
KS = kmeans(COR, centers = 10)
?kmeans
KS$cluster
table(KS$cluster)
KS$cluster[KS$cluster == 1]
Geneset1 = c(10182,1352,2663,6258,1932,13393,7682,13899,8399,15525,13906,6874,9603,2228,39,1694,10546,6713)
Geneset2 = c(5792,15713,399,16681,17367,17229,2858,16546,7903,14025,440,18332,12364,10298,12450)
which(PARA[Order.hclust] %in% Geneset1)
which(PARA[Order.hclust] %in% Geneset2)
PARA[Order.hclust][490]
COR.order = COR[Order.hclust, Order.hclust]
myImagePlot(COR.order[which(PARA[Order.hclust] %in% Geneset1), which(PARA[Order.hclust] %in% Geneset2)])
Geneset
Geneset = c(Geneset1, Geneset2)
myImagePlot(COR.order[which(PARA[Order.hclust] %in% Geneset), which(PARA[Order.hclust] %in% Geneset)])
which(PARA[Order.hclust] %in% Geneset)
colnames(COR)[1:10]
ref = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Reference_SD.txt", header = T)
ref[PARA[1:10],1]
colnames(COR) = rownames(COR) = PARA
COR.order = COR[Order.hclust, Order.hclust]
myImagePlot(COR.order[which(PARA[Order.hclust] %in% Geneset), which(PARA[Order.hclust] %in% Geneset)])
myImagePlot(COR.order[which(rownames(COR.order) %in% Geneset), which(colnames(COR.order) %in% Geneset)])
myImagePlot(COR.order[490:650, 490:650])
myImagePlot(COR.order[3665:4000, 3665:4000])
myImagePlot(COR.order[1:1000, 1:1000])
myImagePlot(COR.order[500:1000, 500:1000])
myImagePlot(COR.order[500:1500, 500:1500])
myImagePlot(COR.order[750:1100, 750:1100])
myImagePlot(COR.order[500:1100, 500:1100])
myImagePlot(COR.order[550:1100, 550:1100])
myImagePlot(COR.order[550:1000, 550:1000])
colnames(COR.order)[550:1000]
myImagePlot(COR.order[550:750, 550:750])
myImagePlot(COR.order[560:750, 560:750])
myImagePlot(COR.order[560:850, 560:850])
myImagePlot(COR.order[560:800, 560:800])
myImagePlot(COR.order[560:770, 560:770])
myImagePlot(COR.order[560:780, 560:780])
myImagePlot(COR.order[560:771, 560:771])
myImagePlot(COR.order[560:772, 560:772])
myImagePlot(COR.order[560:773, 560:773])
colnames(COR.order)[560:772]
AA = as.numeric(colnames(COR.order)[773:100])
AA
BB = as.numeric(colnames(COR.order)[773:110])
BB
myImagePlot(COR.order[773:1000, 773:1000])
myImagePlot(COR.order[773:1000, 50:772])
myImagePlot(COR.order[773:1000, 560:772])
BB = as.numeric(colnames(COR.order)[773:1000])
BB
AA
memb.comun = data.frame("Ref" = AA)
memb.comun.para = merge(memb.comun, paralogous, by = "Ref")
memb.comun.para
length(intersect(BB, memb.comun.para$Anch))
BB
AA
memb.comun.para$Anch
(intersect(BB, memb.comun.para$Anch))
AA
BB
write.table(AA, file = "/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/development/GeneScoreCorrelationSet1", sep = "\t", row.names = F, col.names= F)
write.table(BB, file = "/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/development/GeneScoreCorrelationSet2", sep = "\t", row.names = F, col.names= F)
myImagePlot(COR.order[2500:3500, 2500:3500])
myImagePlot(COR.order[2750:3500, 2750:3500])
myImagePlot(COR.order[2800:3500, 2800:3500])
myImagePlot(COR.order[2850:3500, 2850:3500])
myImagePlot(COR.order[2825:3600, 2825:300])
myImagePlot(COR.order[2830:3500, 2830:3500])
myImagePlot(COR.order[2835:3500, 2835:3500])
myImagePlot(COR.order[2835:3600, 2835:3600])
myImagePlot(COR.order[2835:3900, 2835:3900])
myImagePlot(COR.order[2835:4200, 2835:4200])
myImagePlot(COR.order[2835:4000, 2835:4000])
myImagePlot(COR.order[2835:3980, 2835:3980])
myImagePlot(COR.order[2835:3100, 2835:3100])
myImagePlot(COR.order[2835:3200, 2835:3200])
myImagePlot(COR.order[2835:3300, 2835:3300])
myImagePlot(COR.order[2835:3500, 2835:3500])
myImagePlot(COR.order[2835:3400, 2835:3400])
myImagePlot(COR.order[2835:3350, 2835:3350])
myImagePlot(COR.order[2835:3380, 2835:3380])
myImagePlot(COR.order[2835:3370, 2835:3370])
myImagePlot(COR.order[2835:3365, 2835:3365])
myImagePlot(COR.order[2835:3360, 2835:3360])
myImagePlot(COR.order[2835:3355, 2835:3355])
colnames(COR.order)[2835:3355]
CC = as.numeric(colnames(COR.order)[2835:3355])
CC
DD = as.numeric(colnames(COR.order)[3356:3980])
DD
memb.comun = data.frame("Ref" = CC)
memb.comun.para = merge(memb.comun, paralogous, by = "Ref")
length(intersect(DD, memb.comun.para$Anch))
write.table(CC, file = "/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/development/GeneScoreCorrelationSet11", sep = "\t", row.names = F, col.names= F)
write.table(DD, file = "/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/development/GeneScoreCorrelationSet22", sep = "\t", row.names = F, col.names= F)
savehistory(file="/ngsprojects/iwt/ehsab/Gene_duplication/Scripts/SA/Merge_module/GenescoreCorrplot.R")
















function(edgelist, vertices, sorted = FALSE, decreasing = FALSE, 
           ordering = NULL, labels = NULL)
  {
    # ======================================================
    # Checking arguments
    # ======================================================
    # edgelist as a two-column matrix
    if (!is.matrix(edgelist) || ncol(edgelist) != 2)
      stop("\nSorry, 'edgelist' must be a two column matrix")
    
    num_edges = nrow(edgelist)
    # get nodes (this could be numeric or character)
    if(hasArg(vertices)){
      #to deal with singleton nodes
      nodes = vertices 
    }else{
      nodes = unique(as.vector(t(edgelist)))  
    }
    num_nodes = length(nodes)
    # check labels (i.e. node names)
    if (!is.null(labels))
    {
      if (length(labels) != num_nodes)
        stop("\nLength of 'labels' differs from number of nodes")
    } else {
      labels = nodes
    }
    
    # auxiliar order (this may change if sorted or ordering required)
    aux_ord = 1:num_nodes  
    
    # If sorted is required, ennumerate nodes
    if (sorted) {
      ordered_nodes = order(nodes, decreasing = decreasing)
      nodes = nodes[ordered_nodes]
      labels = labels[ordered_nodes]
      # auxiliar order
      aux_ord = ordered_nodes
    }
    
    # If ordering is provided, re-ennumerate nodes
    if (!is.null(ordering)) 
    {
      if (length(ordering) != num_nodes) {
        stop("\nLength of 'ordering' differs from number of nodes")      
      }
      
      if (is.character(ordering)) {
        # make sure labels contains elements in ordering
        unmatched_ordering <- !(ordering %in% labels)
        if (any(unmatched_ordering)) {
          undetected = ordering[unmatched_ordering]
          stop(sprintf("\nUnrecognized values in ordering: '%s'", undetected))
        }
        ordering = match(ordering, labels)
      }
      
      nodes = nodes[ordering]
      labels = labels[ordering]
      # auxiliar order
      aux_ord = ordering
    }
    
    ## output
    list(
      nodes = nodes,
      labels = labels,
      num_nodes = num_nodes,
      num_edges = num_edges,
      aux_ord = aux_ord
    )
  }
