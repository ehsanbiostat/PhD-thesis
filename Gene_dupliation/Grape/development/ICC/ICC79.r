



library(foreach)
library(doMC)
registerDoMC(3)
options(width = 200)

# Correlation between two rows of two different matrices
source("/ngsprojects/iwt/ehsab/Gene_duplication/Scripts/Functions/CORR.R")



data.grape = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Vitis_vinifera/GE.mean.SD.txt", nrows = 23629, comment.char = "", header = T)
gene.names <- read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Vitis_vinifera/gene.names.txt")
gene.names = as.matrix(gene.names)

condition.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/condition_name.txt", header = T)
NAME = names(table(condition.names))

## Arabidopsis
data.arab = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Integrated_SD.txt", colClasses = "numeric", nrows = 19285, comment.char = "", header = T)
data.arab = as.matrix(data.arab[, condition.names == "development"])
Ref = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Reference_SD.txt", header = T)
rownames(data.arab) = Ref[,1]

###############################################################################################################################################################
# Remove non-uniqued mapped PLAZA ID to gene names and the other way around
###############################################################################################################################################################
map = read.delim("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Orthology/Grape-Arabidopsis/Map.12x.plaza.txt")
# selection = which(sapply(1:sum(duplicated(map$name)), function(i){length(unique(map[map$name == map$name[duplicated(map$name)][i],"PLAZA.ID"]))}) > 1)
# map = map[!(map$name %in% map$name[duplicated(map$name)][selection]),]
# selection = which(sapply(1:sum(duplicated(map$PLAZA.ID)), function(i){length(unique(map[map$PLAZA.ID == map$PLAZA.ID[duplicated(map$PLAZA.ID)][i],"name"]))}) > 1)
# map = map[!(map$PLAZA.ID %in% map$PLAZA.ID[duplicated(map$PLAZA.ID)][selection]),]
map = map[!(duplicated(map[, c("name", "PLAZA.ID")])),]

###############################################################################################################################################################

Ref = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Reference_SD.txt", header = T)
rownames(Ref) = Ref[,1]
onetone = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Orthology/Grape-Arabidopsis/ZhenProcessed.txt", header = T)
onetone = onetone[!(duplicated(onetone[, "PLAZA.ID"]) | duplicated(onetone[, "gene_id.y"])),]

source("/ngsprojects/iwt/ehsab/Gene_duplication/Scripts/Functions/CORR.R")
Ref = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Reference_SD.txt", header = T)
rownames(Ref) = Ref[,1]


LeafRoot = read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/AT3G19710_Leaf_90.txt")
Root = read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/AT1G50090_Root_90.txt")
AnchRoot = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Anchorpoint/AT1G50090_AT3G19710_52_90.txt")



Leaf = read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/Leaf_90.txt")
Pollen = read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/Pollen_90.txt")
Anch = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Anchorpoint/Geneset1_2_53_90.txt.txt")

map = read.delim("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Orthology/Grape-Arabidopsis/Map.12x.plaza.txt")
map = map[!(duplicated(map[, c("name", "PLAZA.ID")])),] # Map PLAZA ID to the gene names in the grapevine expression data
onetone = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Orthology/Grape-Arabidopsis/ZhenProcessed.txt", header = T)
maponetone = merge(map, onetone, by = "PLAZA.ID")
maponetone = maponetone[maponetone$gene_id.y %in% Ref$Gene,] # Remove arabidopsis genes which are not in the compendium
maponetone = maponetone[maponetone$name %in% rownames(data.grape),] # Remoe grapevine genes which are not in the development compendium

k = 79
load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Permutation/10Fold/Folds.RData")
# RR = foreach(k = 1:length(flds)) %dopar% {
	print(k)
	COR.grape = cor(t(data.grape[as.character(maponetone$name), flds[k][[1]]]))
	COR.arab = cor(t(data.arab[as.character(maponetone$gene_id.y),]))
	onetoneCorr = CORR(COR.grape, COR.arab, w = runif(nrow(COR.grape)))

	threshold = 1
	i = 0

	while (threshold > 0.1) {
		i = i +1
		# print(i)
		Exclude = which(onetoneCorr < 0)
		COR.grape = COR.grape[-Exclude,-Exclude]
		COR.arab = COR.arab[-Exclude,-Exclude]
		onetoneCorrold = onetoneCorr[-Exclude]
		onetoneCorr = CORR(COR.grape, COR.arab, w = onetoneCorrold)
		threshold = sum((onetoneCorr - onetoneCorrold)^2)
	}

	save(COR.grape, COR.arab, onetoneCorr, file = paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Permutation/10Fold/",k,".RData", sep = ""))

# }

