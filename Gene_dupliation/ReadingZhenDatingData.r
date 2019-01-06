
library(foreach)

a = read.delim("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/genefamily.noncore.dup.txt", header = F)
b = read.delim("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/genefamily.noncore.dup.txt", header = F)
a = rbind(a, b)

# Find Arabidopsis Genes
FoundArab = function(Row, Col){sum(unlist(lapply(strsplit(as.character(a[Row, Col]), " "), strsplit, "[|]")) == "Atha")}

Result1 = mapply(FoundArab, Row = 1:nrow(a), Col = 4)
Result2 = mapply(FoundArab, Row = 1:nrow(a), Col = 5)
Result = cbind(Result1, Result2)


DupliNumber = 1

OneToMany = foreach(m = 0:9, .combine = rbind) %dopar% {
	print(m)
	Index = which(Result[,1] == DupliNumber + m & Result[,2] == DupliNumber | Result[,1] == DupliNumber & Result[,2] == DupliNumber + m)
	OneTo = foreach(i = Index, .combine = rbind) %do% {
		Q = unlist(lapply(strsplit(as.character(a[i, 4]), " "), strsplit, "[|]"))
		QQ = unlist(lapply(strsplit(as.character(a[i, 5]), " "), strsplit, "[|]"))
		data.frame(a[i, 1:3], "Gene" = Q[which(Q == "Atha") + 1], "Paralogous" = QQ[which(QQ == "Atha") + 1], "Type" = DupliNumber + m)
	}
	OneTo
}


Dup = foreach(DupliNumber = 1:9, .combine = rbind) %dopar% {
	print(DupliNumber)
	TO = foreach(m = 0:9, .combine = rbind) %do% {
		Index = which(Result[,1] == DupliNumber + m & Result[,2] == DupliNumber | Result[,1] == DupliNumber & Result[,2] == DupliNumber + m)
		OneTo = foreach(i = Index, .combine = rbind) %do% {
			Q = unlist(lapply(strsplit(as.character(a[i, 4]), " "), strsplit, "[|]"))
			QQ = unlist(lapply(strsplit(as.character(a[i, 5]), " "), strsplit, "[|]"))
			T = Q[which(Q == "Atha") + 1]
			TT = QQ[which(QQ == "Atha") + 1]
			data.frame(a[i, 1:3], "Gene" = rep(T, length(TT)), "Paralogous" = rep(TT, length(T)), "Type" = paste(DupliNumber, "TO", DupliNumber + m, sep = "."))
		}
		if(length(OneTo) > 0) OneTo
	}
	TO
}

write.table(Dup, file = "/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/ZhenDuplicationEvent.txt", row.names = F, sep = "\t", quote = F)

ZhenDup[(ZhenDup$Gene == as.character(Inter_name[i,2]) & ZhenDup$Paralogous == as.character(Inter_name[i,4]) | ZhenDup$Gene == as.character(Inter_name[i,4]) & ZhenDup$Paralogous == as.character(Inter_name[i,2])),]



################################################################################################################################################################################################
#
#
################################################################################################################################################################################################

zhen = read.delim("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/GeneFamilyZhenFeb16.txt", header = F)
Freq = table(zhen$V1)
genefamily = names(Freq)[Freq > 1]

GeneFamily = foreach (i = 1:length(genefamily), .combine = rbind) %do% {
	print(i)
	a = t(combn(as.character(zhen[zhen$V1 %in% genefamily[i], "V2"]), 2))
	data.frame(rbind(a, cbind(a[,2], a[,1])), "GeneFamily" = genefamily[i])
}

head(GeneFamily)


################################################################################################################################################################################################

library(ggplot2)


data = data.frame("y" = Freq[Freq > 1], "x" = "GeneFamily", "Group" = 0)
ggplot(data, aes(x = x)) + geom_histogram(binwidth = 10) + xlab("Gene family size")
ggplot(data, aes(x = x, y = y)) + geom_violin() + geom_jitter(height = 0) + ylab("Gene family size")


GeneFamilyDup = read.delim("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/GeneFamilyZhenFeb16SSDLSD.txt")
Pollen = as.character(read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/Pollen_90.txt")[,1])
PollenLeaf = as.character(read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/Leaf_90.txt")[,1])
Root = as.character(read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/AT1G50090_Root_90.txt")[,1])
RootLeaf = as.character(read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/AT3G19710_Leaf_90.txt")[,1])

DupRootLeaf = GeneFamilyDup[GeneFamilyDup$Var1 %in% RootLeaf & GeneFamilyDup$Var2 %in% Root, ]
DupRootLeaf = DupRootLeaf[!duplicated(DupRootLeaf[,1:2]),]
colnames(DupRootLeaf)[1:2] = c("Leaf", "Root")

DupPollenLeaf = GeneFamilyDup[GeneFamilyDup$Var1 %in% PollenLeaf & GeneFamilyDup$Var2 %in% Pollen, ]
DupPollenLeaf = DupPollenLeaf[!duplicated(DupPollenLeaf[,1:2]),]
colnames(DupPollenLeaf)[1:2] = c("Leaf", "Pollen")


GeneFamilyRootLeaf = unique(as.character(DupRootLeaf[, "GeneFamily"]))
GeneFamilyPollenLeaf = unique(as.character(DupPollenLeaf[, "GeneFamily"]))
SSDGeneFamily = unique(as.character(DupRootLeaf[DupRootLeaf$Duplication == "SSD", "GeneFamily"]))
LSDGeneFamily = unique(as.character(DupRootLeaf[DupRootLeaf$Duplication == "LSD", "GeneFamily"]))
LSDSSDGeneFamily = intersect(unique(as.character(DupRootLeaf[DupRootLeaf$Duplication == "SSD", "GeneFamily"])), unique(as.character(DupRootLeaf[DupRootLeaf$Duplication == "LSD", "GeneFamily"])))



data$Group[rownames(data) %in% SSDGeneFamily] = "Root_Leaf_SSD"
data$Group[rownames(data) %in% LSDGeneFamily] = "Root_Leaf_LSD"
data$Group[rownames(data) %in% LSDSSDGeneFamily] = "Root_Leaf_LSDSSD"
data$Group[rownames(data) %in% GeneFamilyRootLeaf] = "Root_Leaf"
data$Group[rownames(data) %in% GeneFamilyPollenLeaf] = "Pollen_Leaf"
data$Group[data$Group == "0"] = "All"

ggplot(data, aes(x = x, y = y), color = Group) + geom_violin() + geom_jitter(height = 20, aes(color = data$Group)) + ylab("Gene family size")

ggplot(data, aes(x = y, color = Group)) + geom_density() 


pdf("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/SSD/GeneFamilySize_RootLeaf.pdf", width = 16, height= 8)
ggplot(data, aes(x = x, y = y), color = Group) + geom_violin() + geom_jitter(height = 10, aes(color = data$Group), size = 1.2) + ylab("Gene family size")
dev.off()



