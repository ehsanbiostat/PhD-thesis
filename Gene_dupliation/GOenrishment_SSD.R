


# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# GO analysis
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(igraph)
library(foreach)
library(doMC)
# registerDoMC(2)
library(ape)
options(width = 200)

go = read.delim("/ngsprojects/iwt/ehsab/Prioritization/prioritization/Combined_Network/New/GO_OCT2015.txt", header = F)
clas = c("C", "F", "P")
k = c(2,3)
go = go[go[, 4] %in% clas[k], ]
colnames(go)[3] = "GO"
go.limited = go[!(go$V2 %in%  rownames(data.frame(sort(table(as.matrix(go$V2)), decreasing=T)[1:20]))),]
goMap = go[!(duplicated(go[,c(2:3)])), -1]
# go = go[go[, 4] == clas[k] & (go[, 5] %in% c("EXP", "IDA","IPI", "IMP", "IGI", "IEP")), ] # only experimental evidence codes


# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------



# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Enrichment function
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
GOenrichment = function(go.terms, terms, Threshold, Direction) {
	if(!(Direction %in% c("both", "over", "under"))) {
		cat("Direction has to be one of 'both', 'over' or 'under'.\n")
		cat("Direction default value has been setup to 'both'.\n")
		Direction = "both"
	}
	freq = table(matrix(go.terms)) # Frequency of GO among
	pvalue.under = pvalue.over = odd = global = loc = c()

	if (length(go.terms) > 0){
		for (i in 1:length(freq)){
			z = matrix(c(length(go.terms) - freq[i], freq[i]), 2)
			x = matrix(table(go[,"GO"] == names(freq[i])))
			x = cbind(x ,z) 
			pvalue.over[i] = fisher.test(x, alternative = "greater")$p
			pvalue.under[i] = fisher.test(x, alternative = "less")$p
			global[i] = x[2,1]
			loc[i] = x[2,2]
			odd[i] = (x[1,1] * x[2,2])/(x[1,2] * x[2,1])
		}
		
		cat("First modules finished\n")
		names(pvalue.under) = names(pvalue.over) = names(odd) = names(global) = names(loc) = names(freq)
		f = data.frame(pvalue.over, pvalue.under, odd, global, loc)
		f["FDR.over"] = p.adjust(f$pvalue.over, method = c("fdr"))
		f["FDR.under"] = p.adjust(f$pvalue.under, method = c("fdr"))
		f = f[order(f$FDR.over), ]
		top = c(rownames(f[f$FDR.over < Threshold,]))
		
		if (Direction == "both") f = f[f$FDR.over < Threshold | f$FDR.under < Threshold,]
		if (Direction == "over") f = f[f$FDR.over < Threshold,]
		if (Direction == "under") f = f[f$FDR.under < Threshold,]

		f["GO"] = rownames(f)
		RRf = merge(terms, f, by = "GO")
		if (Direction == "both") RRf = RRf[order(RRf$FDR.over, decreasing = F), ]
		if (Direction == "over") RRf = RRf[order(RRf$FDR.over, decreasing = F), ]
		if (Direction == "under") RRf = RRf[order(RRf$FDR.under, decreasing = F), ]

	} else {
			RRf = data.frame(NA)
	}	
	
	return(RRf)
}
# ---------------------------------------------------------------------------------------------


library(igraph)
library(foreach)
library(doMC)

ref = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Reference_SD.txt", header = T)
paralogous = read.table("~/Gene_duplication/Seed_genes/paralogous_NoFilter.txt", header = T)
GeneFamilyDup = read.delim("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/GeneFamilyZhenFeb16SSDLSD.txt")
GeneFamily = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/GeneFamilyZhenFeb16.txt")

Pollen = as.character(read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/Pollen_90.txt")[,1])
PollenLeaf = as.character(read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/Leaf_90.txt")[,1])
Root = as.character(read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/AT1G50090_Root_90.txt")[,1])
RootLeaf = as.character(read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/AT3G19710_Leaf_90.txt")[,1])

DupRootLeaf = GeneFamilyDup[GeneFamilyDup$Var1 %in% RootLeaf & GeneFamilyDup$Var2 %in% Root, ]
DupRootLeaf = DupRootLeaf[!duplicated(DupRootLeaf[,1:2]),]
colnames(DupRootLeaf)[1:2] = c("Leaf", "Root")
SSDRootLeaf = unique(c(as.matrix(DupRootLeaf[DupRootLeaf$Duplication == "SSD", "Leaf"])))
SSDRoot = unique(c(as.matrix(DupRootLeaf[DupRootLeaf$Duplication == "SSD", "Root"])))

go.termsRoot = go[go$V1 %in% SSDRoot, "GO"]
go.termsRootLeaf = go[go$V1 %in% SSDRootLeaf, "GO"]
GORootLeaf = GOenrichment(go.termsRootLeaf, goMap, 0.1, "both")
GORoot = GOenrichment(go.termsRoot, goMap, 0.1, "both")

write.table(DupRootLeaf, file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/DupRootLeaf.txt", row.names = F, col.names = T, quote = F, sep = "\t")
write.table(GeneFamily[GeneFamily$V2 %in% SSDRootLeaf,], file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/SSDRootLeafGenes.txt", row.names = F, col.names = F, quote = F, sep = "\t")
write.table(GeneFamily[GeneFamily$V2 %in% SSDRoot,], file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/SSDRootGenes.txt", row.names = F, col.names = F, quote = F, sep = "\t")


DupPollenLeaf = GeneFamilyDup[GeneFamilyDup$Var1 %in% PollenLeaf & GeneFamilyDup$Var2 %in% Pollen, ]
DupPollenLeaf = DupPollenLeaf[!duplicated(DupPollenLeaf[,1:2]),]
colnames(DupPollenLeaf)[1:2] = c("Leaf", "Pollen")
SSDPollenLeaf = unique(c(as.matrix(DupPollenLeaf[DupPollenLeaf$Duplication == "SSD", "Leaf"])))
SSDPollen = unique(c(as.matrix(DupPollenLeaf[DupPollenLeaf$Duplication == "SSD", "Pollen"])))
go.termsPollenLeaf = go[go$V1 %in% SSDPollenLeaf, "GO"]
go.termsPollen = go[go$V1 %in% SSDPollen, "GO"]
GOPollenLeaf = GOenrichment(go.termsPollenLeaf, goMap, 0.1, "both")
GOPollen = GOenrichment(go.termsPollen, goMap, 0.1, "both")
DuplicatePollenLeaf = unique(c(as.matrix(DupPollenLeaf[, "Leaf"])))
DupPollen = unique(c(as.matrix(DupPollenLeaf[, "Pollen"])))

write.table(DupPollenLeaf, file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/DupPollenLeaf.txt", row.names = F, col.names = T, quote = F, sep = "\t")
write.table(GeneFamily[GeneFamily$V2 %in% SSDPollenLeaf,], file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/SSDPollenLeafGenes.txt", row.names = F, col.names = F, quote = F, sep = "\t")
write.table(GeneFamily[GeneFamily$V2 %in% SSDPollen,], file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/SSDPollenGenes.txt", row.names = F, col.names = F, quote = F, sep = "\t")

Duplicates = unique(c(as.character(GeneFamilyDup$Var1), as.character(GeneFamilyDup$Var2)))
write.table(Duplicates, file = "/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/AllDuplicateGenes_ZhenGeneFamily.txt", row.names = F, col.names = F, quote = F, sep = "\t")



GORootLeaf$Module = "RootLeaf"
GORoot$Module = "Root"
GOPollenLeaf$Module = "PollenLeaf"
GOPollen$Module = "Pollen"

GOAll = rbind(GORootLeaf, GORoot, GOPollenLeaf, GOPollen)

write.table(GOAll, file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/GO/SSD.txt", sep = "\t", row.names = F, quote = F)


















