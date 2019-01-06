
library(igraph)
library(foreach)
library(doMC)
# registerDoMC(2)
library(ape)


go = read.table("/ngsprojects/iwt/ehsab/Prioritization/prioritization/Combined_Network/New/GO_ID.txt",sep = "\t")
go = read.table("/ngsprojects/iwt/ehsab/Prioritization/prioritization/Combined_Network/New/GO_experimental_terms2015.txt",sep = "\t") # Experimental GO terms
go = read.delim("/ngsprojects/iwt/ehsab/Prioritization/prioritization/Combined_Network/New/GO_OCT2015.txt", header = F)
clas = c("C", "F", "P")
k = c(2,3)
go = go[go[, 4] %in% clas[k], ]
colnames(go)[3] = "GO"
go.limited = go[!(go$V2 %in%  rownames(data.frame(sort(table(as.matrix(go$V2)), decreasing=T)[1:20]))),]



## Old GO terms analysis 2012
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
goo = readLines("/ngsprojects/iwt/ehsab/Prioritization/prioritization/Combined_Network/TAIR/GO_terms.txt")
RR = matrix(NA, 1, 5)
for (i in 1:length(goo)){
	RR = rbind(RR, strsplit(goo[i], "\t")[[1]])
}
RR = RR[-1,]
RR = data.frame(RR)
colnames(RR)[1] = "GO"
rownames(RR) = RR$GO
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

## New GO terms analysis 2015
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# goTerms = go[!duplicated(go$V2),]
goMap = go[!(duplicated(go[,c(2:3)])), -1]
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

GO.terms = unique(go[,2])
AA = table(go[,2])
ExcludGO = names(AA[AA < 10 | AA > 2000])

k = 6
condition.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/condition_name.txt", header = T)
NAME = names(table(condition.names))
NAME.spe = NAME[k]
load(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/",NAME.spe,"/Graph_V2_90.RData", sep = ""))
ref = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Reference_SD.txt", header = T)
# paralogous = read.table("~/Gene_duplication/Seed_genes/paralogous_total.txt", header = T)
paralogous = read.table("~/Gene_duplication/Seed_genes/paralogous_NoFilter.txt", header = T)


j = 2858
jpar = 9603

memb.comun = neighborhood(g, 1, j, mode = c("out"))[[1]]
memb.comun.para = neighborhood(g, 1, jpar, mode = c("out"))[[1]]

memb.comun = which(as.numeric(unlist(gene_score_module[1])) == 1)

load("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/development/Permutation/V2/95/G13393_13899_15525.RData")
memb.comun.para = which(as.numeric(unlist(gene_score_module[1])) == 1)

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Anchorpoints GO terms analysis
# ---------------------------------------------------------------------------------------------
Inter1 = intersect(memb.comun.para, paralogous[paralogous$Ref %in% memb.comun, "Anch"])
Inter2 = intersect(memb.comun, paralogous[paralogous$Ref %in% memb.comun.para, "Anch"])
aa = paralogous[paralogous$Ref %in% Inter2, ]
bb = aa[aa$Anch %in% Inter1,]
colnames(bb)[1] = "References"
Inter_name = merge(ref, bb, by = "References")
colnames(Inter_name)[c(1,3)] = c("Ref", "References")
Inter_name = merge(ref, Inter_name, by = "References")
memb.comun = Inter_name$Ref
memb.comun.para = Inter_name$References
# ---------------------------------------------------------------------------------------------




# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Enrichment function
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
GOenrichment = function(go.terms, Threshold, terms, Direction) {
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




# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Reading genes and anchorpoint inside the modules
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
AnchPL = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Anchorpoint/Geneset1_2_53_90.txt.txt")
AnchRL = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Anchorpoint/AT1G50090_AT3G19710_52_90.txt")

members = AnchPL$Leaf
members = AnchPL$Pollen

members = AnchRL$Leaf
members = AnchRL$Root

go.terms = go[go$V1 %in% members, "GO"]

AnchPollenGO = GOenrichment(go[go$V1 %in% AnchPL$Pollen, "GO"], 0.1, goMap, "over")
AnchPollenLeafGO = GOenrichment(go[go$V1 %in% AnchPL$Leaf, "GO"], 0.1, goMap, "over")
AnchRootLeafGO = GOenrichment(go[go$V1 %in% AnchRL$Leaf, "GO"], 0.1, goMap, "over")
AnchRootGO = GOenrichment(go[go$V1 %in% AnchRL$Root, "GO"], 0.1, goMap, "over")

AnchPollenGO$Module = "Pollen"
AnchPollenLeafGO$Module = "PollenLeaf"
AnchRootGO$Module = "Root"
AnchRootLeafGO$Module = "RootLeaf"

AnchAll = rbind(AnchPollenGO, AnchPollenLeafGO, AnchRootGO, AnchRootLeafGO)
AnchAll$Duplication = "WGD"

GOAll = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/GO/SSD.txt")
GOAll$Duplication = "SSD"

write.table(rbind(AnchAll, GOAll), file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/GO/WGD_SSD.txt", sep = "\t", quote = F, row.names = F)

PollenLeaf = as.character(read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/Leaf_90.txt")[,1])
Pollen = as.character(read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/Pollen_90.txt")[,1])

Root = as.character(read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/AT1G50090_Root_90.txt")[,1])
RootLeaf = as.character(read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/AT3G19710_Leaf_90.txt")[,1])
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------





# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# GO enrichment analysis
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
go.terms = go[go$V1 %in% Pollen, "GO"]
PollenGO = GOenrichment(go.terms, 0.1, goMap, "over")
go.terms = go[go$V1 %in% PollenLeaf, "GO"]
PollenLeafGO = GOenrichment(go.terms, 0.1, goMap, "over")


go.terms = go[go$V1 %in% Root, "GO"]
RootGO = GOenrichment(go.terms, 0.1, goMap, "over")
go.terms = go[go$V1 %in% RootLeaf, "GO"]
RootLeafGO = GOenrichment(go.terms, 0.1, goMap, "over")

PollenGO$Module = "Pollen"
PollenLeafGO$Module = "PollenLeaf"
RootGO$Module = "Root"
RootLeafGO$Module = "RootLeaf"

write.table(rbind(PollenGO, PollenLeafGO, RootGO, RootLeafGO), file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/GO/Module.txt", sep = "\t", row.names = F, quote = F)

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------





# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Comparison dot-plot between different groups
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library("reshape")
GOterms = merge(RR, merge(RootGO, LeafGO, by = "GO", all.x = TRUE, all.y = TRUE)[, c("GO", "FDR.over.x", "FDR.over.y")], by = "GO")[, c("X3", "FDR.over.x", "FDR.over.y")]
colnames(GOterms)[2:3] = c("Flower","Leaf")
colnames(GOterms)[2:3] = c("Root","Leaf")
GO<-melt(GOterms)
colnames(GO)<-c("GOcat","Group","Pvalue")

GO$GOcat<- factor(GO$GOcat,levels=rev(unique(GO$GOcat)))

underover<-vector(length = length(GO$Pvalue))
underover[which(GO$Pvalue<0)]<-"Under"
underover[which(GO$Pvalue>=0)]<-"Over"
underover[which(is.na(GO$Pvalue))]<-"No"
GO<-cbind(GO,underover)
colnames(GO)<-c("GOcat","Group","Pvalue","Direction")
GO$Pvalue[which(GO$Pvalue<0)]<- GO$Pvalue[which(GO$Pvalue<0)]*-1
GO$Pvalue[which(GO$Pvalue<1e-10)]<-(1e-10)
GO$Pvalue<-(-log(GO$Pvalue))

png("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/GO_Leaf_Root_90_Anch.png", width = 10, height = 17, units = "in", res = 1000, pointsize=8)
ggplot(GO,aes(y=GOcat,x=Group))+
  geom_point(aes(colour = Group, size = Pvalue)) +
  scale_size_continuous(breaks= c(4.6, 9.2, 23, 46, 69), labels= c("< 0.01", "< 1e-4", "< 1e-10", "< 1e-20", "< 1e-30"), range = c(2, 5)) +
  scale_colour_manual(values = c("blue","darkorange2")) +
  theme_bw() +
  theme(legend.text = element_text(size = 8, face = "bold"),
        axis.text.y = element_text(size = 14,color="black"),#,angle = 45,),
        axis.text.x = element_text(size = 15,color="black",vjust = 0.5),
        axis.title.y = element_blank(),
        axis.title.x = element_blank() )
dev.off()
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------





# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# GO dendrogram
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load("/ngsprojects/iwt/ehsab/Prioritization/prioritization/Combined_Network/New/GODistanceGeneBase.RData")
load("/ngsprojects/iwt/ehsab/Prioritization/prioritization/Combined_Network/New/ExperimentalGODistanceGeneBase.RData") # Experimental GO terms
RRf["Group"] = 1
RRf2["Group"] = 2

RRF = rbind(RRf, RRf2)
RRF[RRF$GO %in% intersect(RRf$GO, RRf2$GO), "Group"] = 3
RRF= RRF[!duplicated(RRF$GO),]
RRF = RRF[RRF$FDR < Threshold,]
RRF = RRF[!(RRF$GO %in% ExcludGO),]
goids = RRF$GO
BB = GO.Dist[rownames(GO.Dist) %in% goids, colnames(GO.Dist) %in% goids]

PLOT = RRF[, c("GO", "Group")]
clus = hclust(as.dist(1-BB), method = "ward")
colors <- c("#006400", "#0000FF", "#BB0000") # define the colors for the labels
mapColors <- merge(PLOT, data.frame("Group" = c(1:3), "col" = colors), by = "Group", all.x = TRUE)
colorsLabels <- mapColors$col[match(clus$labels, mapColors$GO)]

jpeg(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/12450_13393_90_GODendrogram.jpeg", sep = ""), width = 11, height = 7, units = "in", res = 200, pointsize=8)
plot(as.phylo(clus), cex = 0.6, no.margin = TRUE, tip.color = as.vector(colorsLabels), font = 8)
dev.off()

for(i in 1:nrow(BB)){
	colnames(BB)[i] = rownames(BB)[i] = as.character(RRF[RRF$GO == colnames(BB)[i],"X3"])
}

clus = hclust(as.dist(1-BB), method = "ward")
mapColors <- merge(RRF[, c("X3", "Group")], data.frame("Group" = c(1:3), "col" = colors), by = "Group", all.x = TRUE)
colorsLabels <- mapColors$col[match(clus$labels, mapColors$X3)]

jpeg(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/12450_13393_95_GODendrogram_Terms_Anchor.jpeg", sep = ""), width = 11, height = 7, units = "in", res = 800, pointsize=8)
plot(as.phylo(clus), cex = 0.2, no.margin = TRUE, tip.color = as.vector(colorsLabels), font = 8)
dev.off()

write.table(RRF, file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/12450_13393_90_GO.txt", row.names = F, quote = F, sep = "\t")


for(i in 1:GO.ID) {
	AA = unique(as.character(GO.exp[as.character(GO.exp$V6) %in% GO.ID[i], 1]))

for(i in 1:GO.ID) {
	sapply(1:length(GO.ID), function(j) sum(abs(AA[,i] - AA[,j]))/sum(AA[,i] + AA[,j]))


dFun <- function(x,y){length(which(x%in%y))/min(length(x),length(y))}
outer(list.df, list.df, Vectorize(dFun))
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
