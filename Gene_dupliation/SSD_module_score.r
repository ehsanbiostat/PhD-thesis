/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   Anch_module_score.r                                :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: ehsab <ehsab@student.42.fr>                +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2015/03/23 18:28:09 by ehsab             #+#    #+#             */
/*   Updated: 2015/08/25 11:28:57 by ehsab            ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */




# library(fields)
library(igraph)
library(foreach)
library(doMC)
# registerDoMC(5)

condition.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/condition_name.txt", header = T)
NAME = names(table(condition.names))
NAME = NAME[-4]
# paralogous = read.table("~/Gene_duplication/Seed_genes/paralogous_total.txt", header = T)
ref = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Reference_SD.txt", header = T)
paralogous = read.table("~/Gene_duplication/Seed_genes/paralogous_NoFilter.txt", header = T)
GeneFamilyDup = read.delim("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/GeneFamilyZhenFeb16SSDLSD.txt")
genefamilysize = read.delim("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/GeneFamilySize_Zhen.txt")


Pollen = as.character(read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/Pollen_90.txt")[,1])
PollenLeaf = as.character(read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/Leaf_90.txt")[,1])
Root = as.character(read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/AT1G50090_Root_90.txt")[,1])
RootLeaf = as.character(read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/AT3G19710_Leaf_90.txt")[,1])
DupRootLeaf = GeneFamilyDup[GeneFamilyDup$Var1 %in% RootLeaf & GeneFamilyDup$Var2 %in% Root, ]
DupRootLeaf = DupRootLeaf[!duplicated(DupRootLeaf[,1:2]),]
DupPollenLeaf = GeneFamilyDup[GeneFamilyDup$Var1 %in% PollenLeaf & GeneFamilyDup$Var2 %in% Pollen, ]
DupPollenLeaf = DupPollenLeaf[!duplicated(DupPollenLeaf[,1:2]),]
DupRootLeafSize = genefamilysize[unique(as.character(DupRootLeaf$GeneFamily)),1]
DupPollenLeafSize = genefamilysize[unique(as.character(DupPollenLeaf$GeneFamily)),1]


k = 5

# Score.total = foreach (k = 1:length(NAME), .combine = rbind) %dopar% {
	# print(k)
	load(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/",NAME[k],"/Graph_V2_90.RData", sep = ""))
	Score = matrix(NA, nrow(paralogous), 17)
	colnames(Score) = c("Gene", "Paralogous", "Paralogous.cluster", "Gene.cluster.size", "Paralogous.cluster.size",
		"Cluster.overlap", "Cluster.overlap.Anch",  "Cluster.withi", "Paralogous.withi", "SSDpairs", "LSDpairs", "GeneFamilyNumber", "SSDGeneFamily",
		"LSDGeneFamily", "LSDSSDGeneFamily", "GeneFamilySizeRoot", "GeneFamilySizePollen")
	Score = as.data.frame(Score)

	for (m in 1:nrow(paralogous)) {
		print(m)
		i = paralogous[m,1]
		j = paralogous[m,2]
		memb.comun = data.frame("Ref" = neighborhood(g, 1, i, mode = c("out"))[[1]])
		memb.comun.para = merge(memb.comun, paralogous, by = "Ref")
		memb.comun.para.size = nrow(memb.comun.para)

		a = as.character(ref[neighborhood(g, 1, i, mode = c("out"))[[1]],1])
		A = as.character(ref[neighborhood(g, 1, j, mode = c("out"))[[1]],1])
		DupAa = GeneFamilyDup[GeneFamilyDup$Var1 %in% a & GeneFamilyDup$Var2 %in% A, ]
		Score[m, "SSDpairs"] = sum(DupAa$Duplication == "SSD")
		Score[m, "LSDpairs"] = sum(DupAa$Duplication == "LSD")
		GeneFamilyAa = unique(as.character(DupAa[, "GeneFamily"]))
		SSDGeneFamily = unique(as.character(DupAa[DupAa$Duplication == "SSD", "GeneFamily"]))
		LSDGeneFamily = unique(as.character(DupAa[DupAa$Duplication == "LSD", "GeneFamily"]))
		LSDSSDGeneFamily = intersect(unique(as.character(DupAa[DupAa$Duplication == "SSD", "GeneFamily"])), unique(as.character(DupAa[DupAa$Duplication == "LSD", "GeneFamily"])))

		GSizeTest = genefamilysize[names(table(as.character(DupAa$GeneFamily))),1]
		Score[m, "GeneFamilySizeRoot"] = ks.test(DupRootLeafSize, GSizeTest)$p.value
		Score[m, "GeneFamilySizePollen"] = ks.test(DupPollenLeafSize, GSizeTest)$p.value

		Score[m, "GeneFamilyNumber"] = length(GeneFamilyAa)
		Score[m, "SSDGeneFamily"] = length(SSDGeneFamily)
		Score[m, "LSDGeneFamily"] = length(LSDGeneFamily)
		Score[m, "LSDSSDGeneFamily"] = length(LSDSSDGeneFamily)


		Score[m, "Gene"] = i
		Score[m, "Paralogous"] = j
		Score[m, "Paralogous.cluster"] = length(intersect(neighborhood(g, 1, j, mode = c("out"))[[1]], memb.comun.para$Anch))
		Score[m, "Paralogous.withi"] = nrow(paralogous[paralogous$Ref %in% neighborhood(g, 1, j, mode = c("out"))[[1]],])
		Score[m, "Cluster.withi"] = nrow(paralogous[paralogous$Ref %in% neighborhood(g, 1, i, mode = c("out"))[[1]],])
		Score[m, "Cluster.overlap"] = length(intersect(neighborhood(g, 1, j, mode = c("out"))[[1]], neighborhood(g, 1, i, mode = c("out"))[[1]]))
		Score[m, "Gene.cluster.size"] = dim(memb.comun)[1]
		Score[m, "Paralogous.cluster.size"] = length(neighborhood(g, 1, j, mode = c("out"))[[1]])
		Score[m, "Cluster.overlap.Anch"] = length(intersect(intersect(neighborhood(g, 1, j, mode = c("out"))[[1]], memb.comun.para$Anch), 
			intersect(neighborhood(g, 1, j, mode = c("out"))[[1]], neighborhood(g, 1, i, mode = c("out"))[[1]])))
	}

	
write.table(Score, file = paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Local/Compendium_wise/",NAME[k],"/Score/Anch_Score90_SSD.txt", sep = ""), row.names = F, quote = F, sep = "\t")




library(ggplot2)
library(gridExtra)
Score = read.delim(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Local/Compendium_wise/",NAME[k],"/Score/Anch_Score90_SSD.txt", sep = ""))
data = data.frame("Conserved" = (Score$Cluster.overlap/(Score$Gene.cluster.size + Score$Paralogous.cluster.size - Score$Cluster.overlap)), "SSD" = ((Score$SSDGeneFamily - Score$LSDSSDGeneFamily)/Score$GeneFamilyNumber))
data$Module = ifelse(Score$GeneFamilySizeRoot > 0.5 & Score$Gene.cluster.size > 1900 & Score$Paralogous.cluster.size > 1900, "Similar", "Other")
data$Module = ifelse(Score$GeneFamilySizePollen > 0.5 & Score$Gene.cluster.size > 1900 & Score$Paralogous.cluster.size > 1900, "Similar", "Other")

# Pollen
pMain = ggplot(data, aes (x = Conserved, y = SSD, color = Module)) + geom_point() + geom_hline(yintercept = 0.7428571, linetype="dashed")
pRight = ggplot(data, aes(x = SSD, fill = Module, color = Module)) + geom_density(alpha = 0.3) + geom_vline(xintercept = 0.7428571, linetype="dashed") + coord_flip()

# Root
pMain = ggplot(data, aes (x = Conserved, y = SSD, color = Module)) + geom_point() + geom_hline(yintercept = 0.6666667, linetype="dashed")
pRight = ggplot(data, aes(x = SSD, fill = Module, color = Module)) + geom_density(alpha = 0.3) + geom_vline(xintercept = 0.6666667, linetype="dashed") + coord_flip()


pTop = ggplot(data, aes(x = Conserved, fill = Module, color = Module)) + geom_density(alpha = 0.3)
pEmpty <- ggplot(mtcars, aes(x = wt, y = mpg)) +
          geom_blank() +
          theme(axis.text = element_blank(),
                axis.title = element_blank(),
                line = element_blank(),
                panel.background = element_blank())


pdf("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/SSD/GeneFamilySSD_Enrichment_RootLeaf.pdf", width = 17, height = 10)
grid.arrange(pTop, pEmpty, pMain, pRight, ncol = 2, nrow = 2, widths = c(3, 1), heights = c(1, 3))
dev.off()





