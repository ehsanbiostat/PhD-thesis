

library(foreach)
library(doMC)
registerDoMC(6)

GeneFamily = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/GeneFamilyZhenFeb16.txt")
Freq = table(GeneFamily$V1)
GeneFamily = GeneFamily[GeneFamily$V1 %in% names(Freq[Freq > 1]),]

ref = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Reference_SD.txt", header = T)

GeneFamilyData = (merge(ref, GeneFamily, by.x = "Gene", by.y = "V2"))[, -2]
colnames(GeneFamilyData)[2] = "GeneFamily"
rownames(GeneFamilyData) = GeneFamilyData$Gene

Pollen = as.character(read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/Pollen_90.txt")[,1])
PollenLeaf = as.character(read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/Leaf_90.txt")[,1])
Root = as.character(read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/AT1G50090_Root_90.txt")[,1])
RootLeaf = as.character(read.table("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Module/AT3G19710_Leaf_90.txt")[,1])

AnchPollen = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Anchorpoint/Geneset1_2_53_90.txt.txt")
AnchRoot = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Anchorpoint/AT1G50090_AT3G19710_52_90.txt")

Pollen = setdiff(Pollen, as.character(AnchPollen$Pollen))
PollenLeaf = setdiff(PollenLeaf, as.character(AnchPollen$Leaf))
Root = setdiff(Root, as.character(AnchRoot$Root))
RootLeaf = setdiff(RootLeaf, as.character(AnchRoot$Leaf))

length(intersect(as.character(GeneFamilyData[Pollen, "GeneFamily"]), as.character(GeneFamilyData[PollenLeaf, "GeneFamily"])))
length(intersect(as.character(GeneFamilyData[Root, "GeneFamily"]), as.character(GeneFamilyData[RootLeaf, "GeneFamily"])))

GeneFamilyPermutePollen = PopulationPollen = rbind(GeneFamilyData[Pollen, ],GeneFamilyData[PollenLeaf, ])
GeneFamilyPermuteRoot = PopulationRoot = rbind(GeneFamilyData[Root, ],GeneFamilyData[RootLeaf, ])

GeneFamilyPermuteRoot = PopulationRoot = GeneFamilyPermutePollen = PopulationPollen = GeneFamilyData

Result = foreach (i = 1:1000, .combine = rbind) %dopar% {
	print(i)
	
	GeneFamilyPermutePollen$GeneFamily = sample(as.character(PopulationPollen$GeneFamily), nrow(PopulationPollen))
	GeneFamilyPermuteRoot$GeneFamily = sample(as.character(PopulationRoot$GeneFamily), nrow(PopulationRoot))
	
	data.frame("PollenLeaf" = length(intersect(as.character(GeneFamilyPermutePollen[Pollen, "GeneFamily"]), as.character(GeneFamilyPermutePollen[PollenLeaf, "GeneFamily"]))),
		"RootLeaf" = length(intersect(as.character(GeneFamilyPermuteRoot[Root, "GeneFamily"]), as.character(GeneFamilyPermuteRoot[RootLeaf, "GeneFamily"]))),
		"WithinPollen" = sum(table(as.character(GeneFamilyPermutePollen[Pollen, "GeneFamily"])) > 1), "WithinPollenLeaf" = sum(table(as.character(GeneFamilyPermutePollen[PollenLeaf, "GeneFamily"])) > 1),
		"WithinRoot" = sum(table(as.character(GeneFamilyPermuteRoot[Root, "GeneFamily"])) > 1), "WithinRootLeaf" = sum(table(as.character(GeneFamilyPermuteRoot[RootLeaf, "GeneFamily"])) > 1),
		"NumberPollen" = sum(table(as.character(GeneFamilyPermutePollen[Pollen, "GeneFamily"])) > 0), "NumberPollenLeaf" = sum(table(as.character(GeneFamilyPermutePollen[PollenLeaf, "GeneFamily"])) > 0),
		"NumberRoot" = sum(table(as.character(GeneFamilyPermuteRoot[Root, "GeneFamily"])) > 0), "NumberRootLeaf" = sum(table(as.character(GeneFamilyPermuteRoot[RootLeaf, "GeneFamily"])) > 0))
		
}

library(ggplot2)
ggplot(Result, aes(x = PollenLeaf)) + geom_histogram() + geom_hline(xintersept = 179)



library(gridexpand)
source("/ngsprojects/iwt/ehsab/Gene_duplication/Scripts/General/Multiple_plot.r")

P1 = ggplot(Result, aes(x = PollenLeaf)) + geom_histogram(colour = "darkgreen", fill = "white", binwidth = 1) + geom_vline(xintercept = 176, color = "red")
P2 = ggplot(Result, aes(x = RootLeaf)) + geom_histogram(colour = "darkgreen", fill = "white", binwidth = 1) + geom_vline(xintercept = 139, color = "red")
P3 = ggplot(Result, aes(x = WithinRoot)) + geom_histogram(colour = "darkgreen", fill = "white", binwidth = 1) + geom_vline(xintercept = 160, color = "red")
P4 = ggplot(Result, aes(x = WithinRootLeaf)) + geom_histogram(colour = "darkgreen", fill = "white", binwidth = 1) + geom_vline(xintercept = 202, color = "red")
P5 = ggplot(Result, aes(x = WithinPollen)) + geom_histogram(colour = "darkgreen", fill = "white", binwidth = 1) + geom_vline(xintercept = 230, color = "red")
P6 = ggplot(Result, aes(x = WithinPollenLeaf)) + geom_histogram(colour = "darkgreen", fill = "white", binwidth = 1) + geom_vline(xintercept = 198, color = "red")
pdf("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Plot/SSD/GeneFamily_PermutationTest_AllGeneFamilies.pdf", width = 20, height = 11)
multiplot(P1,P2,P3,P4,P5,P6, cols = 3)
dev.off()








PollenGeneFamily = sum(table(as.character(GeneFamilyData[Pollen, "GeneFamily"])) > 0)
RootGeneFamily = sum(table(as.character(GeneFamilyData[Root, "GeneFamily"])) > 0)

Population = c(names(table(as.character(GeneFamilyData[Pollen, "GeneFamily"]))), names(table(as.character(GeneFamilyData[PollenLeaf, "GeneFamily"]))))
Population = c(names(table(as.character(GeneFamilyData[Root, "GeneFamily"]))), names(table(as.character(GeneFamilyData[RootLeaf, "GeneFamily"]))))

Result = foreach (i = 1:10000, .combine = rbind) %dopar% {
	print(i)
	sample = sample(1:length(Population), RootGeneFamily, replace = F)
	sample1 = Population[setdiff(1:length(Population), sample)]
	sample2 = Population[sample]

	data.frame("Cross" = length(intersect(sample1, sample2)), "WhitinRootLeaf" = sum(table(sample1) > 1), "WhitinRoot" = sum(table(sample2) > 1),
		"RootLeaf" = sum(table(sample1) > 0), "Root" = sum(table(sample2) > 0))
}





Result = foreach (i = 1:10000, .combine = rbind) %dopar% {
	print(i)
	sample = sample(1:length(Population), RootGeneFamily, replace = F)
	sample1 = Population[setdiff(1:length(Population), sample)]
	sample2 = Population[sample]

	data.frame("Cross" = length(intersect(sample1, sample2)), "WhitinRootLeaf" = sum(table(sample1) > 1), "WhitinRoot" = sum(table(sample2) > 1),
		"RootLeaf" = sum(table(sample1) > 0), "Root" = sum(table(sample2) > 0))
}


