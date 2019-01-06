
load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/development/Score99/Result/Total1.RData")

# Select paralogous modules and their ancesteral modules from grapvine
# A pair of paralogous modules should have a same grapevine modules which has a significant enrichment score
A = sapply(g, unlist)
AA = sapply(A, length)

# Load paralogous modules from Arabidopsis data
load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/development/Score99/Result.RData")

# Load significant grape modules modules with arabidopsis modules
load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/development/Score99/G.score.RData")


## Grape
	library(igraph)
	condition.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/condition_name.txt", header = T)
	NAME = names(table(condition.names))

	k = 6
	
	## Grape modules
	load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Graph_V2.R")
	g.grape = g
	
	## Arabidopsis modules
	load(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/",NAME[k],"/Graph_V2.RData", sep = ""))
	
	orthologous = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Orthology/Grape-Arabidopsis/All.txt", header = T)
	gene.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Vitis_vinifera/gene.names.txt")
	map = read.delim("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Orthology/Grape-Arabidopsis/Map.12x.plaza.txt")
	Ref = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Reference_SD.txt", header = T)
	paralogous = read.table("~/Gene_duplication/Seed_genes/paralogous_total.txt", header = T)


	
library(foreach)
library(doMC)
registerDoMC(20)	
# INTER.Tot = c()

INTER.Tot1 = foreach (m = 1:length(which(AA > 1))) %dopar% {
	print(m)
	LEN = length(unlist(A[which(AA > 1)[m]]))
	INTER = ORTH1 = ORTH2 = c()
	for (n in 1:LEN)	{
		
		i = unlist(A[which(AA > 1)[m]])[n]

		memb.comun.grape = data.frame("name" = gene.names[neighborhood(g.grape, 1, i, mode = c("out"))[[1]],])
		memb.comun.grape.para = merge(memb.comun.grape, map, by = "name")
		## Remove duplicates 
		memb.comun.grape.para = memb.comun.grape.para[!(duplicated(memb.comun.grape.para$name) 
							| duplicated(memb.comun.grape.para$name, fromLast = TRUE)), ]
		memb.comun.grape.para = merge(memb.comun.grape.para, orthologous, by = "PLAZA.ID")
		
		j = Result[which(AA > 1)[m], "Gene"]
		memb.comun.Arab = data.frame("gene_id.y" = Ref[neighborhood(g, 1, j, mode = c("out"))[[1]],1])
		Merge1 = merge(memb.comun.Arab, memb.comun.grape.para, by = "gene_id.y")
		Orth1 = unique(Merge1$name)
		
		jpar = Result[which(AA > 1)[m], "Paralogous"]
		memb.comun.Arab = data.frame("gene_id.y" = Ref[neighborhood(g, 1, jpar, mode = c("out"))[[1]],1])
		Merge2 = merge(memb.comun.Arab, memb.comun.grape.para, by = "gene_id.y")
		Orth2 = unique(Merge2$name)

		Inter.orth = intersect(Orth2, Orth1)
		Inter.Para1 = Merge1[Orth1 %in% Inter.orth, "gene_id.y"]
		Inter.Para2 = Merge2[Orth2 %in% Inter.orth, "gene_id.y"]

		INTER[n] = length(intersect(Orth2, Orth1))
		ORTH1[n] = length(Orth1)
		ORTH2[n] = length(Orth2)
				
	}
	
	cbind(INTER, ORTH1, ORTH2)
}	
	
save(INTER.Tot1, file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/development/Score99/Result/Intersect_Orth.RData")
	
j = Result[which(Result[, "Paralogous"] == 1633 & Result[, "Gene"] == 8732), "Gene"]
TT.filter[TT.filter$orthologous == 1633 & TT.filter$orthologous.1 == 8732,]


# Extract shared paralogous and orthologous pairs in a triangile fashion
TT.para.orth = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Validation/Compendium_wise/development/5.Para.Orth.list.txt")

Inter = list()
for (n in 1:nrow(TT.para.orth))	{
	print(n)
	i = TT.para.orth$Gene.1[n]

	memb.comun.grape = data.frame("name" = gene.names[neighborhood(g.grape, 1, i)[[1]],])
	memb.comun.grape.para = merge(memb.comun.grape, map, by = "name")
	## Remove duplicates 
	memb.comun.grape.para = memb.comun.grape.para[!(duplicated(memb.comun.grape.para$name) 
						| duplicated(memb.comun.grape.para$name, fromLast = TRUE)), ]
	memb.comun.grape.para = merge(memb.comun.grape.para, orthologous, by = "PLAZA.ID")
	
	j = TT.para.orth$orthologous[n]
	memb.comun.Arab = data.frame("gene_id.y" = Ref[neighborhood(g, 1, j)[[1]],1])
	Merge1 = merge(memb.comun.Arab, memb.comun.grape.para, by = "gene_id.y")
	Orth1 = Merge1$name
	
	jpar = TT.para.orth$orthologous.1[n]
	memb.comun.Arab = data.frame("gene_id.y" = Ref[neighborhood(g, 1, jpar)[[1]],1])
	Merge2 = merge(memb.comun.Arab, memb.comun.grape.para, by = "gene_id.y")
	Orth2 = Merge2$name

	Inter.orth = intersect(Orth2, Orth1)
	Inter.Para1 = Merge1[Orth1 %in% Inter.orth, c("gene_id.y", "name")]
	Inter.Para2 = Merge2[Orth2 %in% Inter.orth, c("gene_id.y", "name")]
	Inter.Para1 = Inter.Para1[order(Inter.Para1$name),]
	Inter.Para2 = Inter.Para2[order(Inter.Para2$name),]

	Inter[n] = list(cbind.data.frame("gene_id.x" = Inter.Para1[,1], Inter.Para2))
	names(Inter)[n] = paste(i,j,jpar, sep = "_")
}
	
save(Inter, file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Validation/Compendium_wise/development/5.Para.Orth.Heatmap.R")
	






###############################################################################################################################################################################
# For two Pollen and Root modules
###############################################################################################################################################################################

Anchorpoint = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Anchorpoint/Geneset1_2_53_90.txt.txt")
Leaf = orthologous[as.character(orthologous$gene_id.y) %in% as.character(Anchorpoint$Leaf),]
Pollen = orthologous[as.character(orthologous$gene_id.y) %in% as.character(Anchorpoint$Pollen),]
MM = merge(Leaf, Pollen, by = "PLAZA.ID")
MM = MM[!duplicated(MM),]
colnames(MM)[2:3] = c("Leaf", "Pollen")

PollenLeaf = merge(map, MM, by = "PLAZA.ID")

INTER = foreach (j = 1:length(degree(g)), .combine = c) %dopar% {
	print(j)
	memb.comun.grape = data.frame("name" = gene.names[neighborhood(g.grape, 1, j, mode = c("out"))[[1]],])
	sum(as.character(memb.comun.grape$name) %in% as.character(unique(PollenLeaf$name)))
}	


Anchorpoint1 = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Anchorpoint/AT1G50090_AT3G19710_52_90.txt")
Leaf = orthologous[as.character(orthologous$gene_id.y) %in% as.character(Anchorpoint1$Leaf),]
Root = orthologous[as.character(orthologous$gene_id.y) %in% as.character(Anchorpoint1$Root),]
MM = merge(Leaf, Root, by = "PLAZA.ID")
MM = MM[!duplicated(MM),]
colnames(MM)[2:3] = c("Leaf", "Root")

RootLeaf = merge(map, MM, by = "PLAZA.ID")

INTERoot = foreach (j = 1:length(degree(g)), .combine = c) %dopar% {
	print(j)
	memb.comun.grape = data.frame("name" = gene.names[neighborhood(g.grape, 1, j, mode = c("out"))[[1]],])
	sum(as.character(memb.comun.grape$name) %in% as.character(unique(RootLeaf$name)))
}	

###############################################################################################################################################################################


INTERootRandomization = foreach (j = 1:10000, .combine = c) %dopar% {
	print(j)
	memb.comun.grape = data.frame("name" = gene.names[sample(1:length(degree(g)), 2103),])
	sum(as.character(memb.comun.grape$name) %in% as.character(unique(RootLeaf$name)))
}	

INTERandomization = foreach (j = 1:10000, .combine = c) %dopar% {
	print(j)
	memb.comun.grape = data.frame("name" = gene.names[sample(1:length(degree(g)), 2414),])
	sum(as.character(memb.comun.grape$name) %in% as.character(unique(PollenLeaf$name)))
}	



memb.comun.grape.para = merge(memb.comun.grape, map, by = "name")
	## Remove duplicates 
	memb.comun.grape.para = memb.comun.grape.para[!(duplicated(memb.comun.grape.para$name) 
						| duplicated(memb.comun.grape.para$name, fromLast = TRUE)), ]
	memb.comun.grape.para = merge(memb.comun.grape.para, orthologous, by = "PLAZA.ID")
	length(unique(memb.comun.grape.para[memb.comun.grape.para$PLAZA.ID %in% PollenLeaf$PLAZA.ID, "name"]))
	# cbind(INTER, ORTH1, ORTH2)




###############################################################################################################################################################################
# For two Pollen and Root modules - one by one and triangile relationships (one ancester and two daughters)
###############################################################################################################################################################################

library(foreach)
orthologous = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Orthology/Grape-Arabidopsis/All.txt", header = T)
map = read.delim("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Orthology/Grape-Arabidopsis/Map.12x.plaza.txt")
Ref = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Reference_SD.txt", header = T)
AnchPollen = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Anchorpoint/Geneset1_2_53_90.txt.txt")
AnchRoot = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/Genes/Anchorpoint/AT1G50090_AT3G19710_52_90.txt")


PollenLeaf = foreach (i = 1:nrow(AnchPollen), .combine = rbind) %do% { 
	Pollen = orthologous[orthologous$gene_id.y %in% AnchPollen$Pollen[i],]
	Leaf = orthologous[orthologous$gene_id.y %in% AnchPollen$Leaf[i],]
	if(nrow(Pollen) == 0) Pollen = data.frame("gene_id.y" = AnchPollen$Pollen[i], "PLAZA.ID" = NA)
	if(nrow(Leaf) == 0) Leaf = data.frame("gene_id.y" = AnchPollen$Leaf[i], "PLAZA.ID" = NA)
	Pollen["Module"] = "Pollen"
	Leaf["Module"] = "Leaf"
	cbind(Pollen, Leaf)
}


RootLeaf = foreach (i = 1:nrow(AnchRoot), .combine = rbind) %do% { 
	Root = orthologous[orthologous$gene_id.y %in% AnchRoot$Root[i],]
	Leaf = orthologous[orthologous$gene_id.y %in% AnchRoot$Leaf[i],]
	if(nrow(Root) == 0) Root = data.frame("gene_id.y" = AnchRoot$Root[i], "PLAZA.ID" = NA)
	if(nrow(Leaf) == 0) Leaf = data.frame("gene_id.y" = AnchRoot$Leaf[i], "PLAZA.ID" = NA)
	Root["Module"] = "Root"
	Leaf["Module"] = "Leaf"
	cbind(Root, Leaf)
}


PollenRootLeaf = rbind(PollenLeaf, RootLeaf)

PollenRootLeafclean = PollenRootLeaf[PollenRootLeaf[,2] == PollenRootLeaf[,5] & !is.na(PollenRootLeaf[,2]) & !is.na(PollenRootLeaf[,5]), -2]

PollenRootLeafclean["Name"] = NA
for(i in 1:nrow(PollenRootLeafclean)) {
	Name = map[map$PLAZA.ID %in% as.character(PollenRootLeafclean[, "PLAZA.ID"])[i],]
	if(nrow(Name) == 1) PollenRootLeafclean[i, "Name"] = as.character(map[map$PLAZA.ID %in% as.character(PollenRootLeafclean[, "PLAZA.ID"])[i],"name"])
}

PollenRootLeafclean = PollenRootLeafclean[!(is.na(PollenRootLeafclean$Name)),]
PollenRootLeafclean = PollenRootLeafclean[!(duplicated(PollenRootLeafclean)),]
write.table(PollenRootLeafclean, file = "/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Orthology/Grape-Arabidopsis/PollenRootLeafclean.txt", row.names = F, sep = "\t", quote = F)


###############################################################################################################################################################################
# For All APs
###############################################################################################################################################################################

library(foreach)
orthologous = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Orthology/Grape-Arabidopsis/All.txt", header = T)
map = read.delim("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Orthology/Grape-Arabidopsis/Map.12x.plaza.txt")
Ref = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/Reference_SD.txt", header = T)
rownames(Ref) = Ref[,1]
paralogous = read.table("~/Gene_duplication/Seed_genes/paralogous_NoFilter.txt", header = T)[1:2952, ]
paralogous["Ref"] = Ref[paralogous$Ref, 1]
paralogous["Anch"] = Ref[paralogous$Anch, 1]



PollenLeaf = foreach (i = 1:nrow(paralogous), .combine = rbind) %do% { 
	Pollen = orthologous[orthologous$gene_id.y %in% paralogous$Ref[i],]
	Leaf = orthologous[orthologous$gene_id.y %in% paralogous$Anch[i],]
	if(nrow(Pollen) == 0) Pollen = data.frame("gene_id.y" = paralogous$Ref[i], "PLAZA.ID" = NA)
	if(nrow(Leaf) == 0) Leaf = data.frame("gene_id.y" = paralogous$Anch[i], "PLAZA.ID" = NA)
	cbind(Pollen, Leaf)
}


PollenLeafclean = PollenLeaf[PollenLeaf[,2] == PollenLeaf[,4] & !is.na(PollenLeaf[,2]) & !is.na(PollenLeaf[,4]), -2]

PollenLeafclean["Name"] = NA
for(i in 1:nrow(PollenLeafclean)) {
	Name = map[map$PLAZA.ID %in% as.character(PollenLeafclean[, "PLAZA.ID"])[i],]
	if(nrow(Name) == 1) PollenLeafclean[i, "Name"] = as.character(map[map$PLAZA.ID %in% as.character(PollenLeafclean[, "PLAZA.ID"])[i],"name"])
}

PollenLeafclean = PollenLeafclean[!(is.na(PollenLeafclean$Name)),]
PollenLeafclean = PollenLeafclean[!(duplicated(PollenLeafclean)),]
write.table(PollenLeafclean, file = "/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Orthology/Grape-Arabidopsis/Triangle.txt", row.names = F, sep = "\t", quote = F)


