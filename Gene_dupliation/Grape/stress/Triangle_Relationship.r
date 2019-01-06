
load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/abiotic/Score99/Result/Total1.RData")

# Select paralogous modules and their ancesteral modules from grapvine
# A pair of paralogous modules should have a same grapevine modules which has a significant enrichment score
A = sapply(g, unlist)
AA = sapply(A, length)

# Load paralogous modules from Arabidopsis data
load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/abiotic/Score99/Result.RData")

# Load significant grape modules modules with arabidopsis modules
load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/abiotic/Score99/G.score.RData")


## Grape
	library(igraph)
	condition.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated CORNET2.0/condition_name.txt", header = T)
	NAME = names(table(condition.names))

	k = 1
	
	## Grape modules
	load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/stress/Graph_V2.R")
	g.grape = g
		
	## Arabidopsis modules
	load(paste("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Gene_score/Compendium_wise/",NAME[k],"/Graph_V2.RData", sep = ""))
	
	orthologous = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Orthology/Grape-Arabidopsis/All.txt", header = T)
	gene.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated CORNET2.0/Vitis_vinifera/Stress/gene.names.txt")
	Ref = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated CORNET2.0/Reference_SD.txt", header = T)
	paralogous = read.table("~/Gene_duplication/Seed_genes/paralogous_total.txt", header = T)


	
library(foreach)
library(doMC)
registerDoMC(16)	
# INTER.Tot = c()

INTER.Tot1 = foreach (m = 1:length(which(AA > 1))) %dopar% {
	print(m)
	LEN = length(unlist(A[which(AA > 1)[m]]))
	INTER = ORTH1 = ORTH2 = c()
	for (n in 1:LEN)	{
		
		i = unlist(A[which(AA > 1)[m]])[n]

		memb.comun.grape.para = data.frame("PLAZA.ID" = gene.names[neighborhood(g.grape, 1, i)[[1]],])
		## Remove duplicates 
		memb.comun.grape.para = data.frame("PLAZA.ID" = memb.comun.grape.para[!(duplicated(memb.comun.grape.para$PLAZA.ID) 
							| duplicated(memb.comun.grape.para$PLAZA.ID, fromLast = TRUE)), ])
		memb.comun.grape.para = merge(memb.comun.grape.para, orthologous, by = "PLAZA.ID")
		
		j = Result[which(AA > 1)[m], "Gene"]
		memb.comun.Arab = data.frame("gene_id.y" = Ref[neighborhood(g, 1, j)[[1]],1])
		Merge1 = merge(memb.comun.Arab, memb.comun.grape.para, by = "gene_id.y")
		Orth1 = Merge1$PLAZA.ID
		
		jpar = Result[which(AA > 1)[m], "Paralogous"]
		memb.comun.Arab = data.frame("gene_id.y" = Ref[neighborhood(g, 1, jpar)[[1]],1])
		Merge2 = merge(memb.comun.Arab, memb.comun.grape.para, by = "gene_id.y")
		Orth2 = Merge2$PLAZA.ID

		Inter.orth = intersect(Orth2, Orth1)
		Inter.Para1 = Merge1[Orth1 %in% Inter.orth, "gene_id.y"]
		Inter.Para2 = Merge2[Orth2 %in% Inter.orth, "gene_id.y"]

		INTER[n] = length(intersect(Orth2, Orth1))
		ORTH1[n] = length(Orth1)
		ORTH2[n] = length(Orth2)
				
	}
	
	cbind(INTER, ORTH1, ORTH2)
}	
	
save(INTER.Tot1, file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/abiotic/Score99/Result/Intersect_Orth.RData")
	

