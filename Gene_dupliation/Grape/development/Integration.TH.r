

library(foreach)
library(doMC)
registerDoMC(cores = 5)

#--------------------------------------------------------------------------------------------------------------------------------
## Read and prepare the file list for reading for Grape
list = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Scripts/General/order.TH.txt")
list = as.vector(unlist(list[1:19589,]))
setwd("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/development/Score99/")
g = lapply(list, read.table, sep= "\t" , header=T)

## Remove those file which has no significant module's score
ROW = lapply(g, nrow)
g = g[-which(unlist(ROW) == 0)]
ROW = ROW[-which(unlist(ROW) == 0)]
G.score = matrix(NA, sum(unlist(ROW)), ncol(data.frame(g[1])))
j = 0
for (i in 1:length(g)) {
	print(i)
	G.score[(j+1):(j+unlist(ROW[i])), ] = as.matrix(data.frame(g[i]))
	j = j+unlist(ROW[i])
}

colnames(G.score) = colnames(data.frame(g[i]))
## Order based on the Normalized.Paralogous.cluster score 
G.score = G.score[order(-G.score[, "orthologous.cluster"]), ]
## Select top 5.000.000 cases
G.score = G.score[1:5000000, ]
save(G.score, file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/development/Score99/G.score.RData")
# --------------------------------------------------------------------------------------------------------------------------------

job47666_I
job35892_I

# --------------------------------------------------------------------------------------------------------------------------------
## Read and prepare the file list for reading for Arabidopsis
list = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Scripts/General/order.TH.txt")
list = as.vector(unlist(list[1:19285,]))
setwd("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Local/Compendium_wise/development/Score99/")
g = lapply(list, read.table, sep= "\t" , header=T)

## Remove those file which has no significant module's score
ROW = lapply(g, nrow)
g = g[-which(unlist(ROW) == 0)]
ROW = ROW[-which(unlist(ROW) == 0)]
A.score = matrix(NA, sum(unlist(ROW)), ncol(data.frame(g[1])))
j = 0
for (i in 1:length(g)) {
	print(i)
	A.score[(j+1):(j+unlist(ROW[i])), ] = as.matrix(data.frame(g[i]))
	j = j+unlist(ROW[i])
}

colnames(A.score) = colnames(data.frame(g[i]))
save(A.score, file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Local/Compendium_wise/development/Score99/A.score.RData")
# --------------------------------------------------------------------------------------------------------------------------------




## Extract orthologous modules in the top 500.000 of grapevine
BB = unique(G.score[, "orthologous"])

## Extract paralogous modules
Result = foreach(1 = 1:length(BB), .combine = rbind) %do% {
	A.score[A.score[, "Gene"] == BB[i] | A.score[, "Paralogous"] == BB[i],]
}

## Order paralogous modules based on the number of shared paralogous
Result = Result[order(-Result[, "Paralogous.cluster"]),]

## Search orthologous modules for a pair of paralogous modules
Res = foreach (i = 1:nrow(Result)) %dopar% { 
	print(i)
	RR = table(G.score[(G.score[, "orthologous"] == Result[i , "Gene"] | 
		G.score[, "orthologous"] == Result[i , "Paralogous"]), "Gene"])
	names(RR[RR > 1])
	

}

## Remove those paralogous pairs without any common ancestral orthologous module
## and combine the results togheter
counter = which(sapply(Res, length) > 0)
Final = list()
for (i in 1:length(counter)) { 

	j = counter[i]
	RR = (G.score[(G.score[, "orthologous"] == Result[j , "Gene"] | 
		G.score[, "orthologous"] == Result[j , "Paralogous"]), ])
	Final[i] = list(RR[RR[, "Gene"] %in% unlist(Res[counter[i]]),])

}




































