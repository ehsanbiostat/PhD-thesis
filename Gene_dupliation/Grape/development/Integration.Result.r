
library(foreach)
library(doMC)
registerDoMC(20)


load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/development/Score99/G.score.RData")
load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Local/Compendium_wise/development/Score99/A.score.RData")
load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/development/Score99/Result.RData")


g = foreach (i = 1:50000) %dopar% {
	print(i)
	RR = table(G.score[(G.score[, "orthologous"] == Result[i , "Gene"] |
	                G.score[, "orthologous"] == Result[i , "Paralogous"]), "Gene"])
	as.numeric(names(RR[RR > 1]))
}

save(g, file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/development/Score99/Result/Total1.RData")

