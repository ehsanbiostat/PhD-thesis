library(foreach)
library(doMC)
registerDoMC(cores = 30)


load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/abiotic/Score99/G.score.RData")
load("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Local/Compendium_wise/abiotic/Score99/A.score.RData")

## Extract orthologous modules in the top 500.000 of grapevine
BB = unique(G.score[, "orthologous"])

## Extract paralogous modules
Result = foreach(i = 1:length(BB), .combine = rbind) %dopar% {
        A.score[A.score[, "Gene"] == BB[i] | A.score[, "Paralogous"] == BB[i],]
}

Result = Result[order(-Result[, "Paralogous.cluster"]),]
Result = Result[!duplicated(Result), ]
save(Result, file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Grape/Compendium_wise/abiotic/Score99/Result.RData")













