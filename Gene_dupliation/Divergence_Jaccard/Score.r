
library(foreach)
# install.packages("doMC", repos = getOption("http://cran.us.r-project.org"))
library(doMC)
registerDoMC(cores = 4)
# install.packages("doMPI", repos = getOption("http://cran.us.r-project.org"))
# library(doMPI)
# cl = startMPIcluster(16)
# registerDoMPI(cl)

load("/ngsprojects/iwt/ehsab/Gene_duplication/Results/SA/Local/Randomization/95/Threshold.RData")
# paralogous = read.table("~/Gene_duplication/Seed_genes/paralogous_total.txt", header = T)
paralogous = read.table("~/Gene_duplication/Seed_genes/paralogous_NoFilter.txt", header = T)
condition.names = read.table("/ngsprojects/iwt/ehsab/Gene_duplication/Integrated_CORNET2.0/condition_name.txt", header = T)
NAME = names(table(condition.names))

PARA = unique(paralogous$Ref)


k = 6

Result = foreach (i = PARA, .combine = rbind) %dopar%{
	print(i)
	a = read.table(paste("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Local/Compendium_wise/",NAME[k],"/Score99/Scores_V2/Score",i,".txt", sep = ""), header = T)
	# Select clusters with more than 3 nodes
	a[paralogous[paralogous$Ref %in% i, "Anch"], ]
}

write.table(Result, file = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Divergence_Jaccard/Score.txt", row.names = F, quote = F, sep = "\t")

















