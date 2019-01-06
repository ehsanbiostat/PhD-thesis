

library(foreach)
library(doMC)
registerDoMC(30)


## GO terms loading
go = read.table("/ngsprojects/iwt/ehsab/Prioritization/prioritization/Combined_Network/New/GO_ID.txt",sep = "\t")
clas = c("C", "F", "P")
k = 3
go = go[go[, 4] == clas[k], ]

GO.terms = unique(go[,2])
# GO.dist = matrix(NA, length(GO.terms), length(GO.terms))
GO.dist = c()
length(GO.terms)
GO.Dist = foreach(i = 1:10, .combine = rbind) %dopar% {
	print(i)
	GO.inter = go[go[,2] == GO.terms[i],1]
	for(j in 1:10) {
		GO.dist[j] = length(intersect(GO.inter, go[go[,2] == GO.terms[j],1]))/length(union(GO.inter, go[go[,2] == GO.terms[j],1]))
	}
	GO.dist
}

colnames(GO.Dist) = rownames(GO.Dist) = GO.terms
save(GO.Dist, file = "/ngsprojects/iwt/ehsab/Prioritization/prioritization/Combined_Network/New/GODistanceGeneBase.RData")
# write.table(GO.Dist, file = "", col.names = F, row.names = F, sep = "\t", qoute = F)