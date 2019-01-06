
library(ggplot2)
library(csbl.go)
set.prob.table(organism=TAXONOMY.ARABIDOPSIS, type="similarity")
GO.Cons = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Divergence_Jaccard/GO.Modu.Cons.abiotic.large.txt")
GO.Div = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Divergence_Jaccard/GO.Modu.Div.abiotic.large.txt")
goids = rbind(data.frame("GO" = as.character(GO.Cons$GO), "Group" = "Cons"), data.frame("GO" = as.character(GO.Div$GO), "Group" = "Div"))
sim.matrix <- term.sim(goids$GO, "JiangConrathGraSM")
simil = data.frame("x" = sim.matrix[upper.tri(sim.matrix, diag = F)])

jpeg("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Divergence_Jaccard/GO/Plot/SemSim_All_abiotic.jpeg", width=8, height=5, units="in", res=800, pointsize = 5)
ggplot(simil, aes(x = x)) + geom_histogram(colour = "darkgreen", fill = "white", binwidth = 0.001) + xlab("GO semantic similarity")
dev.off()

sim.cons = sim.matrix[rownames(sim.matrix) %in% goids[goids$Group == "Cons", "GO"], colnames(sim.matrix) %in% goids[goids$Group == "Cons", "GO"]]
sim.div = sim.matrix[rownames(sim.matrix) %in% goids[goids$Group == "Div", "GO"], colnames(sim.matrix) %in% goids[goids$Group == "Div", "GO"]]
sim.div.cons = sim.matrix[rownames(sim.matrix) %in% goids[goids$Group == "Cons", "GO"], colnames(sim.matrix) %in% goids[goids$Group == "Div", "GO"]]
sim.all = rbind(data.frame("x" = sim.div[upper.tri(sim.div, diag = F)], "Group" = "Div"), 
	data.frame("x" = sim.cons[upper.tri(sim.cons, diag = F)], "Group" = "Cons"), 
	data.frame("x" = c(sim.div.cons), "Group" = "Cons.div"))

jpeg("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Divergence_Jaccard/GO/Plot/SemSim_sep_abiotic.jpeg", width=8, height=5, units="in", res=800, pointsize = 5)
ggplot(sim.all, aes(x = x, fill = Group)) + geom_density(alpha = 0.5) + xlab("GO semantic similarity")
dev.off()

# MDS Plot
library(MASS)
fit = isoMDS(as.dist(sim.matrix), k=2) # k is the number of dim

# plot solution 
x <- fit$points[,1]
y <- fit$points[,2]
xy = data.frame(x,y)


ggplot(xy, aes(x = x, y = y ,fill = Group)) + geom_density(alpha = 0.5) + xlab("GO semantic similarity")
