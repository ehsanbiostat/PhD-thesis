
library(preprocessCore)
library(foreach)
library(circlize)
source("/ngsprojects/iwt/ehsab/Gene_duplication/Scripts/General/iwanthue.R")

# Result from BINGO
PollenLeaf = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/GO/BINGO/PollenLeaf.txt")
Pollen = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/GO/BINGO/Pollen.txt")
RootLeaf = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/GO/BINGO/RootLeaf.txt")
Root = read.delim("/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/GO/BINGO/Root.txt")

# Remove GO terms with more than 335 and less than 20 genes
PollenLeaf = PollenLeaf[PollenLeaf$n < 335 & PollenLeaf$n > 20 & PollenLeaf$corr.p.value < 0.0001, ]
Pollen = Pollen[Pollen$n < 335 & Pollen$n > 20 & Pollen$corr.p.value < 0.001, ]
RootLeaf = RootLeaf[RootLeaf$n < 335 & RootLeaf$n > 20 & RootLeaf$corr.p.value < 0.0001, ]
Root = Root[Root$n < 335 & Root$n > 20 & Root$corr.p.value < 0.01, ]

All = unique(c(as.character(Pollen[, "Description"]), as.character(RootLeaf[, "Description"]), as.character(Root[, "Description"]), as.character(PollenLeaf[, "Description"])))

mat = matrix(0, length(All), 4)
colnames(mat) = c("Pollen", "PollenLeaf", "Root", "RootLeaf")
rownames(mat) = All
mat[as.character(Pollen$Description),"Pollen"] = Pollen$corr.p.value
mat[as.character(PollenLeaf$Description),"PollenLeaf"] = PollenLeaf$corr.p.value
mat[as.character(RootLeaf$Description),"RootLeaf"] = RootLeaf$corr.p.value
mat[as.character(Root$Description),"Root"] = Root$corr.p.value


matoriginal = mat
xx = log(mat)
xx = ifelse(xx == "-Inf", 0, xx)
mat = -xx


mat = t(normalize.quantiles(t(matoriginal)))
rownames(mat) = All



address = "/scratch/iwt/ehsab/Gene_duplication/Results/SA/Compare_SA__permut_thresold/development/GO/BINGO/Circular.pdf"

CirculHeatmap = function(GenExp, address) {
	# mat = as.matrix(GenExp)
	# if(Reorder) mat = as.matrix(ReorderMatrix(GenExp, condReverse, Row = F, Column = T)) # clustered matrixed on and column
	# mat = mat[c(seq(1,nrow(mat)/2),seq(nrow(mat),(nrow(mat)/2)+1)),]
	ColorRampPollen = colorRampPalette(c("white", "green"))(10)
	ColorRampPollenLeaf = colorRampPalette(c("white", "blue"))(10)
	ColorRampRootLeaf = colorRampPalette(c("white", "red"))(10)
	ColorRampRoot = colorRampPalette(c("white", "darkorange4"))(10)
	ColorRampAll = data.frame(ColorRampPollen, ColorRampPollenLeaf, ColorRampRootLeaf, ColorRampRoot)
	# ColorRamp = rgb(seq(0,1,length=256),seq(0,1,length=256),seq(1,0,length=256))
	# ColorLevels = seq(-max(mat), max(mat), length=length(ColorRamp))
	col_mat = foreach (i = 1:ncol(mat), .combine = cbind) %do% {
		ColorLevels = seq(0, max(mat[,i]), length = length(ColorRampPollen))
		col_fun = colorRamp2(ColorLevels, as.character(ColorRampAll[,i]))
		col_fun(mat[,i])
	}
	color = iwanthue(nrow(mat)/2)
	circos.clear()
	pdf(address, width = 10, height = 10)
	par(mar = c(1, 1, 1, 1))
	circos.par(cell.padding = c(0, 0, 0, 0), gap.degree = 1, start.degree = 90, clock.wise = TRUE)
	factors = letters[1]
	circos.par(points.overflow.warning = FALSE)
	circos.initialize(factors = factors, xlim = c(0, nrow(mat)))
	# Drawing the heatmap
	circos.trackPlotRegion(factors = NULL, ylim = c(0.05, ncol(mat) + 0.05), track.height = 0.7, bg.border = NA, panel.fun = function(x, y) {
		nr = nrow(mat)
		nc = ncol(mat)
		for(i in 1:nc) {
	    	for(j in 1:nr) {
	        	circos.rect(j-1, nc-i, j, nc-i+1, border = col_mat[j, i], col = col_mat[j, i])
	    	}
		}
	})
	circos.axis(# Writting gene names on the outside of the plot
		sector.index = "a", labels = rownames(mat),	major.at = seq(0.5, nrow(mat), by = 1),
		labels.cex = 0.25, major.tick = FALSE,  direction = "outside", labels.facing = "clockwise",
		labels.away.percentage = 0, major.tick.percentage = 0.01
	)
	# ylabels = sapply(1:length(rev(colnames(mat))), function(i) {strsplit(rev(colnames(mat))[i], "[*]")[[1]][1]})
	# for(i in 1:length(ylabels)) { #y axis labels, in this case it is the condition names
	# 	circos.text(0,i-0.5, ylabels[i], sector = "a", facing = "bending.inside", cex = 0.4, col = "black")
	# }
	dev.off()
}


CirculHeatmap(mat, address)
