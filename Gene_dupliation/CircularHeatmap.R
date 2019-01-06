
pdf("/home/ehsab/Gene_duplication/Plot/test.pdf", width = 10, height = 10)
ggplot(nba.m, aes(x=Name, y=var2, fill=value)) +
     geom_tile() +
     scale_fill_gradient(low = "blue", high = "yellow") +
     ylim(c(0, max(nba.m$var2) + 0.5)) +
     scale_y_discrete(breaks=y_breaks, labels=y_labels) +
     coord_polar(theta="x") +
     theme(panel.background=element_blank(),
           axis.title=element_blank(),
           panel.grid=element_blank(),
           axis.text.x=element_blank(),
           axis.ticks=element_blank(),
           axis.text.y=element_text(size=5))
dev.off()


set.seed(999)
n = 1000
a = data.frame(factor = sample(letters[1:8], n, replace = TRUE), x = rnorm(n), y = runif(n))
library(circlize)

pdf("/home/ehsab/Gene_duplication/Plot/test1.pdf", width = 10, height = 10)
par(mar = c(1, 1, 1, 1), lwd = 0.1, cex = 0.7)
circos.par("track.height" = 0.1)
circos.initialize(factors = a$factor, x = a$y)
circos.trackPlotRegion(factors = a$factor, y = a$y, panel.fun = function(x, y) {
	circos.axis()
})
dev.off()

col = rep(c("#FF0000", "#00FF00"), 4)
circos.trackPoints(a$factor, a$x, a$y, col = col, pch = 16, cex = 0.5)
circos.text(-1, 0.5, "left", sector.index = "a", track.index = 1)
circos.text(1, 0.5, "right", sector.index = "a")
bgcol = rep(c("#EFEFEF", "#CCCCCC"), 4)
circos.trackHist(a$factor, a$x, bg.col = bgcol, col = NA)
dev.off()














op = par(no.readonly = TRUE)

library(circlize)
mat = matrix(rnorm(100*10), nrow = 10, ncol = 100)
col_fun = colorRamp2(c(-2, 0, 2), c("green", "black", "red"))
factors = rep(letters[1:2], 15) # which condition should be up and down on the plot
par(mar = c(1, 1, 1, 1))
circos.par(cell.padding = c(0, 0, 0, 0), gap.degree = 0	)
circos.initialize(factors, xlim = c(0, 15))
maxy = 0
circos.trackPlotRegion(ylim = c(0, 106), bg.border = NA, panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
    m = mat[, factors == sector.index]
    
    dend.col = as.dendrogram(hclust(dist(t(m))))

    maxy = ifelse(maxy > attr(dend.col, "height"), maxy, attr(dend.col, "height"))
    assign("maxy", maxy, envir = .GlobalEnv)

    m2 = m[, labels(dend.col)]
	col_mat = col_fun(m2)
    nr = nrow(m2)
    nc = ncol(m2)
    for(i in 1:nr) {
        for(j in 1:nc) {
            circos.rect(j-1, nr-i, j, nr-i+1, border = col_mat[i, j], col = col_mat[i, j])
        }
    }
    
})
circos.trackPlotRegion(ylim = c(0, maxy), bg.border = NA, track.height = 0.3, panel.fun = function(x, y) {
    sector.index = get.cell.meta.data("sector.index")
    m = mat[, factors == sector.index]
    
    dend.col = as.dendrogram(hclust(dist(t(m))))

    circos.dendrogram(dend.col, max_height = maxy)
    
})
circos.clear()










op = par(no.readonly = TRUE)
library(circlize)

mat = as.matrix(t(PollenLeaf))
ColorRamp = rgb(seq(0,1,length=256),seq(0,1,length=256),seq(1,0,length=256))
ColorLevels = seq(-max(mat), max(mat), length=length(ColorRamp))
col_fun = colorRamp2(ColorLevels, ColorRamp)
factors = rep(letters[1:2], ncol(mat)/2) # which condition should be up and down on the plot
factors = rep(letters[1], ncol(mat)) # which condition should be up and down on the plot

pdf("/home/ehsab/Gene_duplication/Plot/test1.pdf", width = 10, height = 10)
par(mar = c(1, 1, 1, 1))
circos.par(cell.padding = c(0, 0, 0, 0), gap.degree = 0	)
circos.initialize(factors, xlim = c(0, ncol(mat)))
maxy = 0
circos.trackPlotRegion(ylim = c(0, nrow(mat)), bg.border = NA, panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
    m = mat[, factors == sector.index]
    
    # dend.col = as.dendrogram(hclust(dist(t(m))))

    # maxy = ifelse(maxy > attr(dend.col, "height"), maxy, attr(dend.col, "height"))
    # assign("maxy", maxy, envir = .GlobalEnv)

    # m2 = m[, labels(dend.col)]
    m2 = m
	col_mat = col_fun(m2)
    nr = nrow(m2)
    nc = ncol(m2)
    for(i in 1:nr) {
        for(j in 1:nc) {
            circos.rect(j-1, nr-i, j, nr-i+1, border = col_mat[i, j], col = col_mat[i, j])
        }
    }
    
})
circos.link("a", 0, "b", 0, h = 0.4)
dev.off()



circos.trackPlotRegion(ylim = c(0, maxy), bg.border = NA, track.height = 0.3, panel.fun = function(x, y) {
    sector.index = get.cell.meta.data("sector.index")
    m = mat[, factors == sector.index]
    
    dend.col = as.dendrogram(hclust(dist(t(m))))

    circos.dendrogram(dend.col, max_height = maxy)
    
})
circos.clear()





ReorderMatrix = function(mydata) { 
	# set the custom distance and clustering functions, per your example
	hclustfunc <- function(x) hclust(x, method="ward.D2")
	distfunc <- function(x) dist(x, method="euclidean")

	# perform clustering on rows and columns
	cl.row <- hclustfunc(distfunc(t(mydata)))
	cl.col <- hclustfunc(distfunc(mydata))
	mydata = mydata[cl.col$order, rev(cl.row$order)]
	return(mydata)
}


AnchPollen$Pollen = sapply(1:nrow(AnchPollen), function(i) {ifelse(duplicated(AnchPollen$Pollen)[i], paste(AnchPollen$Pollen[i],1,sep = ""), as.character(AnchPollen$Pollen[i]))})
AnchPollen$Leaf = sapply(1:nrow(AnchPollen), function(i) {ifelse(duplicated(AnchPollen$Leaf)[i], paste(AnchPollen$Leaf[i],1,sep = ""), as.character(AnchPollen$Leaf[i]))})

AnchRoot$Root = sapply(1:nrow(AnchRoot), function(i) {ifelse(duplicated(AnchRoot$Root)[i], paste(AnchRoot$Root[i],1,sep = ""), as.character(AnchRoot$Root[i]))})
AnchRoot$Leaf = sapply(1:nrow(AnchRoot), function(i) {ifelse(duplicated(AnchRoot$Leaf)[i], paste(AnchRoot$Leaf[i],1,sep = ""), as.character(AnchRoot$Leaf[i]))})


library(circlize)

CirculHeatmap = function(GenExp, Reorder, condReverse, Module) {
	mat = as.matrix(GenExp)
	if(Reorder) mat = as.matrix(ReorderMatrix(GenExp, condReverse))
	ColorRamp = rgb(seq(0,1,length=256),seq(0,1,length=256),seq(1,0,length=256))
	ColorLevels = seq(-max(mat), max(mat), length=length(ColorRamp))
	col_fun = colorRamp2(ColorLevels, ColorRamp)
	col_mat = col_fun(mat)
	color = rainbow(nrow(mat)/2)

	circos.clear()
	# pdf(address, width = 10, height = 10)
	par(mar = c(1, 1, 1, 1))
	circos.par(cell.padding = c(0, 0, 0, 0), gap.degree = 3, start.degree = 90, clock.wise = TRUE)
	factors = letters[1]
	circos.par(points.overflow.warning = FALSE)
	circos.initialize(factors = factors, xlim = c(0, nrow(mat)))

	# Drawing the heatmap
	circos.trackPlotRegion(factors = NULL, ylim = c(0, ncol(mat)), track.height = 0.7, bg.border = NA, panel.fun = function(x, y) {
		nr = nrow(mat)
		nc = ncol(mat)
		for(i in 1:nc) {
	    	for(j in 1:nr) {
	        	circos.rect(j-1, nc-i, j, nc-i+1, border = col_mat[j, i], col = col_mat[j, i])
	    	}
		}
	})

	if(Module == "Pollen") {
		for (i in 1:c(nrow(mat)/2)) { # Anchorpoint connections in the middle of plot
			circos.link("a", which(rownames(mat) %in% as.character(AnchPollen$Pollen[i])) - 0.5,
				"a", which(rownames(mat) %in% as.character(AnchPollen$Leaf[i])) - 0.5, rou = 0.29, col = color[i])
		}
	}

	if(Module == "Root") {
		for (i in 1:c(nrow(mat)/2)) { # Anchorpoint connections in the middle of plot
			circos.link("a", which(rownames(mat) %in% as.character(AnchRoot$Root[i])) - 0.5,
				"a", which(rownames(mat) %in% as.character(AnchRoot$Leaf[i])) - 0.5, rou = 0.29, col = color[i])
		}
	}


	circos.axis(# Writting gene names on the outside of the plot
		sector.index = "a", labels = rownames(mat),	major.at = seq(0.5, nrow(mat), by = 1),
		labels.cex = 0.25, major.tick = FALSE,  direction = "outside", labels.facing = "outside",
		labels.away.percentage = 0, major.tick.percentage = 0.01
	)

	ylabels = sapply(1:length(rev(colnames(mat))), function(i) {strsplit(rev(colnames(mat))[i], "[*]")[[1]][1]})
	for(i in 1:length(ylabels)) { #y axis labels, in this case it is the condition names
		circos.text(0,i-0.5, ylabels[i], sector = "a", facing = "bending.inside", cex = 0.4, col = "black")
	}

	dev.off()
}



# circos.yaxis(side = c("right"), sector.index = "a", labels = 1:ncol(mat), at = seq(0.5, ncol(mat), by = 1), labels.cex = 0.4, tick = FALSE)
mat1 = mat[1:10, 1:10]
mat2 = mat[54:63, 1:10]

circos.clear()
# pdf("/home/ehsab/Gene_duplication/Plot/test2.pdf", width = 10, height = 10)
par(mar = c(1, 1, 1, 1))
circos.par(cell.padding = c(0, 0, 0, 0), gap.degree = 3)
factors = letters[1:2]
circos.par(points.overflow.warning = FALSE)
circos.initialize(factors = factors, xlim = c(0, nrow(mat1)))

col_mat = col_fun(mat1)
circos.trackPlotRegion(factors = factors, ylim = c(0, ncol(mat1)), track.height = 0.7, bg.border = NA, panel.fun = function(x, y) {
	col_mat = col_fun(mat1)
	nr = nrow(mat1)
	nc = ncol(mat1)
	for(i in 1:nc) {
    	for(j in 1:nr) {
        	circos.rect(j-1, nc-i, j, nc-i+1, border = col_mat[j, i], col = col_mat[j, i], sector.index = "a")
    	}
	}

	col_mat = col_fun(mat2)
	nr = nrow(mat2)
	nc = ncol(mat2)
	for(i in 1:nc) {
    	for(j in 1:nr) {
        	circos.rect(j-1, nc-i, j, nc-i+1, border = col_mat[j, i], col = col_mat[j, i], sector.index = "b")
    	}
	}


})

mat = mat[1:10, 1:10]
circos.clear()
# pdf("/home/ehsab/Gene_duplication/Plot/test2.pdf", width = 10, height = 10)
par(mar = c(1, 1, 1, 1))
circos.par(cell.padding = c(0, 0, 0, 0), gap.degree = 3)
factors = rep(letters[1:2], 5)
circos.par(points.overflow.warning = FALSE)
circos.initialize(factors = factors, xlim = c(0, nrow(mat1)/2))
circos.trackPlotRegion(ylim = c(0, ncol(mat)), track.height = 0.7, bg.border = NA, panel.fun = function(x, y) {
	
	sector.index = get.cell.meta.data("sector.index")
    m1 = mat[factors == sector.index,]
    m = m1[1:5,]
    

    col_mat = col_fun(m)
	nr = nrow(m)
	nc = ncol(m)
	for(i in 1:nc) {
    	for(j in 1:nr) {
        	circos.rect(j-1, nc-i, j, nc-i+1, border = col_mat[j, i], col = col_mat[j, i])
    	}
	}
})




# color = rgb(seq(0,1,length=256),seq(0,1,length=256),seq(1,0,length=256))
color = rainbow(nrow(mat)/2)

for (i in 1:c(nrow(mat)/2)) {
	circos.link("a", which(rownames(mat) %in% as.character(AnchPollen$Pollen[i])) - 0.5,
		"a", which(rownames(mat) %in% as.character(AnchPollen$Leaf[i])) - 0.5, rou = 0.25, col = color[i])
}
 circos.axis(sector.index = "a", labels = rownames(mat),
 	major.at = seq(0.5, nrow(mat), by = 1), labels.cex = 0.25, major.tick = FALSE,  direction = "outside", labels.facing = "outside", labels.away.percentage = 0,
 	major.tick.percentage = 0.01)

circos.yaxis(side = c("right"), sector.index = "a", labels = colnames(mat), at = seq(0.5, ncol(mat), by = 1), labels.cex = 0.25, tick = FALSE)
dev.off()



































# k = 0
# for(i in 1:26) {
# 	k = k + 0.015
# 	circos.link("a", i-0.5, "a", 106-i+0.5, rou = 0.5, h = k)
# }
# circos.link("a", 25.5, "a", 80.5, rou1 = 0.5, rou2 = 0.5, col = "black", border = "black", h = 0.49)
# circos.link("a", 26.5, "a", 79.5, rou1 = 0.5, rou2 = 0.5, col = "black", border = "black", h = 0.52)
# circos.link("a", 27.5, "a", 78.5, rou1 = 0.5, rou2 = 0.5, col = "black", border = "black", h = 0.54)
# circos.link("a", 28.5, "a", 77.5, rou1 = 0.5, rou2 = 0.5, col = "black", border = "black", h = 0.55)
# k = 0.55
# for(i in 31:53) {
# 	k = k + 0.01
# 	circos.link("a", i-0.5, "a", 106-i+0.5, rou = 0.5, h = k)
# }


labels.facing = "reverse.clockwise",

circos.trackPlotRegion(factors = factors, ylim = c(0, 30))
circos.text(5, 10, get.cell.meta.data("sector.index"))

circos.axis(sector.index = "a")
circos.rect(1,1,2,3, border = "red", col = "black")
circos.rect(1,3,2,5, border = "red", col = "blue")
circos.rect(1,5,1.1,7, border = "red", col = "green")
circos.axis(sector.index = "b", direction = "inside", labels.facing = "outside")
circos.axis(sector.index = "c", h = "bottom")
circos.axis(sector.index = "d", h = "bottom", direction = "inside", labels.facing = "reverse.clockwise")
circos.axis(sector.index = "e", h = 5, major.at = c(1, 3, 5, 7, 9))
dev.off()
