combine all AT codes from different datasources



# text mining data set
# 3 colums (8803 edges)
# 1- Remove all AT codes related to DNA ribosom such as "ArthMp008" order and manually remove them (8717 edges) manually was done
t=read.table("~/../../group/biocomp/users/ehsab/Text_mining/Text_mining.txt")
# 2- Remove duplicated connection (7756 edges)
t=t[!duplicated(t),]
# 3- Remove loops (6499 edges)
t=t[!mapply(identical,as.vector(t[,1]),as.vector(t[,2])),]
# Number of genes (2335)
length(unique(c(as.vector(t[,1]),as.vector(t[,2]))))
write.table(t, file="~/../../group/biocomp/users/ehsab/Text_mining/Text_mining.txt",sep="\t",col.names=F,row.names=F, quote=F)




# Genie3 data set
# 3 columns(148095 edges)
g=read.table("~/../../group/biocomp/users/ehsab/Genie3/development_weights.txt_0.02.txt")
g1=read.table("~/../../group/biocomp/users/ehsab/Genie3/weights_inhouse_leaf2percent.txt")
g=rbind(g,g1)
# There is no duplication (148095 edges)
g=g[!duplicated(g),]
# There is no loop (148095 edges)
g=g[!mapply(identical,as.vector(g[,1]),as.vector(g[,2])),]
# Number of genes (21414)
length(unique(c(as.vector(g[,1]),as.vector(g[,2]))))



# GeneMANIA dataset (version 25-10-2012)
# 3 columns (10244303) edges
G=read.table("~/../../group/biocomp/users/ehsab/Genemania_25_10_2012/COMBINED.DEFAULT_NETWORKS.BP_COMBINING.txt",header=T)
# No duplicated connection (10244303 edges)
G=G[!duplicated(G),]
# No loops (10244303 edges)
G=G[!mapply(identical,as.vector(G[,1]),as.vector(G[,2])),]
# Number of genes (24815)
length(unique(c(as.vector(G[,1]),as.vector(G[,2]))))




# Mutant
# Mutant 5 columns (142783 edges)
m=read.table("~/../../group/biocomp/users/ehsab/Mutants/mutants_single.txt")
# Removing duplicated edges (142434 edges)
m=m[!duplicated(m),]
# Remove loops (142299 edges)
m=m[!mapply(identical,as.vector(m[,1]),as.vector(m[,2])),]


# PCC
# PCC 3 columns (12526356 edges)
p=read.table("~/../../group/biocomp/users/ehsab/PCC/inhouse_leaf_pcc.txt")
# No duplicated edges (12526356 edges)
p=p[!duplicated(p),]
# No loops (12526356 edges)
p=p[!mapply(identical,as.vector(p[,1]),as.vector(p[,2])),]


# PPI 13 columns (1232284 edges)
pp=read.table("~/../../group/biocomp/users/ehsab/Cornet_PPI/cornet_all_ppi_table_17012012.txt",sep="\t")
# select the three columns related to network
pp=pp[,c(2,3,4)]
# Remove duplicated edges (1149785 edges)
pp=pp[!duplicated(pp),]
# Remove loops (1147589 edges)
pp=pp[!mapply(identical,as.vector(pp[,1]),as.vector(pp[,2])),]




# AGRIS
# AGRIS 2 columns (13042 edges)
a=read.table("~/../../group/biocomp/users/ehsab/AGRIS/REG_NET_2col.txt")
# Remove duplicated edges (13041 edges)
a=a[!duplicated(a),]
# Remove loops (13033 edges)
a=a[!mapply(identical,as.vector(a[,1]),as.vector(a[,2])),]






# Reference file
# Number of genes (27376 genes)
length(unique(c(as.vector(g[,1]),as.vector(g[,2]),as.vector(a[,1]),as.vector(a[,2]),as.vector(t[,1]),as.vector(t[,2]),as.vector(G[,1]),as.vector(G[,2])
,as.vector(m[,1]),as.vector(m[,2]),as.vector(p[,1]),as.vector(p[,1]),as.vector(pp[,1]),as.vector(pp[,2]))))

new_ref=unique(c(as.vector(g[,1]),as.vector(g[,2]),as.vector(a[,1]),as.vector(a[,2]),as.vector(t[,1]),as.vector(t[,2]),as.vector(G[,1]),as.vector(G[,2])
,as.vector(m[,1]),as.vector(m[,2]),as.vector(p[,1]),as.vector(p[,1]),as.vector(pp[,1]),as.vector(pp[,2])))

new_ref=as.matrix(new_ref)
new_ref=new_ref[order(new_ref[,1]),]
new_ref=cbind(new_ref,c(1:27376))
write.table(new_ref, file="~/prioritization/Combined_Network/Ref_new.txt",sep="\t",col.names=F,row.names=F, quote=F)



