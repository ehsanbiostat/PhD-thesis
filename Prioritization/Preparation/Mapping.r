

# Reference ID
ref=read.table("~/prioritization/Combined_Network/Ref_new.txt") # Reference ID 2 columns
colnames(ref)=c("ID","NO")



# Genie3 Total
# whole dataset 3 columns (148095 edges)
d=read.table("~/../../group/biocomp/users/ehsab/Genie3/development_weights.txt_0.02.txt")
# leaf development Network 3 columns (133232 edges)
t=read.table("~/../../group/biocomp/users/ehsab/Genie3/weights_inhouse_leaf2percent.txt")
# 281327 edges
d=rbind(d,t)
# No duplication (281327 edges)
d=d[!duplicated(d),]
# No loops (281327 edges)
d=d[!mapply(identical,as.vector(d[,1]),as.vector(d[,2])),]
colnames(d)[1]=c("ID")
# mapping the first column of d
k = merge(ref, d, by = "ID")
k = k[,-1]
colnames(k)[2]=c("ID")
# mapping the second column of d
k = merge(ref, k, by = "ID")
k = k[,-1]
write.table(k,file="~/../../group/biocomp/users/ehsab/Combined_Network/Indivitual_Networks/Genie3.txt",col.names=F,row.names=F,sep="\t")



# Text mining
# text mining 3 columns
t=read.table("~/../../group/biocomp/users/ehsab/Text_mining/Text_mining.txt")
colnames(t)[1]=c("ID")
# mapping the first column of t
k = merge(ref, t, by = "ID")
k = k[,-1]
colnames(k)[2]=c("ID")
# mapping the second column of t
k = merge(ref, k, by = "ID")
k = k[,-1]
write.table(k,file="~/../../group/biocomp/users/ehsab/Combined_Network/Indivitual_Networks/text.txt",col.names=F,row.names=F,sep="\t")


# Mutant
# Mutant 5 columns (142783 edges)
t=read.table("~/../../group/biocomp/users/ehsab/Mutants/mutants_single.txt")
# Removing duplicated edges (142434 edges)
t=t[!duplicated(t),]
# Remove loops (142299 edges)
t=t[!mapply(identical,as.vector(t[,1]),as.vector(t[,2])),]
# Number of genes (17145)
length(unique(c(as.vector(t[,1]),as.vector(t[,2]))))
# Removing 3 last columns
t=t[,c(1,2)]
colnames(t)[1]=c("ID")
# mapping the first column of t
k = merge(ref, t, by = "ID")
k = k[,-1]
colnames(k)[2]=c("ID")
# mapping the second column of t
k = merge(ref, k, by = "ID")
k = k[,-1]
write.table(k,file="~/../../group/biocomp/users/ehsab/Combined_Network/Indivitual_Networks/Mutant.txt",col.names=F,row.names=F,sep="\t")


# PCC
# PCC 3 columns (12526356 edges)
t=read.table("~/../../group/biocomp/users/ehsab/PCC/inhouse_leaf_pcc.txt")
# No duplicated edges (12526356 edges)
t=t[!duplicated(t),]
# No loops (12526356 edges)
t=t[!mapply(identical,as.vector(t[,1]),as.vector(t[,2])),]
# Number of genes (11763)
length(unique(c(as.vector(t[,1]),as.vector(t[,2]))))
colnames(t)[1]=c("ID")
# mapping the first column of t
k = merge(ref, t, by = "ID")
k = k[,-1]
colnames(k)[2]=c("ID")
# mapping the second column of t
k = merge(ref, k, by = "ID")
k = k[,-1]
write.table(k,file="~/../../group/biocomp/users/ehsab/Combined_Network/Indivitual_Networks/PCC.txt",col.names=F,row.names=F,sep="\t")


# PPI
# PPI 13 columns (1232284 edges)
t=read.table("~/../../group/biocomp/users/ehsab/Cornet_PPI/cornet_all_ppi_table_17012012.txt",sep="\t")
# select the three columns related to network
t=t[,c(2,3,4)]
# Remove duplicated edges (1149785 edges)
t=t[!duplicated(t),]
# Remove loops (1147589 edges)
t=t[!mapply(identical,as.vector(t[,1]),as.vector(t[,2])),]
# Number of genes (22372)
length(unique(c(as.vector(t[,1]),as.vector(t[,2]))))
colnames(t)[1]=c("ID")
# mapping the first column of t
k = merge(ref, t, by = "ID")
k = k[,-1]
colnames(k)[2]=c("ID")
# mapping the second column of t
k = merge(ref, k, by = "ID")
k = k[,-1]
write.table(k,file="~/../../group/biocomp/users/ehsab/Combined_Network/Indivitual_Networks/PPI.txt",col.names=F,row.names=F,sep="\t",quote=F)


# GeneMANIA
# Genemania 3 columns (10244303 edges)
t=read.table("~/../../group/biocomp/users/ehsab/Genemania_25_10_2012/COMBINED.DEFAULT_NETWORKS.BP_COMBINING.txt",header=T)
# Remove duplicated edges (10244303 edges)
t=t[!duplicated(t),]
# Remove loops (10244303 edges)
t=t[!mapply(identical,as.vector(t[,1]),as.vector(t[,2])),]
# Number of genes (24815)
length(unique(c(as.vector(t[,1]),as.vector(t[,2]))))
colnames(t)[1]=c("ID")
# mapping the first column of t
k = merge(ref, t, by = "ID")
k = k[,-1]
colnames(k)[2]=c("ID")
# mapping the second column of t
k = merge(ref, k, by = "ID")
k = k[,-1]
write.table(k,file="~/../../group/biocomp/users/ehsab/Combined_Network/Indivitual_Networks/Genemania.txt",col.names=F,row.names=F,sep="\t",quote=F)


# AGRIS
# AGRIS 2 columns (13042 edges)
t=read.table("~/../../group/biocomp/users/ehsab/AGRIS/REG_NET_2col.txt")
# Remove duplicated edges (13041 edges)
t=t[!duplicated(t),]
# Remove loops (13033 edges)
t=t[!mapply(identical,as.vector(t[,1]),as.vector(t[,2])),]
# Number of genes (9443)
length(unique(c(as.vector(t[,1]),as.vector(t[,2]))))
colnames(t)[1]=c("ID")
# mapping the first column of t
k = merge(ref, t, by = "ID")
k = k[,-1]
colnames(k)[2]=c("ID")
# mapping the second column of t
k = merge(ref, k, by = "ID")
k = k[,-1]
write.table(k,file="~/../../group/biocomp/users/ehsab/Combined_Network/Indivitual_Networks/AGRIS.txt",col.names=F,row.names=F,sep="\t",quote=F)



# IYG & REG genes
IYG=read.table("~/prioritization/IYG_atcodes_New1.txt")
REG=read.table("~/prioritization/nanostring_regulators.txt")
cregiyg=unique(rbind(IYG,REG))
colnames(cregiyg)[1]=c("ID")
colnames(IYG)[1]=c("ID")
colnames(REG)[1]=c("ID")
k = merge(ref, cregiyg, by = "ID")
write.table(k,file="~/prioritization/Combined_Network/New/Filtering_list/cregiyg.txt",col.names=F,row.names=F,sep="\t",quote=F)
k = merge(ref, IYG, by = "ID")
write.table(k,file="~/prioritization/Combined_Network/New/Filtering_list/iyg.txt",col.names=F,row.names=F,sep="\t",quote=F)
k = merge(ref, REG, by = "ID")
write.table(k,file="~/prioritization/Combined_Network/New/Filtering_list/reg.txt",col.names=F,row.names=F,sep="\t",quote=F)





# GO ID
g=read.table("~/prioritization/GO/ATH_GO_GOSLIM_modified.txt",sep="\t")
g=g[!duplicated(g),]
colnames(g)=c("ID","GO")
k = merge(ref, g, by = "ID")
d=data.frame(k[,1],k[,3],k[,2])
write.table(d,file="~/prioritization/Combined_Network/New/GO_ID.txt",col.names=F,row.names=F,sep="\t",quote=F)





























































