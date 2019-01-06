load("~/prioritization/Combined_Network/New/Total/Total_graph.RData")
iyg=read.table("~/prioritization/Combined_Network/New/Filtering_list/cregiyg.txt") # IYG and REG
iyg=as.vector(unlist(iyg[,2]))
w=list()
# Selecting GO term type: C=Cellular Components | F=Molecular Function | P=Biological Pathways
clas=c("C","F","P")
for (k in 1:3){
 go=read.table("~/prioritization/Combined_Network/New/GO_ID.txt",sep="\t")
 go=go[,-1]
 go=go[go[,3]==clas[k],]
 GOiyg=subset(go,go[,2]%in%iyg)[,1] # GO of IYG & REG
 freq=table(matrix(GOiyg)) # Frequency of GO among IYG & REG
 pvalue=odd=global=loc=c()

 for (i in 1:length(freq)){
  z=matrix(c(length(GOiyg)-freq[i],freq[i]),2)
  x=matrix(table(go[,1]==names(freq[i])))
  x=cbind(x,z) 
  pvalue[i]=fisher.test(x,alternative = "greater")$p
  global[i]=x[2,1]
  loc[i]=x[2,2]
  odd[i]=(x[1,1]*x[2,2])/(x[1,2]*x[2,1])
 }
 names(pvalue)=names(odd)=names(global)=names(loc)=names(freq)
 f=data.frame(pvalue,odd,global,loc)
 f=f[order(f[,1]),]
 top=c(rownames(f[1:10,]))

 library(igraph)
 e=pvalue_TopTen=odd_TopTen=Jaccard_N=Jaccard=Jaccard_Q=pvalue_first=odd_first=c()
 pvalue_5=odd_5=c()
 a1=Total_graph


 for (i in 1:27376){
	# Neighbors
	ne=unlist(neighborhood(a1,1,i))[-1]
	# GO terms of neighbors
	e=matrix(subset(go,go[,2]%in%ne)[,1])
	
	# GO enrichment analysis
	# 1. Top ten most significant terms within cregiyg
	z=matrix(table(e==top[1]|e==top[2]|e==top[3]|e==top[4]|e==top[5]|e==top[6]|e==top[7]|e==top[8]|e==top[9]|e==top[10]))
	if(dim(z)[1]==0) z=matrix(0,2,1) 
	if(dim(z)[1]==1) z=rbind(z,0)
	x=matrix(table(go[,1]==top[1]|go[,1]==top[2]|go[,1]==top[3]|go[,1]==top[4]|go[,1]==top[5]|go[,1]==top[6]|go[,1]==top[7]|go[,1]==top[8]|go[,1]==top[9]|go[,1]==top[10]))
	x=cbind(x,z)
	pvalue_TopTen[i]=fisher.test(x,alternative = "greater")$p
	odd_TopTen[i]=(x[1,1]*x[2,2])/(x[1,2]*x[2,1])
	
	# 2. The most significant term within cregiyg
	z=matrix(table(e==top[1]))
	if(dim(z)[1]==0) z=matrix(0,2,1) 
	if(dim(z)[1]==1) z=rbind(z,0)
	x=matrix(table(go[,1]==top[1]))
	x=cbind(x,z)
	pvalue_first[i]=fisher.test(x,alternative = "greater")$p
	odd_first[i]=(x[1,1]*x[2,2])/(x[1,2]*x[2,1])


	# 3. Top 5 most significant terms within cregiyg
        z=matrix(table(e==top[1]|e==top[2]|e==top[3]|e==top[4]|e==top[5]))
        if(dim(z)[1]==0) z=matrix(0,2,1)
        if(dim(z)[1]==1) z=rbind(z,0)
        x=matrix(table(go[,1]==top[1]|go[,1]==top[2]|go[,1]==top[3]|go[,1]==top[4]|go[,1]==top[5]))        
		x=cbind(x,z)
        pvalue_5[i]=fisher.test(x,alternative = "greater")$p
        odd_5[i]=(x[1,1]*x[2,2])/(x[1,2]*x[2,1])
	
	
	# Jaccard Similarity Coefficient
	# 1. The query gene
	q=matrix(subset(go,go[,2]%in%i)[,1])
	Jaccard_Q[i]=length(intersect(q,GOiyg))/length(union(q,GOiyg))
	
	# 2. Neighbors of query gene
	Jaccard_N[i]=length(intersect(e,GOiyg))/length(union(e,GOiyg))
	
	# 3. The query gene and its neighbors 
	t=rbind(q,e)
	Jaccard[i]=length(intersect(t,GOiyg))/length(union(t,GOiyg))
	
	print(i)
 }

 c=data.frame(pvalue_TopTen,odd_TopTen,pvalue_first,odd_first,pvalue_5,odd_5,Jaccard_Q,Jaccard_N,Jaccard)
 w[k]=list(c)
 #write.table(c, file="~/prioritization/Leave-one-out/New/GO/GO.txt",col.names=T,row.names=F,sep="\t")
}

save(w, file="~/prioritization/Leave-one-out/New/GO/GO.RData")













