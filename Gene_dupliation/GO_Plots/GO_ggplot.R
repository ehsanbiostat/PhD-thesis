library("reshape")
library("ggplot2")

GOterms <- as.matrix(read.table("//psb.ugent.be/shares/biocomp/groups/group_biocomp/users/jofoo/ERC/Zhenset_080615/GO/GOslimSummaryRef.txt", header=T,row.names = 1,sep = "\t" ))
colnames(GOterms)<-c("Group 1","Group 2","Group 3")
GO<-melt(GOterms)
colnames(GO)<-c("GOcat","Group","Pvalue")

GO$GOcat<- factor(GO$GOcat,levels=rev(unique(GO$GOcat)))

underover<-vector(length = length(GO$Pvalue))
underover[which(GO$Pvalue<0)]<-"Under"
underover[which(GO$Pvalue>=0)]<-"Over"
underover[which(is.na(GO$Pvalue))]<-"No"
GO<-cbind(GO,underover)
colnames(GO)<-c("GOcat","Group","Pvalue","Direction")
GO$Pvalue[which(GO$Pvalue<0)]<- GO$Pvalue[which(GO$Pvalue<0)]*-1
GO$Pvalue[which(GO$Pvalue<1e-10)]<-(1e-10)
GO$Pvalue<-(-log(GO$Pvalue))

ggplot(GO,aes(y=GOcat,x=Group))+
  geom_point(aes(colour = Direction,size = Pvalue))+
  scale_size_continuous(breaks= c(2,6.6,10,16.7,20,25),
                        labels= c("0","< 0.05","< 0.01","< 1e-3","< 1e-5","<0.0001"), range = c(2, 4))+
  scale_colour_manual(breaks = c("Over","Under","No"),
                      values = c("white","darkgreen","red"),
                      labels=c("Over","Under","No"))+
  theme_bw()+
  theme(legend.text = element_text(size = 9, face = "bold"),
        axis.text.y = element_text(size = 7,color="black"),#,angle = 45,),
        axis.text.x = element_text(size = 9,color="black",angle = 90,vjust = 0.5),
        axis.title.y = element_blank(),
        axis.title.x = element_blank() )

