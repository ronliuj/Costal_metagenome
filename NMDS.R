library(vegan)
bin<-read.table("bin_abundance_table.tab",header=T,sep="\t")
row.names(bin)<-bin$Genomic.bins
bin.abund<-bin[,2:117]

#bin.tr<- decostand(bin, method = 'hellinger')
dist <- vegdist(t(bin.abund), method = 'bray',na.rm=T)
nmds<-metaMDS(dist,2)

gene<-read.table("KO.sum.TPM",header=T,sep="\t")

#KO.tr<- decostand(gene[,2:117], method = 'hellinger')
dist2 <- vegdist(t(gene[,2:117]), method = 'bray',na.rm=T)
nmds2<-metaMDS(dist2,2)

para<-read.table("meta_251214.txt",header=T,sep="\t")
row.names(para)<-para$ID

#获得应力值（stress）
stress <- nmds$stress
#将绘图数据转化为数据框
df <- as.data.frame(nmds$points)
#与分组数据合并
df <- merge(para[,c(1,3)], df,by="row.names")

stress2 <- nmds2$stress
#将绘图数据转化为数据框
df2 <- as.data.frame(nmds2$points)
#与分组数据合并

df2 <- merge(para[,c(1,3)], df2,by="row.names")


ggplot(df, aes(MDS1, MDS2,color = Category))+
  geom_point(aes(color = Category), size = 2.5)+
  scale_color_manual(values=c(MF = "#B22271", SM = "#2271B2", SG = "#71B222" ))+
  theme_bw()+
  #xlim(c(-0.03,0.045))+
  stat_ellipse(level=0.6) +
  ggtitle(paste("stress: ",signif(stress,3)," Shared MAGs between two sites: Max 89.1% Min 87.8%",sep=""))+
  xlim(c(-3,5))+
  #geom_encircle(aes(fill=metcls), alpha = 0.1, show.legend = F) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black"))

ggplot(df2, aes(MDS1, MDS2,color = Category))+
  geom_point(aes(color = Category), size = 2.5)+
  scale_color_manual(values=c(MF = "#B22271", SM = "#2271B2", SG = "#71B222" ))+
  theme_bw()+
  #xlim(c(-0.03,0.045))+
  stat_ellipse(level=0.6) +
  ggtitle(paste("stress: ",stress2,sep=""))+
  #geom_encircle(aes(fill=metcls), alpha = 0.1, show.legend = F) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black"))


bin.abund<-t(bin.abund)
bin2 <- merge(para[,c(2:3)], bin.abund,by="row.names")
row.names(bin2)<-bin2$Row.names
bin.sum<-aggregate(.~Category, data=bin2[,c(3:1956)],sum)
row.names(bin.sum)<-bin.sum$Category
bin.sum<-bin.sum[,2:1954]
bin.sum<-data.frame(t(bin.sum))

MF<-row.names(bin.sum[bin.sum$MF!=0,])
SG<-row.names(bin.sum[bin.sum$SG!=0,])
SM<-row.names(bin.sum[bin.sum$SM!=0,])
intersect1<-length(intersect(MF,SG))/length(unique(c(MF,SG)))
intersect2<-length(intersect(MF,SM))/length(unique(c(MF,SM)))
intersect3<-length(intersect(SG,SM))/length(unique(c(SG,SM)))


#Dom.t<-data.frame(t(D.abun))
bin.t<-data.frame(t(bin2[,4:1956]))
bin.tf<-bin.t!=0

prop<-vector()
n=1
for (i in 1:115){
  aa<-names(bin.tf[,i][bin.tf[,i]==T])
  for(j in (i+1):116){
    bb<-names(bin.tf[,j][bin.tf[,j]==T])
    prop[n]<-length(intersect(aa,bb))/length(unique(c(aa,bb)))
    n=n+1
  }
}
min(prop)#0.14
max(prop)#0.72


bin.mean<-aggregate(.~Category, data=df[,c(3:5)],mean)

p1<-ggplot(df, aes(MDS1, MDS2,color = Category))+
  geom_point(aes(color = Category), size = 2.5)+
  geom_point(data=bin.mean,aes(MDS1, MDS2,color = Category), size = 2.5)+
  scale_color_manual(values=c(MF = "#B22271", SM = "#2271B2", SG = "#71B222" ))+
  theme_bw()+
  #xlim(c(-0.03,0.045))+
  stat_ellipse(level=0.6) +
  ggtitle(paste("MAGs stress: ",signif(stress,3),sep=""))+
  #xlim(c(-3,5))+
  geom_segment(data=df[df$Category=="MF",],aes(x=bin.mean[bin.mean$Category=="MF",2],y=bin.mean[bin.mean$Category=="MF",3],xend=MDS1,yend=MDS2),color="#B22271",arrow = arrow(length=unit(0.1,"inches")))+
  geom_segment(data=df[df$Category=="SG",],aes(x=bin.mean[bin.mean$Category=="SG",2],y=bin.mean[bin.mean$Category=="SG",3],xend=MDS1,yend=MDS2),color="#71B222",arrow = arrow(length=unit(0.1,"inches")))+
  geom_segment(data=df[df$Category=="SM",],aes(x=bin.mean[bin.mean$Category=="SM",2],y=bin.mean[bin.mean$Category=="SM",3],xend=MDS1,yend=MDS2),color="#2271B2",arrow = arrow(length=unit(0.1,"inches")))+
#geom_encircle(aes(fill=metcls), alpha = 0.1, show.legend = F) +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black"))

KO.mean<-aggregate(.~Category, data=df2[,c(3:5)],mean)

p2<-ggplot(df2, aes(MDS1, MDS2,color = Category))+
  geom_point(aes(color = Category), size = 2.5)+
  geom_point(data=KO.mean,aes(MDS1, MDS2,color = Category), size = 2.5)+
  scale_color_manual(values=c(MF = "#B22271", SM = "#2271B2", SG = "#71B222" ))+
  theme_bw()+
  #xlim(c(-0.03,0.045))+
  stat_ellipse(level=0.6) +
  ggtitle(paste("Genes stress: ",signif(nmds2$stress,3),sep=""))+
  #xlim(c(-3,5))+
  geom_segment(data=df2[df2$Category=="MF",],aes(x=KO.mean[KO.mean$Category=="MF",2],y=KO.mean[KO.mean$Category=="MF",3],xend=MDS1,yend=MDS2),color="#B22271",arrow = arrow(length=unit(0.1,"inches")))+
  geom_segment(data=df2[df2$Category=="SG",],aes(x=KO.mean[KO.mean$Category=="SG",2],y=KO.mean[KO.mean$Category=="SG",3],xend=MDS1,yend=MDS2),color="#71B222",arrow = arrow(length=unit(0.1,"inches")))+
  geom_segment(data=df2[df2$Category=="SM",],aes(x=KO.mean[KO.mean$Category=="SM",2],y=KO.mean[KO.mean$Category=="SM",3],xend=MDS1,yend=MDS2),color="#2271B2",arrow = arrow(length=unit(0.1,"inches")))+
  #geom_encircle(aes(fill=metcls), alpha = 0.1, show.legend = F) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black"))


ggpubr::ggarrange(p1,p2,nrow=2,ncol=1)

save.image("nmds.Rdata")

write.table(df,"nmds.mag.txt",sep="\t",quote=F,row.names = F)
write.table(df2,"nmds.gene.txt",sep="\t",quote=F,row.names = F)

bin2<-bin2[,4:1579]
para<-para[row.names(bin2),]
ene.test <- adonis2(bin2 ~ Met.cls, data = para, permutations = 999, method="bray",na.rm=T)
# adonis r2=0.037 p= 0.034
ene.test.mf <- adonis2(bin2[1:19,] ~ Met.cls, data = para[1:19,], permutations = 999, method="bray",na.rm=T)
# adonis r2=0.140 p= 0.092
ene.test.sm <- adonis2(bin2[20:56,] ~ Met.cls, data = para[20:56,], permutations = 999, method="bray",na.rm=T)
#adonis r2=0.067 p= 0.042
ene.test.sg <- adonis2(bin2[57:73,] ~ Met.cls, data = para[57:73,], permutations = 999, method="bray",na.rm=T)
#adonis r2=0.171 p= 0.071

bin.tr<- decostand(bin2, method = 'hellinger')
ene.test <- adonis2(bin.tr ~ Met.cls, data = para, permutations = 999, method="bray",na.rm=T)
# adonis r2=0.046 p= 0.019
ene.test.mf <- adonis2(bin.tr[1:19,] ~ Met.cls, data = para[1:19,], permutations = 999, method="bray",na.rm=T)
# adonis r2=0.156 p= 0.039
ene.test.sm <- adonis2(bin.tr[20:56,] ~ Met.cls, data = para[20:56,], permutations = 999, method="bray",na.rm=T)
#adonis r2=0.077 p= 0.048
ene.test.sg <- adonis2(bin.tr[57:73,] ~ Met.cls, data = para[57:73,], permutations = 999, method="bray",na.rm=T)
#adonis r2=0.203 p= 0.039