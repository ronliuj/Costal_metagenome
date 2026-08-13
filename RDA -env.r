library(ggplot2)
library(tidyr)
library(reshape2)
library(vegan)
library(dplyr)
library(ggrepel)
library(ggnewscale)
##adapt from tree.Rdata
meta<-read.table("meta_116-250124.txt",header = T,sep="\t")

selected<-c("TOC", "POC",  "MOC", "DOC", "Ca.OC", "Fe.OC","Lignin.phenol_real", "Microbial.necromass", 
            "pH" , "EC" , "TN", "C_N","C.P",  "P" ,"Fed",   "Feo", "Fep" , "Ald", "Alo" , "Alp","X.4", "X4.63" ,"X.63" ,
            "MAP" ,"MAT",
            "Lon", "Lat",
            "Density", "GDP","Agriculture")


carbon<-selected[1:8]


soil<-selected[9:23]

climate<-selected[24:25]
geograph<-selected[26:27]
human<-selected[28:30]



#is.element(selected,colnames(meta))

meta_sel<-meta[,selected]  
library(mice)
imp <- mice(meta_sel, method = 'rf', m = 5)

meta_filled <- complete(imp, 1)


abund<-bins[,3:118]

abundT<-data.frame(t(abund))

abundT<-abundT[row.names(meta),]

abundT.h<-decostand(abundT,method = "hellinger")

ene.test <- adonis2(abundT.h ~ Category, data = meta, permutations = 999, method="bray",na.rm=T)
# adonis r2=0.120 p= 0.001

ene.test2 <- adonis2(abundT.h ~ Area, data = meta, permutations = 999, method="bray",na.rm=T)
# adonis r2=0.117 p= 0.001

ene.test3 <- adonis2(abundT.h ~ Category+Area, data = meta, permutations = 999, method="bray",by="margin")
# adonis r2=0.03 p= 0.474 (area)





#NMDS
#D.abun3<- decostand(D.abun2, method = 'hellinger')
dist <- vegdist(abundT.h, method = 'bray',na.rm=T)
nmds<-metaMDS(dist,2)

nmds2<-metaMDS(dist,k=3)

stress <- nmds$stress
#0.14 应力小于0.2 可接受
stress2 <- nmds2$stress
#0.08
#将绘图数据转化为数据框
df <- as.data.frame(nmds2$points)
#与分组数据合并
df <- merge(meta[,c("Category","Area")], df,by="row.names")
row.names(df)<-df$Row.names
df<-df[,2:6]


eig<-nmds2$engine


##RDA or CCA


decorana(abundT.h)#all <4

#如果前4个轴的所有轴长度均小于3则选择RDA分析，如果4个轴中存在一个轴的长度大于4则选择CCA分析。

#差异膨胀因子分析


#差异膨胀因子分析,删除大于10的factor
 test.cca<-cca(abundT.h~.,meta_filled)
 vif.cca(test.cca)

##permutest检验模型显著性，enfit检验单个环境因子显著性
anova(test.cca)#p 0.001

rda.per<-permutest(test.cca,permu=999)

fit1<-envfit(nmds2,meta_filled,permu=999,na.rm=T)
fit_val1 <- scores(fit1, display = c("vectors"))
fit_val1 <- fit_val1*vegan::ordiArrowMul(fit_val1, fill = 2.5)

fit_env <- data.frame(cbind(fit_val1, fit1$vectors$r, fit1$vectors$pvals))
colnames(fit_env)[3:4]<-c("r2","pval")

fit_env<-mutate(fit_env,pd = cut(pval, breaks = c(-Inf, 0.01, 0.05, Inf),
                                 labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
fit_env$type<-c(rep("carbon",8),rep("soil",15),rep("climate",2),rep("geograph",2),rep("human",3))

ggplot() +
  geom_point(data=df, aes(x = MDS1, y = MDS2, color = Category,shape=Area),size = 2, alpha = 1) +
  scale_color_manual(values=c(MF = "#B22271", SM = "#2271B2", SG = "#71B222" )) +
  #new_scale("color") +
  geom_segment(data=data.frame(fit_env), arrow=arrow(length=unit(0.1,"cm"),type = "closed"),
               aes(x=0,y=0,xend=NMDS1, yend=NMDS2,linetype=pd), color="#C9CACA",alpha=1,linewidth=0.4)  + 
  #scale_linetype_manual(values=c("solid","dashed")) +
  new_scale("color")+
  geom_text_repel(data=data.frame(fit_env), aes(NMDS1, NMDS2, label=rownames(fit_env),color=type),
                   alpha=1,size = 4,
                   segment.color = 'grey35',
                   point.padding = unit(0.2,"lines"),max.overlaps = getOption("ggrepel.max.overlaps", default = 50))+
  scale_color_manual(values=c(carbon="#2EA7E0",soil="#C9A063",climate="#B5B5B6",geograph="#6C7A89",human="#6B5B95")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black"))+
  annotate("text",x=-2,y=-2,label="stress 0.11")

meta.v2<-read.table("meta_116-250411.txt",header = T,sep="\t")

meta.use<-meta.v2[,c("Category","Area",colnames(meta_filled))]

write.csv(meta.use,"Parameters.csv")

df$Area<-meta.v2[row.names(df),]$Area

df<-df %>% mutate(lab=case_when(Area=="BHS" ~ 'B',
                                Area=="SCS" ~ 'S',
                                Area=="ECS" ~ 'E',
                                Area=="HHS" ~ 'Y',
                                TRUE ~ 'Unkown'))

ggplot() +
  stat_ellipse(data=df, aes(x = MDS1, y = MDS2, color = Category),
               type = "norm", linetype = 2,level = 0.9)+
  stat_ellipse(data=df, aes(x = MDS1, y = MDS2, fill = Category),
               type = "norm", geom ="polygon",level = 0.85,alpha=0.3)+
  geom_text(data=df, aes(x = MDS1, y = MDS2, color = Category,label=lab),size = 3, alpha = 1,fontface="bold") +
  scale_color_manual(values=c(MF = "#B22271", SM = "#2271B2", SG = "#71B222" )) +
  scale_fill_manual(values=c(MF = "#B22271", SM = "#2271B2", SG = "#71B222" )) +
  #new_scale("color") +
  geom_segment(data=data.frame(fit_env), arrow=arrow(length=unit(0.1,"cm"),type = "closed"),
               aes(x=0,y=0,xend=NMDS1, yend=NMDS2,linetype=pd), color="#C9CACA",alpha=1,linewidth=0.4)  + 
  #scale_linetype_manual(values=c("solid","dashed")) +
  new_scale("color")+
  geom_text_repel(data=data.frame(fit_env), aes(NMDS1, NMDS2, label=rownames(fit_env),color=type),
                  alpha=1,size = 4,
                  segment.color = 'grey35',
                  point.padding = unit(0.2,"lines"),max.overlaps = getOption("ggrepel.max.overlaps", default = 50))+
  scale_color_manual(values=c(carbon="#2EA7E0",soil="#C9A063",climate="#B5B5B6",geograph="#6C7A89",human="#6B5B95")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black"))+
  annotate("text",x=-2,y=-2,label="stress 0.11")

#export 8x10

write.table(fit_env,"fit_env.txt",sep="\t",quote=F)

save.image("CCA.Rdata")

env2<-met.sel[,c("Cinnamyl","Syringyl","Vanillyl","Bacterial.necromass","Fungal.necromass")]

fit2<-envfit(nmds2,env2,permu=999,na.rm=T)
fit_val2 <- scores(fit2, display = c("vectors"))
fit_val2 <- fit_val2*vegan::ordiArrowMul(fit_val2, fill = 0.25)

fit_env2 <- data.frame(cbind(fit_val2, fit2$vectors$r, fit2$vectors$pvals))
colnames(fit_env2)[3:4]<-c("r2","pval")

fit_env2<-mutate(fit_env2,pd = cut(pval, breaks = c(-Inf, 0.01, 0.05, Inf),
                                   labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
fit_env2$type<-c(rep("plant",3),rep("micro",2))

ggplot() +
  geom_point(data=df2, aes(x = MDS1, y = MDS2, color = Category),size = 2, alpha = 1,shape=16) +
  scale_color_manual(values=c(MF = "#B22271", SM = "#2271B2", SG = "#71B222" )) +
  #new_scale("color") +
  geom_segment(data=data.frame(fit_env2), arrow=arrow(length=unit(0.1,"cm"),type = "closed"),
               aes(x=0,y=0,xend=NMDS1, yend=NMDS2,linetype=pd), color="#C9CACA",alpha=1,linewidth=0.4)  + 
  scale_linetype_manual(values=c("solid","dashed")) +
  new_scale("fill")+
  geom_label_repel(data=data.frame(fit_env2), aes(NMDS1, NMDS2, label=rownames(fit_env2),fill=type),
                   color='white',alpha=1,size = 4,
                   segment.color = 'grey35',
                   point.padding = unit(0.2,"lines"),max.overlaps = getOption("ggrepel.max.overlaps", default = 50))+
  scale_fill_manual(values=c(plant="#2EA7E0",micro="#C9A063")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black"))



env2<-met.sel[,c("Cinnamyl","Syringyl","Vanillyl","Bacterial.necromass","Fungal.necromass",
                 "POC","MOC","TOC","Ca.OC","Fe.OC","EC","Feo","Fep","Fed","X.4","X4.63","X.63","invsim.bin",
                 "shannon.bin","PCoA1.bin","shannon.cazy","invsim.cazy","PCoA1.cazy")]

fit2<-envfit(nmds2,env2,permu=999,na.rm=T)
fit_val2 <- scores(fit2, display = c("vectors"))
fit_val2 <- fit_val2*vegan::ordiArrowMul(fit_val2, fill = 0.25)

fit_env2 <- data.frame(cbind(fit_val2, fit2$vectors$r, fit2$vectors$pvals))
colnames(fit_env2)[3:4]<-c("r2","pval")

fit_env2<-mutate(fit_env2,pd = cut(pval, breaks = c(-Inf, 0.01, 0.05, Inf),
                                   labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
fit_env2$type<-c(rep("carbon",10),rep("soil",8),rep("mic",5))

ggplot() +
  geom_point(data=df2, aes(x = MDS1, y = MDS2, color = Category),size = 2, alpha = 1,shape=16) +
  scale_color_manual(values=c(MF = "#B22271", SM = "#2271B2", SG = "#71B222" )) +
  #new_scale("color") +
  geom_segment(data=data.frame(fit_env2), arrow=arrow(length=unit(0.1,"cm"),type = "closed"),
               aes(x=0,y=0,xend=NMDS1, yend=NMDS2,linetype=pd), color="#C9CACA",alpha=1,linewidth=0.4)  + 
  scale_linetype_manual(values=c("solid","dashed")) +
  new_scale("fill")+
  geom_label_repel(data=data.frame(fit_env2), aes(NMDS1, NMDS2, label=rownames(fit_env2),fill=type),
                   color='white',alpha=1,size = 4,
                   segment.color = 'grey35',
                   point.padding = unit(0.2,"lines"),max.overlaps = getOption("ggrepel.max.overlaps", default = 50))+
  scale_fill_manual(values=c(carbon="#2EA7E0",soil="#C9A063",mic="#B5B5B6")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black"))


bin<-read.table("../metagenome/bin.anno.txt",header=T,sep="\t")
row.names(bin)<-bin$Genomic.bins
bin<-bin[,2:74]
bin<-data.frame(t(bin))
bin<-bin[row.names(D.abun),]

bc.DOM<-vegdist(D.abun,"bray",na.rm=T) 
pcoa.DOM <- cmdscale(bc.DOM, k = (nrow(D.abun) - 1), eig = TRUE)
# 提取样本点坐标（points记录了各样本在各排序轴中的坐标值）
# 前5轴
DOM.pc <- data.frame(pcoa.DOM$point)[,1:5]
names(DOM.pc)[1:5] <- c('PCoA1', 'PCoA2', 'PCoA3', 'PCoA4', 'PCoA5')
# eig记录了PCoA排序结果中，主要排序轴的特征值（再除以特征值总和就是各轴的解释量）
eig.DOM = pcoa.DOM$eig


bc.bin<-vegdist(bin,"bray",na.rm=T) 
pcoa.bin <- cmdscale(bc.bin, k = (nrow(bin) - 1), eig = TRUE)
# 提取样本点坐标（points记录了各样本在各排序轴中的坐标值）
# 前5轴
bin.pc <- data.frame(pcoa.bin$point)[,1:5]
names(bin.pc)[1:5] <- c('PCoA1.bin', 'PCoA2.bin', 'PCoA3.bin', 'PCoA4.bin', 'PCoA5.bin')
# eig记录了PCoA排序结果中，主要排序轴的特征值（再除以特征值总和就是各轴的解释量）
eig.bin = pcoa.bin$eig

bin.pc<-bin.pc[row.names(DOM.pc),]


cazy<-read.table("../metagenome/cazy.sum.TPM",header=T,sep="\t")
cazy<-aggregate(.~Name,cazy,sum)
row.names(cazy)<-cazy$Name
cazy<-cazy[,2:74]
cazy<-data.frame(t(cazy))
cazy<-cazy[row.names(D.abun),]

bc.cazy<-vegdist(cazy,"bray",na.rm=T) 
pcoa.cazy <- cmdscale(bc.cazy, k = (nrow(cazy) - 1), eig = TRUE)
# 提取样本点坐标（points记录了各样本在各排序轴中的坐标值）
# 前5轴
cazy.pc <- data.frame(pcoa.cazy$point)[,1:5]
names(cazy.pc)[1:5] <- c('PCoA1.cazy', 'PCoA2.cazy', 'PCoA3.cazy', 'PCoA4.cazy', 'PCoA5.cazy')
# eig记录了PCoA排序结果中，主要排序轴的特征值（再除以特征值总和就是各轴的解释量）
eig.cazy = pcoa.cazy$eig


env2<-met.sel[,c("Cinnamyl","Syringyl","Vanillyl","Bacterial.necromass","Fungal.necromass")]
aa<-cbind(bin.pc,cazy.pc)
env2<-cbind(env2,aa[,c(2:3,7:8)])

fit2<-envfit(nmds2,env2,permu=999,na.rm=T)
fit_val2 <- scores(fit2, display = c("vectors"))
fit_val2 <- fit_val2*vegan::ordiArrowMul(fit_val2, fill = 0.8)

fit_env2 <- data.frame(cbind(fit_val2, fit2$vectors$r, fit2$vectors$pvals))
colnames(fit_env2)[3:4]<-c("r2","pval")

fit_env2<-mutate(fit_env2,pd = cut(pval, breaks = c(-Inf, 0.05, Inf),
                                   labels = c("< 0.05", ">= 0.05")))
fit_env2$type<-c(rep("carbon",5),rep("mic",4))

ggplot() +
  geom_point(data=df2, aes(x = MDS1, y = MDS2, color = Category),size = 2, alpha = 1,shape=16) +
  scale_color_manual(values=c(MF = "#B22271", SM = "#2271B2", SG = "#71B222" )) +
  #new_scale("color") +
  geom_segment(data=data.frame(fit_env2), arrow=arrow(length=unit(0.1,"cm"),type = "closed"),
               aes(x=0,y=0,xend=NMDS1, yend=NMDS2,linetype=pd), color="#C9CACA",alpha=1,linewidth=0.4)  + 
  scale_linetype_manual(values=c("solid","dashed")) +
  new_scale("fill")+
  geom_label_repel(data=data.frame(fit_env2), aes(NMDS1, NMDS2, label=rownames(fit_env2),fill=type),
                   color='white',alpha=1,size = 4,
                   segment.color = 'grey35',
                   point.padding = unit(0.2,"lines"),max.overlaps = getOption("ggrepel.max.overlaps", default = 50))+
  scale_fill_manual(values=c(carbon="#2EA7E0",mic="#B5B5B6")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black"))






Med<-aggregate(.~Category,DOM[,14:87],median,na.rm=TRUE, na.action=NULL)
ave<-aggregate(.~Category,DOM[,14:87],mean,na.rm=TRUE, na.action=NULL)
Sum<-aggregate(.~Category,DOM[,14:87],sum,na.rm=TRUE, na.action=NULL)

Sum<-Sum[c(1:6,8),]
for (i in 1:73){
  Sum[,i+1]<-Sum[,i+1]/sum(Sum[,i+1])
}


med<-melt(Med,id.vars = "Category",variable.name = "sample",value.name = "med")
ave<-melt(ave,id.vars = "Category",variable.name = "sample",value.name = "ave")
Sum<-melt(Sum,id.vars = "Category",variable.name = "sample",value.name = "sum")

sta<-merge(med,ave,by=c("Category","sample"))
sta<-merge(sta,Sum,by=c("Category","sample"))
sta$sample<-as.character(sta$sample)

para<-read.table("para-data2.txt",header = T,sep="\t")
row.names(para)<-para$ID
area<-read.table("Area.txt",header = T,sep="\t")
row.names(area)<-area$ID

sta$group<-para[sta$sample,]$Category
sta$area<-area[sta$sample,]$Area

Sum$sample<-as.character(Sum$sample)
Sum$group<-para[Sum$sample,]$Category


#div.roup<-aggregate(.~Category+group,sta[,c(1,3:6)],median,na.rm=TRUE, na.action=NULL)
#a.bin<-dcast(a.bin, Sample+Category+metcls ~ alpha, mean)

div.roup<-aggregate(.~Category+group,Sum[,c(1,3,4)],median,na.rm=TRUE, na.action=NULL)

med.group<-dcast(div.roup, Category~group )
# ave.group<-dcast(div.roup[,c(1,2,4)], Category~group )
# med.group<-dcast(div.roup[,c(1,2,3)], Category~group )

div.25<-aggregate(.~Category+group,Sum[,c(1,3,4)],function(x){quantile(x,0.25)}, na.action=NULL)
div.75<-aggregate(.~Category+group,Sum[,c(1,3,4)],function(x){quantile(x,0.75)}, na.action=NULL)
di25.group<-dcast(div.25, Category~group )
di75.group<-dcast(div.75, Category~group )


# div.group<-merge(med.group,ave.group,by="Category")
# div.group<-merge(div.group,med.group,by="Category")
# colnames(div.group)<-c("Category", "MF.sum","SM.sum","SG.sum","MF.ave","SM.ave",
#                        "SG.ave","MF.med","SM.med","SG.med")
# aa<-colSums(div.group[,2:10])

digoup2<-merge(med.group,di25.group,by="Category")
digoup2<-merge(digoup2,di75.group,by="Category")
colnames(digoup2)<-c("Category", "MF.med","SG.med","SM.med","MF.25","SG.25",
                       "SM.25","MF.75","SG.75","SM.75")
# digoup2<-digoup2[c(1:6,8),]
# for (i in 1:9){
#   digoup2[,i+1]<-digoup2[,i+1]/sum(digoup2[,i+1])
# }
# 
# 
# 
# sum.group2<-sum.group[c(1:6,8),]
# for (i in 1:3){
#   sum.group2[,i+1]<-sum.group2[,i+1]/sum(sum.group2[,i+1])
# }
write.table(digoup2,"mol-portion2.txt",row.names = F,sep="\t")
#por<-sum.group2
write.table(sum.group2,"mol-portion.txt",row.names = F,sep="\t")
por<-read.table("mol-portion.txt",header = T,sep="\t")

sum<-aggregate(.~Category,data=DOM[,c(14:87)],sum,na.rm=T,na.action = NULL)
sum<-sum[c(1:6,8),]
row.names(sum)<-sum$Category
sum<-sum[,-1]
aa<-apply(sum,1,median,na.rm=T)
aa<-data.frame(aa)
sum.por<-aa[,1]/sum(aa[,1])
sum.por<-data.frame(sum.por)
mol.por<-cbind(por,sum.por)
mol.por<-mol.por[,c(1,6,2:4)]

write.table(mol.por,"mol-portion.txt",row.names = F,sep="\t")

pie(mol.por$sum.por,border="white",col="orange",lwd=0.5,labels=mol.por$Category)

sum.group2$count<-round(sum.group2$MF*100)
#install.packages("ggpie")
library(ggpie)
?ggpie
ggpie(sum.group2, "count",group_key = "Category",fill_color = "grey")

pie(sum.group2$MF,border="white",col="#B22271",lwd=0.5,labels=sum.group2$Category)
pie(sum.group2$SM,border="white",col="#2271B2",lwd=0.5,labels=sum.group2$Category)
#pie(sum.group2$SG,labels=NA,border="white",col="#71B222",lwd=0.5)
pie(sum.group2$SG,border="white",col="#71B222",lwd=0.5,labels=sum.group2$Category)

ggpie3D(sum.group2, "count",group_key = "Category",tilt_degrees = -5)
?ggpie3D

row.names(meta)<-meta$ID
meta2<-meta[colnames(DOM)[15:87],]
plan.mean<-aggregate(Lignin.phenol~Category,data=meta2,mean,na.rm=T,na.action = NULL)
plan.median<-aggregate(Lignin.phenol~Category,data=meta2,median,na.rm=T,na.action = NULL)
plan.min<-aggregate(Lignin.phenol~Category,data=meta2,min,na.rm=T,na.action = NULL)
plan.max<-aggregate(Lignin.phenol~Category,data=meta2,max,na.rm=T,na.action = NULL)

mic.mean<-aggregate(Microbial.necromass~Category,data=meta2,mean,na.rm=T,na.action = NULL)
mic.median<-aggregate(Microbial.necromass~Category,data=meta2,median,na.rm=T,na.action = NULL)
mic.min<-aggregate(Microbial.necromass~Category,data=meta2,min,na.rm=T,na.action = NULL)
mic.max<-aggregate(Microbial.necromass~Category,data=meta2,max,na.rm=T,na.action = NULL)

aa<-bind_cols(plan.mean,plan.median,plan.min,plan.max,mic.mean,mic.median,mic.min,mic.max,.name_repair = "minimal")
aa<-aa[,c(1,2,4,6,8,10,12,14,16)]
colnames(aa)<-c("Category","plan.mean","plan.median","plan.min","plan.max","mic.mean","mic.median","mic.min","mic.max")
write.table(aa,"carbons.txt",row.names = F,sep="\t",quote=F)
