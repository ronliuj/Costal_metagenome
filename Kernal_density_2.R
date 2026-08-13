library(ggplot2) 
library(gcookbook)
library(RColorBrewer)
library(scales)
library(geosphere)
library(reshape2)
library(vegan)
library(dplyr)
library(viridis) # 提供美观的颜色方案
library(gridExtra) # 用于组合多个图形

abund<-read.table("bin_abundance_table.tab",header = T,sep="\t")
row.names(abund)<-abund$Genomic.bins
KO<-read.table("KO.sum.TPM",header = T,sep="\t")
row.names(KO)<-KO$Name
meta<-read.table("meta_250715.txt",header = T,sep="\t")
row.names(meta)<-meta$ID
meta$ind<-seq(1:116)

bc.MAG<-vegan::vegdist(t(abund[,row.names(meta)]),method = "bray")
bc.MAG<-as.matrix(bc.MAG)

bc.fun<-vegan::vegdist(t(KO[,row.names(meta)]),method = "bray")
bc.fun<-as.matrix(bc.fun)


bcML <- melt(bc.MAG, 
            varnames = c("Sample1", "Sample2"), value.name = "Distance")

bcML$Ind1<-as.numeric(bcML$Sample1)
bcML$Ind2<-as.numeric(bcML$Sample2)

bcML <- bcML[bcML$Ind1 < bcML$Ind2, ]


bcML$Sample1<-as.character(bcML$Sample1)
bcML$Sample2<-as.character(bcML$Sample2)

bcML$Ind1<-meta[bcML$Sample1,]$ind
bcML$Ind2<-meta[bcML$Sample2,]$ind

bcML<-bcML%>%arrange(Ind1,Ind2)
colnames(bcML)[3]<-"distMAG"

bcFL <- melt(bc.fun, 
             varnames = c("Sample1", "Sample2"), value.name = "Distance")

bcFL$Ind1<-as.numeric(bcFL$Sample1)
bcFL$Ind2<-as.numeric(bcF$Sample2)

bcFL <- bcFL[bcFL$Ind1 < bcFL$Ind2, ]


bcFL$Sample1<-as.character(bcFL$Sample1)
bcFL$Sample2<-as.character(bcFL$Sample2)

bcFL$Ind1<-meta[bcFL$Sample1,]$ind
bcFL$Ind2<-meta[bcFL$Sample2,]$ind

bcFL<-bcFL%>%arrange(Ind1,Ind2)
colnames(bcFL)[3]<-"distFun"

bc<-left_join(bcFL,bcML,by=c("Sample1", "Sample2","Ind1","Ind2"))

write.table(bc,"bc.txt",sep="\t",quote=F)

bc.a<-ggplot(data = bc, aes(y=distFun, x=distMAG)) +
  geom_hex(bins = 60, # 调整分箱数量以控制粒度
  #binwidth = c(0.01, 0.005), # 替代bins参数，直接设置宽度
  aes(fill = after_stat(density)), # 使用密度值着色
  alpha = 0.9 )+    # 设置透明度
  scale_fill_viridis(option = "turbo", # 使用viridis配色方案
    name = "Density", 
    trans = "log", # 对数变换使低密度区域更明显
    breaks = scales::log_breaks())+
  theme_bw()+
  theme(panel.grid.major = element_line(color="lightgrey"), panel.grid.minor = element_line(color="lightgrey"))
  
  
MAG.presence<-ifelse(abund[,row.names(meta)] > 0, 1, 0)
MAG.sor <- vegan::betadiver(t(MAG.presence), method = "sor")  
MAG.sor<-as.matrix(MAG.sor)  



Fun.presence<-ifelse(KO[,row.names(meta)] > 0, 1, 0)
Fun.sor <- vegan::betadiver(t(Fun.presence), method = "sor")  
Fun.sor<-as.matrix(Fun.sor)  

SorML <- melt(MAG.sor, 
            varnames = c("Sample1", "Sample2"), value.name = "Distance")  
  
SorML$Ind1<-as.numeric(SorML$Sample1)
SorML$Ind2<-as.numeric(SorML$Sample2)

SorML <- SorML[SorML$Ind1 < SorML$Ind2, ]


SorML$Sample1<-as.character(SorML$Sample1)
SorML$Sample2<-as.character(SorML$Sample2)

SorML$Ind1<-meta[SorML$Sample1,]$ind
SorML$Ind2<-meta[SorML$Sample2,]$ind

SorML<-SorML%>%arrange(Ind1,Ind2)
colnames(SorML)[3]<-"distMAG"  
  
  
SorFL <- melt(Fun.sor, 
              varnames = c("Sample1", "Sample2"), value.name = "Distance")  



SorFL$Ind1<-as.numeric(SorFL$Sample1)
SorFL$Ind2<-as.numeric(SorFL$Sample2)

SorFL <- SorFL[SorFL$Ind1 < SorFL$Ind2, ]


SorFL$Sample1<-as.character(SorFL$Sample1)
SorFL$Sample2<-as.character(SorFL$Sample2)

SorFL$Ind1<-meta[SorFL$Sample1,]$ind
SorFL$Ind2<-meta[SorFL$Sample2,]$ind

SorFL<-SorFL%>%arrange(Ind1,Ind2)
colnames(SorFL)[3]<-"distFun"  
  
  
sor<-left_join(SorFL,SorML,by=c("Sample1", "Sample2","Ind1","Ind2"))

write.table(sor,"sor.txt",sep="\t",quote=F)


cor.test(sor$distFun,sor$distMAG,method = "spearman")
cor.test(bc$distFun,bc$distMAG,method = "spearman")

sor.a<-ggplot(data = sor, aes(y=distFun, x=distMAG)) +
  geom_hex(bins = 60, # 调整分箱数量以控制粒度
           #binwidth = c(0.009, 0.003), # 替代bins参数，直接设置宽度
           aes(fill = after_stat(density)), # 使用密度值着色
           alpha = 0.9 )+    # 设置透明度
  scale_fill_viridis(option = "turbo", # 使用viridis配色方案
                     name = "Density", 
                     trans = "log", # 对数变换使低密度区域更明显
                     breaks = scales::log_breaks())+
  theme_bw()+
  theme(panel.grid.major = element_line(color="lightgrey"), panel.grid.minor = element_line(color="lightgrey"))



mean(bc$distFun)
median(bc$distFun)  # 0.2064894
quantile(bc$distFun,0.25) #0.1687355 
quantile(bc$distFun,0.75) #0.2414069 


median(bc$distMAG)  # 0.965006
quantile(bc$distMAG,0.25) # 0.9215389 
quantile(bc$distMAG,0.75) # 0.9821617 

##carbon
#carbons<-c("K01601","K00855","K15230","K15231","K01962","K02160","K01961","K01963","K00198","K14138","K00197",
#           "K00194","K15023","K14470","K14534")


carbon.fix<-read.table("D:\\Onedrive\\德州学院\\Reasearch\\肖雷雷合作\\Fix_gene-碳固定.txt",header = T,sep="\t")


KO.c<-KO[carbons,]


bc.cab<-vegan::vegdist(t(KO.c[,row.names(meta)]),method = "bray")
bc.cab<-as.matrix(bc.cab)



bcCabL <- melt(bc.cab, 
             varnames = c("Sample1", "Sample2"), value.name = "Distance")

bcCabL$Ind1<-as.numeric(bcCabL$Sample1)
bcCabL$Ind2<-as.numeric(bcCabL$Sample2)

bcCabL <- bcCabL[bcCabL$Ind1 < bcCabL$Ind2, ]


bcCabL$Sample1<-as.character(bcCabL$Sample1)
bcCabL$Sample2<-as.character(bcCabL$Sample2)

bcCabL$Ind1<-meta[bcCabL$Sample1,]$ind
bcCabL$Ind2<-meta[bcCabL$Sample2,]$ind

bcCabL<-bcCabL%>%arrange(Ind1,Ind2)
colnames(bcCabL)[3]<-"distCarbon"

bc2<-left_join(bcCabL,bcML,by=c("Sample1", "Sample2","Ind1","Ind2"))



ggplot(data = bc2, aes(y=distCarbon, x=distMAG)) +
  geom_hex(bins = 150, # 调整分箱数量以控制粒度
           binwidth = c(0.01, 0.005), # 替代bins参数，直接设置宽度
           aes(fill = after_stat(density)), # 使用密度值着色
           alpha = 0.9 )+    # 设置透明度
  scale_fill_viridis(option = "turbo", # 使用viridis配色方案
                     name = "Density", 
                     trans = "log", # 对数变换使低密度区域更明显
                     breaks = scales::log_breaks())+
  theme_bw()+
  theme(panel.grid.major = element_line(color="lightgrey"), panel.grid.minor = element_line(color="lightgrey"))

Cab.presence<-ifelse(KO.c[,row.names(meta)] > 0, 1, 0)
Cab.sor <- vegan::betadiver(t(Cab.presence), method = "sor")  
Cab.sor<-as.matrix(Cab.sor) 

SorCab <- melt(Cab.sor, 
              varnames = c("Sample1", "Sample2"), value.name = "Distance")  

SorCab$Ind1<-as.numeric(SorCab$Sample1)
SorCab$Ind2<-as.numeric(SorCab$Sample2)

SorCab <- SorCab[SorCab$Ind1 < SorCab$Ind2, ]


SorCab$Sample1<-as.character(SorCab$Sample1)
SorCab$Sample2<-as.character(SorCab$Sample2)

SorCab$Ind1<-meta[SorCab$Sample1,]$ind
SorCab$Ind2<-meta[SorCab$Sample2,]$ind

SorCab<-SorCab%>%arrange(Ind1,Ind2)
colnames(SorCab)[3]<-"distCarbon"  


sor2<-left_join(SorCab,SorML,by=c("Sample1", "Sample2","Ind1","Ind2"))

cor.test(sor$distFun,sor$distMAG,method = "pearson")
# r= 0.647 p< 0.0001
cor.test(bc$distFun,bc$distMAG,method = "pearson")

cor.test(sor2$distCarbon,sor2$distMAG,method = "pearson")
# r= 0.647 p< 0.0001
cor.test(bc2$distCarbon,bc$distMAG,method = "pearson")

ggplot(data = sor2, aes(y=distCarbon, x=distMAG)) +
  geom_hex(bins = 150, # 调整分箱数量以控制粒度
           binwidth = c(0.009, 0.003), # 替代bins参数，直接设置宽度
           aes(fill = after_stat(density)), # 使用密度值着色
           alpha = 0.9 )+    # 设置透明度
  scale_fill_viridis(option = "turbo", # 使用viridis配色方案
                     name = "Density", 
                     trans = "log", # 对数变换使低密度区域更明显
                     breaks = scales::log_breaks())+
  theme_bw()+
  theme(panel.grid.major = element_line(color="lightgrey"), panel.grid.minor = element_line(color="lightgrey"))


car2<-read.table("Carbons.TPM",header = T,sep="\t")

bc.cab2<-vegan::vegdist(t(car2[,row.names(meta)]),method = "bray")
bc.cab2<-as.matrix(bc.cab2)



bcCabL2 <- melt(bc.cab2, 
               varnames = c("Sample1", "Sample2"), value.name = "Distance")

bcCabL2$Ind1<-as.numeric(bcCabL2$Sample1)
bcCabL2$Ind2<-as.numeric(bcCabL2$Sample2)

bcCabL2 <- bcCabL2[bcCabL2$Ind1 < bcCabL2$Ind2, ]


bcCabL2$Sample1<-as.character(bcCabL2$Sample1)
bcCabL2$Sample2<-as.character(bcCabL2$Sample2)

bcCabL2$Ind1<-meta[bcCabL2$Sample1,]$ind
bcCabL2$Ind2<-meta[bcCabL2$Sample2,]$ind

bcCabL2<-bcCabL2%>%arrange(Ind1,Ind2)
colnames(bcCabL2)[3]<-"distCarbon"

bc2.2<-left_join(bcCabL2,bcML,by=c("Sample1", "Sample2","Ind1","Ind2"))



ggplot(data = bc2.2, aes(y=distCarbon, x=distMAG)) +
  geom_hex(bins = 150, # 调整分箱数量以控制粒度
           binwidth = c(0.01, 0.005), # 替代bins参数，直接设置宽度
           aes(fill = after_stat(density)), # 使用密度值着色
           alpha = 0.9 )+    # 设置透明度
  scale_fill_viridis(option = "turbo", # 使用viridis配色方案
                     name = "Density", 
                     trans = "log", # 对数变换使低密度区域更明显
                     breaks = scales::log_breaks())+
  theme_bw()+
  theme(panel.grid.major = element_line(color="lightgrey"), panel.grid.minor = element_line(color="lightgrey"))

Cab.presence2<-ifelse(car2[,row.names(meta)] > 0, 1, 0)
Cab.sor2 <- vegan::betadiver(t(Cab.presence2), method = "sor")  
Cab.sor2<-as.matrix(Cab.sor2) 

SorCab2 <- melt(Cab.sor2, 
               varnames = c("Sample1", "Sample2"), value.name = "Distance")  

SorCab2$Ind1<-as.numeric(SorCab2$Sample1)
SorCab2$Ind2<-as.numeric(SorCab2$Sample2)

SorCab2 <- SorCab2[SorCab2$Ind1 < SorCab2$Ind2, ]


SorCab2$Sample1<-as.character(SorCab2$Sample1)
SorCab2$Sample2<-as.character(SorCab2$Sample2)

SorCab2$Ind1<-meta[SorCab2$Sample1,]$ind
SorCab2$Ind2<-meta[SorCab2$Sample2,]$ind

SorCab2<-SorCab2%>%arrange(Ind1,Ind2)
colnames(SorCab2)[3]<-"distCarbon"  


sor2.2<-left_join(SorCab2,SorML,by=c("Sample1", "Sample2","Ind1","Ind2"))

cor.test(sor$distFun,sor$distMAG,method = "pearson")
# r= 0.647 p< 0.0001
cor.test(bc$distFun,bc$distMAG,method = "pearson")

cor.test(sor2$distCarbon,sor2$distMAG,method = "pearson")
# r= 0.647 p< 0.0001
cor.test(bc2$distCarbon,bc2$distMAG,method = "pearson")

cor.test(sor2.2$distCarbon,sor2.2$distMAG,method = "pearson")
# r= 0.647 p< 0.0001
cor.test(bc2.2$distCarbon,bc2.2$distMAG,method = "pearson")

ggplot(data = sor2.2, aes(y=distCarbon, x=distMAG)) +
  geom_hex(bins = 150, # 调整分箱数量以控制粒度
           binwidth = c(0.009, 0.003), # 替代bins参数，直接设置宽度
           aes(fill = after_stat(density)), # 使用密度值着色
           alpha = 0.9 )+    # 设置透明度
  scale_fill_viridis(option = "turbo", # 使用viridis配色方案
                     name = "Density", 
                     trans = "log", # 对数变换使低密度区域更明显
                     breaks = scales::log_breaks())+
  theme_bw()+
  theme(panel.grid.major = element_line(color="lightgrey"), panel.grid.minor = element_line(color="lightgrey"))



carbon.mag<-read.table("carbon_MAG.txt",header = F,sep="\t")

MAGs.carb<-intersect(unique(carbon.mag$V1),row.names(abund))


bc.MAG.carb<-vegan::vegdist(t(abund[MAGs.carb,row.names(meta)]),method = "bray",na.rm=T)
bc.MAG.carb<-as.matrix(bc.MAG.carb)

bcML.carb <- melt(bc.MAG.carb, 
             varnames = c("Sample1", "Sample2"), value.name = "Distance")

bcML.carb$Ind1<-as.numeric(bcML.carb$Sample1)
bcML.carb$Ind2<-as.numeric(bcML.carb$Sample2)

bcML.carb <- bcML.carb[bcML.carb$Ind1 < bcML.carb$Ind2, ]


bcML.carb$Sample1<-as.character(bcML.carb$Sample1)
bcML.carb$Sample2<-as.character(bcML.carb$Sample2)

bcML.carb$Ind1<-meta[bcML.carb$Sample1,]$ind
bcML.carb$Ind2<-meta[bcML.carb$Sample2,]$ind

bcML.carb<-bcML.carb%>%arrange(Ind1,Ind2)
colnames(bcML.carb)[3]<-"distMAG"

bc2.3<-left_join(bcCabL2,bcML.carb,by=c("Sample1", "Sample2","Ind1","Ind2"))



bc.b<-ggplot(data = bc2.3, aes(y=distCarbon, x=distMAG)) +
  geom_hex(bins = 60, # 调整分箱数量以控制粒度
           #binwidth = c(0.01, 0.005), # 替代bins参数，直接设置宽度
           aes(fill = after_stat(density)), # 使用密度值着色
           alpha = 0.9 )+    # 设置透明度
  scale_fill_viridis(option = "turbo", # 使用viridis配色方案
                     name = "Density", 
                     trans = "log", # 对数变换使低密度区域更明显
                     breaks = scales::log_breaks())+
  theme_bw()+
  theme(panel.grid.major = element_line(color="lightgrey"), panel.grid.minor = element_line(color="lightgrey"))


MAG.carb.presence<-ifelse(abund[MAGs.carb,row.names(meta)] > 0, 1, 0)
MAG.carb.sor <- vegan::betadiver(t(MAG.carb.presence), method = "sor")  
MAG.carb.sor<-as.matrix(MAG.carb.sor)  



SorML.carb <- melt(MAG.carb.sor, 
              varnames = c("Sample1", "Sample2"), value.name = "Distance")  

SorML.carb$Ind1<-as.numeric(SorML.carb$Sample1)
SorML.carb$Ind2<-as.numeric(SorML.carb$Sample2)

SorML.carb <- SorML.carb[SorML.carb$Ind1 < SorML.carb$Ind2, ]


SorML.carb$Sample1<-as.character(SorML.carb$Sample1)
SorML.carb$Sample2<-as.character(SorML.carb$Sample2)

SorML.carb$Ind1<-meta[SorML.carb$Sample1,]$ind
SorML.carb$Ind2<-meta[SorML.carb$Sample2,]$ind

SorML.carb<-SorML.carb%>%arrange(Ind1,Ind2)
colnames(SorML.carb)[3]<-"distMAG"  


sor2.3<-left_join(SorCab2,SorML.carb,by=c("Sample1", "Sample2","Ind1","Ind2"))

cor.test(sor$distFun,sor$distMAG,method = "pearson")
# r= 0.647 p< 0.0001
cor.test(bc$distFun,bc$distMAG,method = "pearson")

cor.test(sor2$distCarbon,sor2$distMAG,method = "pearson")
# r= 0.647 p< 0.0001
cor.test(bc2$distCarbon,bc2$distMAG,method = "pearson")

cor.test(sor2.2$distCarbon,sor2.2$distMAG,method = "pearson")
# r= 0.647 p< 0.0001
cor.test(bc2.2$distCarbon,bc2.2$distMAG,method = "pearson")

sor.b<-ggplot(data = sor2.3, aes(y=distCarbon, x=distMAG)) +
  geom_hex(bins = 60, # 调整分箱数量以控制粒度
           #binwidth = c(0.003, 0.003), # 替代bins参数，直接设置宽度
           aes(fill = after_stat(density)), # 使用密度值着色
           alpha = 1 )+    # 设置透明度
  scale_fill_viridis(option = "turbo", # 使用viridis配色方案
                     name = "Density", 
                     trans = "log", # 对数变换使低密度区域更明显
                     breaks = scales::log_breaks())+
  theme_bw()+
  theme(panel.grid.major = element_line(color="lightgrey"), panel.grid.minor = element_line(color="lightgrey"))

write.table(sor2.3,"sorCarbon.txt",sep="\t",quote=F)

ggpubr::ggarrange(sor.a,sor.b)
ggpubr::ggarrange(bc.a,bc.b)

#8X18


median(sor$distFun,na.rm=T)#0.880
median(sor$distMAG,na.rm=T)#0.538
median(sor2.3$distCarbon,na.rm=T)#0.178
median(sor2.3$distMAG,na.rm=T)#0.547

quantile(sor$distFun,probs = c(0.25,0.75),na.rm=T)
0.867-0.892
quantile(sor$distMAG,probs = c(0.25,0.75),na.rm=T)
0.477-0.600
quantile(sor2.3$distCarbon,probs = c(0.25,0.75),na.rm=T)
0.112-0.268
quantile(sor2.3$distMAG,probs = c(0.25,0.75),na.rm=T)
0.484-0.607

aa<-data.frame(a1=sor2.3$distMAG,a2=sor$distMAG)

save.image("Density.bc.Rdata")
