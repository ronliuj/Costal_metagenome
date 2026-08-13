##根据网络，重新构建类似相关系数的cor矩阵用于ggnClusterNet分析, 重建的网络分析会和之前的不太一样。


library(igraph)
library(ggClusterNet)
library(ggplot2)
library(ggraph)

#载入预置脚本
source("edToCor.R")
source("model_maptreeLJ.R")

Oneyear<-read_graph("MF.graphml",format = "graphml")
Oneyear<-as_undirected(Oneyear,mode = "each")
plot(Oneyear)

bin.cyto<-read.table("MF.node.txt",header=T,sep="\t",quote="")
anno<-read.table("anno.txt",header = T,sep="\t") 
node.fam<-merge(bin.cyto,anno,by="ID")

binMF<-read.table("bin.MF.txt",header = T,sep="\t")
binMF$mean<-rowMeans(binMF,na.rm=T)



aa<-edge.attributes(Oneyear)

edges<-data.frame(aa)

sum(edges$weight<0)# 0

edges$weight[edges$weight==0]<-5e-16

## 将互斥的weight改为负值，作为下边cor里边的负相关系数
edges[edges$interactionType =="mutualExclusion",]$weight<-edges[edges$interactionType =="mutualExclusion",]$weight*(-1)

##根据构建好的网络边重建类距离矩阵用于ggClusterNet

corMF<-edToCor(edges)



## map_tree2 布局, use modified model_maptreeLJ cluster_edge_betweenness
result2 = model_maptreeLJ(cor = corMF,
                        method = "cluster_edge_betweenness")

# result2 = PolygonRrClusterG (cor = cor,nodeGroup =group2 )
node = result2[[1]]##节点信息
netClu = result2[[2]]##模块信息
cc<-result2[[3]]##绘图按照模块聚焦点信息
#this is the new One year
Mygraph<-result2[[4]]##重构网络
fc<-result2[[5]]##网络模块化分析结果

netClu<-arrange(netClu,ID)#16个module

membership<- membership(fc)
modularity <- modularity(Mygraph, membership = membership)   #0.618
# ---node节点注释#-----------
#nodes = nodeadd(plotcord =node,otu_table = otu_table,tax_table = tax_table)
#-----计算边#--------
nodes2 = node %>% inner_join(netClu,by = c("elements" = "ID"))
nodes2$mean<-binMF[nodes2$elements,]$mean
nodes2<-data.frame(nodes2)
row.names(nodes2)<-nodes2$elements
#nodes2$group<-paste("model_",nodes2$group,sep="")#nodes2 点的坐标

length(unique(nodes2$group))#16
table(nodes2$group)

edge = edgeBuild(cor = corMF,node = node)#获得边的坐标


#颜色设定根据实际情况来

cols <-  colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(16)
names(cols)<-unique(nodes2$group)
cols<-c(cols,"others"="lightgrey")

edge$model1<-nodes2[edge$OTU_1,]$group
edge$model2<-nodes2[edge$OTU_2,]$group


edge2 = edge %>% mutate(color = ifelse(model1 == model2,model1,"others"))

write.table(edge2,"networkmf.txt",sep="\t",quote=F)


### 出图
pnet <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor)),
                                data = edge, size = 0.5) +
  geom_point(aes(X1, X2,fill = group,size = log(mean)),pch = 21, data = nodes2) +
  scale_fill_manual(values=cols) +
  scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
  # labs( title = paste(layout,"network",sep = "_"))+
  # geom_text_repel(aes(X1, X2,label=Phylum),size=4, data = plotcord)+
  # discard default grid + titles in ggplot2
  theme(panel.background = element_blank()) +
  # theme(legend.position = "none") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.background = element_rect(colour = NA)) +
  theme(panel.background = element_rect(fill = "white",  colour = NA)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
pnet

pnet2 <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = color),
                                         data = edge2[edge2$color=="others",], size = 0.5) +
  geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = color),
               data = edge2[edge2$color!="others",], size = 0.5) +
  geom_point(aes(X1, X2,fill = group,color="lightgrey",size = log(mean)),pch = 21, data = nodes2) +
  scale_fill_manual(values = cols)+scale_color_manual(values = cols)+
  scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
  # labs( title = paste(layout,"network",sep = "_"))+
  # geom_text_repel(aes(X1, X2,label=Phylum),size=4, data = plotcord)+
  # discard default grid + titles in ggplot2
  theme(panel.background = element_blank()) +
  # theme(legend.position = "none") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.background = element_rect(colour = NA)) +
  theme(panel.background = element_rect(fill = "white",  colour = NA)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
pnet2

ggpubr::ggarrange(pnet,pnet2,nrow=1,ncol=2)
#export 6X15

## map_tree2 布局, fast greedy 聚类

result3 = model_maptreeLJ(cor = corMF,
                          method = "cluster_fast_greedy")

# result3 = PolygonRrClusterG (cor = cor,nodeGroup =group2 )
node2 = result3[[1]]##节点信息
netClu2 = result3[[2]]##模块信息
dd<-result3[[3]]##绘图按照模块聚焦点信息
#this is the new One year
Mygraph2<-result3[[4]]##重构网络
fc2<-result3[[5]]##网络模块化分析结果

netClu2<-arrange(netClu2,ID)#37个module

length(unique(netClu2$group))#37
modnum<-data.frame(table(netClu2$group))
modnum<-arrange(modnum,desc(Freq))

##只显示前15个点多于6个的moldule

modsel<-modnum[1:15,]$Var1

membership2<- membership(fc2)
modularity2 <- modularity(Mygraph2, membership = membership2)   #0.38
# ---node2节点注释#-----------
#nodes = node2add(plotcord =node2,otu_table = otu_table,tax_table = tax_table)
#-----计算边#--------
nodes3 = node2 %>% inner_join(netClu2,by = c("elements" = "ID"))
nodes3$mean<-binMF[nodes3$elements,]$mean
nodes3<-data.frame(nodes3)
row.names(nodes3)<-nodes3$elements
#nodes3$group<-paste("model_",nodes3$group,sep="")#nodes3 点的坐标

length(unique(nodes3$group))#37
table(nodes3$group)

edge3 = edgeBuild(cor = corMF,node = node2)#获得边的坐标


#颜色设定根据实际情况来

cols2 <-  colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(15)
names(cols2)<-modsel
cols2<-c(cols2,"others"="lightgrey")



##只显示前15个，其他均改为others

nodes3<-nodes3 %>% mutate(color = ifelse(is.element(group,modsel),group,"others"))

edge3$model1<-nodes3[edge3$OTU_1,]$color
edge3$model2<-nodes3[edge3$OTU_2,]$color

edge4 = edge3 %>% mutate(color = ifelse(model1 == model2,model1,"others"))  
     

### 出图
pnet3 <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor)),
                                data = edge3, size = 0.5) +
  geom_point(aes(X1, X2,fill = color,size = log(mean)),pch = 21, data = nodes3) +
  scale_fill_manual(values=cols2) +
  scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
  # labs( title = paste(layout,"network",sep = "_"))+
  # geom_text_repel(aes(X1, X2,label=Phylum),size=4, data = plotcord)+
  # discard default grid + titles in ggplot2
  theme(panel.background = element_blank()) +
  # theme(legend.position = "none") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.background = element_rect(colour = NA)) +
  theme(panel.background = element_rect(fill = "white",  colour = NA)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
pnet3

pnet4 <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = color),
                                 data = edge4[edge4$color=="others",], size = 0.5) +
  geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = color),
               data = edge4[edge4$color!="others",], size = 0.5) +
  geom_point(aes(X1, X2,fill = color,color="lightgrey",size = log(mean)),pch = 21, data = nodes3) +
  scale_fill_manual(values = cols2)+scale_color_manual(values = cols2)+
  scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
  # labs( title = paste(layout,"network",sep = "_"))+
  # geom_text_repel(aes(X1, X2,label=Phylum),size=4, data = plotcord)+
  # discard default grid + titles in ggplot2
  theme(panel.background = element_blank()) +
  # theme(legend.position = "none") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.background = element_rect(colour = NA)) +
  theme(panel.background = element_rect(fill = "white",  colour = NA)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
pnet4

pnet5 <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = color),
                                 data = edge4[edge4$color!="others",], size = 0.5) + ##不显示无关连接
  geom_point(aes(X1, X2,fill = color,color="lightgrey",size = log(mean)),pch = 21, data = nodes3) +
  scale_fill_manual(values = cols2)+scale_color_manual(values = cols2)+
  scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
  # labs( title = paste(layout,"network",sep = "_"))+
  # geom_text_repel(aes(X1, X2,label=Phylum),size=4, data = plotcord)+
  # discard default grid + titles in ggplot2
  theme(panel.background = element_blank()) +
  # theme(legend.position = "none") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.background = element_rect(colour = NA)) +
  theme(panel.background = element_rect(fill = "white",  colour = NA)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
pnet5


ggpubr::ggarrange(pnet3,pnet4,pnet5,nrow=3,ncol=1)
#export 6X12

save.image("try2.Rdata")
