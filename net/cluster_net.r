#（二）MF 2024.11.7
#1.Cytoscape运行结果数据导入与整理——核心数据框Oneyear与att
#install.packages("igraph")
#install.packages("RColorBrewer")

library(igraph)
library(ggClusterNet)
library(ggplot2)
library(ggraph)

Oneyear<-read_graph("MF.graphml",format = "graphml")
Oneyear<-as.undirected(Oneyear,mode = "each")
plot(Oneyear)

bin.cyto<-read.table("MF.node.txt",header=T,sep="\t",quote="")
anno<-read.table("anno.txt",header = T,sep="\t") 
node.fam<-merge(bin.cyto,anno,by="ID")
#整理数据，给node.fam添加key(ID)  Color (Phylum)两列 
node.fam$Color<-node.fam$Phylum
node.fam$key<-node.fam$ID
#增加bin丰度均值信息到node.fam数据框
abund<-bin.MF
abund$mean<-rowMeans(abund)        
abund[,36]<-row.names(abund) 
colnames(abund)[36]<-"ID"
node.fam<-merge(node.fam,abund,by="ID")
node.fam<-node.fam[,-c(24:57)]
row.names(node.fam)<-node.fam$ID

att<-node.fam
row.names(att)<-att$key 
#整理数据索引列key
ver<-data.frame(V(Oneyear))
ver$name<-row.names(ver)
library(tidyr)
ver<-separate(ver,name,into=c("o","ID"),sep=4)
ver<-ver[,c(1,3)]
#ver$ID<-sub("-","_",ver$ID) #将ID列中所有的-字符替换成_字符。    
key<-ver$ID   

#2. 将domian  phylum  key;  abund 添加到Oneyear，只能一列一列加(注意行名对应起来)
###特别注意：att需要索引的全部列必须在索引key创建前整理好，否则后来添加的列无效；
V(Oneyear)$Yanse<-att[key,]$Color
V(Oneyear)$domian<-att[key,]$Domain
V(Oneyear)$phylum<-att[key,]$Color
V(Oneyear)$abundance<-att[key,]$mean
V(Oneyear)$key<-att[key,]$key
V(Oneyear)$abundance     #检查合并结果
V(Oneyear)$domian        #检查合并结果
V(Oneyear)$phylum


#2.1 颜色设置 phylum
unique(att$Phylum) #检查门的数量，MF共12个
library(RColorBrewer)
colors<-brewer.pal(12,"Paired")    #该调色板最多12种颜色
colors.point<-c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C",
                "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928")  
point.type<-c("Pseudomonadota",  "Bacteroidota" , "Bacillota", "Gemmatimonadota",  
              "Campylobacterota","Actinomycetota","Chloroflexota","Desulfobacterota", 
              "Acidobacteriota", "Patescibacteria","Verrucomicrobiota","Bdellovibrionota")

colors.point<-c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#AB7519",
                "#CAB2D6","#6A3D9A","#FFFF99","#B15928","#E64E00","#E6EB00","#65B48E","#3E5CC5","#FB1299",
                "#F3C846","#D882AD","#DDDEDE","#AAAB20")
point.type<-c("Acidobacteriota","Actinobacteriota","Armatimonadota","Bacteroidota","Bdellovibrionota",
              "Chloroflexota", "Cyanobacteria","Desulfobacterota","Elusimicrobiota","Fibrobacterota",
              "Firmicutes","Gemmatimonadota", "Halobacteriota","Margulisbacteria","Methanobacteriota",
              "Methylomirabilota","Myxococcota","Nitrospirota","Patescibacteria","Planctomycetota",
              "Proteobacteria","Thermoproteota")

point.type<-as.factor(point.type)
names(colors.point)<-levels(point.type)
colors.point        #检查对应关系   色板设置，仅运行一次即可

V(Oneyear)$color <- colors.point[V(Oneyear)$phylum]   
V(Oneyear)$color      


#2.2形状 domian
shapes<-c("circle","square")
Xingz<-att[key,]$Domain
Xingz<-as.factor(Xingz)
names(shapes)<-levels(Xingz)
shapes       #检查对应关系
V(Oneyear)$shape <- shapes[V(Oneyear)$domian] 
V(Oneyear)$Label <- V(Oneyear)$key    

#2.3 连线颜色设置
colors.ed<-c("#66C2A5","#FC8D62")
ed.type<-E(Oneyear)$interactionType
ed.type<-as.factor(ed.type)
names(colors.ed)<-levels(ed.type)
E(Oneyear)$color <- colors.ed[E(Oneyear)$interactionType]
colors.ed        #检查对应关系
#copresence mutualExclusion 
#"#66C2A5"       "#FC8D62" 

#2.4 节点大小及连线宽度设置
V(Oneyear)$size <- log(V(Oneyear)$abundance*100)
E(Oneyear)$width <- abs(log10(E(Oneyear)$weight+0.0001))/2


#3. 作图
circle<- layout_in_circle(Oneyear)   
plot(Oneyear, layout=circle, vertex.label="",edge.curved=.2)   #单独画一个图
plot(Oneyear, layout=circle, vertex.label="")   #单独画一个图

sphere<- layout_on_sphere(Oneyear)
plot(Oneyear, layout=sphere, vertex.label="",edge.curved=.2)

gem<-layout_with_gem(Oneyear)
plot(Oneyear, layout=gem, vertex.label="",edge.curved=.2) #8*8出图

library(igraph)
#基于点链接的社群划分
clusters(Oneyear,mode="strong")
components(Oneyear,mode="strong")
#随机游走-社群划分
member<-walktrap.community(Oneyear,weights=E(Oneyear)$weight,step=4)

#中间中心度-社群划分
# 检查权重的非零最小值
min_nonzero_value <- min(E(Oneyear)$weight[E(Oneyear)$weight != 0], na.rm = TRUE)    #查找非零最小值
#2.81e-15     添加最小值2.2e-16
# 如果权重中有非正数，将其替换为正数

E(Oneyear)$weight<-abs(E(Oneyear)$weight)
E(Oneyear)$weight[E(Oneyear)$weight==0]<-2.2e-16

fc<-cluster_edge_betweenness(Oneyear,weight=E(Oneyear)$weight,directed=F)   #directed=T时就代表有向线。
print(fc)  #groups: 19, mod: 0.61  模块度modularity
#警告：在执行基于边介数中心性的社区划分时，算法计算了多个可能的社区划分，
#并从中选择了模块度（modularity）最高的一个作为最终结果。
#模块度是衡量社区划分质量的一个指标，它反映了社区内部节点连接的紧密程度相对于社区间节点连接的紧密程度。
#membership <- cluster_edge_betweenness(Oneyear, directed = FALSE, weights = E(Oneyear)$weight)
#cluster_fast_greedy 马尔科夫
#模块化指标Q——modularity


V(Oneyear)$modularity <- membership(fc)%>% as.numeric()
V(Oneyear)$label <- V(Oneyear)$name
V(Oneyear)$label <- NA
modu_sort <- V(Oneyear)$modularity %>% table() %>% sort(decreasing = T)
Top_M = 20
top_num <- Top_M
modu_name <- names(modu_sort[1:Top_M])
cols <-  colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(Top_M)
#modu_cols <- cols[1:length(modu_name)]
#names(modu_cols) <- modu_name
V(Oneyear)$color <- V(Oneyear)$modularity
col_g <- "#C1C1C1"
#vertex_attr(Oneyear)查看全部性质

V(Oneyear)$color[!(V(Oneyear)$color %in% modu_name)] <- col_g
V(Oneyear)$color[(V(Oneyear)$color %in% modu_name)] <- modu_cols[match(V(Oneyear)$color[(V(Oneyear)$color %in% modu_name)],modu_name)]
V(Oneyear)$frame.color <- V(Oneyear)$color

E(Oneyear)$color <- col_g
for ( i in modu_name){
  col_edge <- cols[which(modu_name==i)]
  otu_same_modu <-V(Oneyear)$name[which(V(Oneyear)$modularity==i)]
  E(Oneyear)$color[(data.frame(as_edgelist(Oneyear))$X1 %in% otu_same_modu)&(data.frame(as_edgelist(Oneyear))$X2 %in% otu_same_modu)] <- col_edge
}

Oneyear.degree<-igraph::degree(Oneyear) %>% as.data.frame()
colnames(Oneyear.degree) = "degree"
Oneyear.degree$ID = row.names(Oneyear.degree)

netClu = data.frame(ID = names(membership(fc)),group = as.vector(membership(fc)))
dim(netClu)
table(netClu$group)
netClu$group = paste("model_",netClu$group,sep = "")

netClu <- netClu %>% full_join(Oneyear.degree,na_matches = "never")
netClu$degree[is.na(netClu$degree)] = 0
netClu <- netClu %>%
  dplyr::arrange(desc(degree))

#netClu = netClu %>% full_join(data.frame(ID = row.names(cor)),na_matches = "never")
#head(netClu)
netClu2 = netClu
netClu2$group[is.na(netClu2$group)] = paste("mother_","no",sep = "")
netClu2$degree[is.na(netClu2$degree)] = 0

##如果存在单独的点情况下运行一下代码（参加ggClusterNet--modle_maptree2
# tem = netClu$group[is.na(netClu$group)]
# netClu$group[is.na(netClu$group)] = paste("other_",1: length(tem),sep = "")
# netClu$degree[is.na(netClu$degree)] = 0

edge =  data.frame(model = paste(netClu$group,sep = ""),OTU = netClu$ID)
head(edge)
colnames(edge) = c("from","to")

vertices_t  <-  data.frame(
  name = unique(c(as.character(edge$from), as.character(edge$to))))
head(vertices_t)
vertices_t$size = sample(1:10,nrow(vertices_t),replace = TRUE)

##根据module进行聚集分布，整体分布在一个圆形上（利用module与node重新构建网络，然后使用circlepack实现！！！！
mygraph <- igraph::graph_from_data_frame(edge, vertices= vertices_t )
#-----------------------------------设置颜色映射参数-------------------------
# ?create_layout
# ,weight = mean, sort.by = NULL, direction = "out"
data = ggraph::create_layout(mygraph, layout = 'circlepack',weight = size)
head(data)

node = data %>% dplyr::filter(leaf == TRUE ) %>%
  dplyr::select(x,y,name)
colnames(node) = c("X1","X2", "elements")
row.names(node) = node$elements

branch = data %>% dplyr::filter(leaf != TRUE ) %>%
  dplyr::select(x,y,name)
colnames(branch) = c("X1","X2", "elements")
row.names(branch) = branch$elements
colnames(branch)[1:2] = c("x","y")
branch$elements = gsub("model_","",branch$elements)
row.names(branch) = branch$elements

##node为分布情况，netClu保存分组信息，branches是分组节点位置，作图不需要




# ---node节点注释#-----------


nodes2 = node %>% inner_join(netClu,by = c("elements" = "ID"))
n=length(unique(nodes2$group))
nodes2<-nodes2 %>% separate(elements,into=c("OTU","ID"),sep="-",remove=F)

cols <-  colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(19)
names(cols)<-unique(nodes2$group)
cols<-c(cols,"others"="lightgrey")

nodes2$mean<-node.fam[nodes2$ID,]$mean

#-----计算边#--------
Edges<-edge_attr(Oneyear)
Edges<-data.frame(Edges)

Edges<-Edges %>% separate(shared.name, c('OTU_1', 'OTU_2'),sep="->",remove=F)

Edges<-Edges %>% dplyr::select(SUID, weight, OTU_1, OTU_2)

tem2 = nodes2 %>% 
  dplyr::select(elements,group,X1,X2) %>%
  dplyr::right_join(Edges,by =c("elements" = "OTU_1" ) ) %>%
  dplyr::rename(OTU_1 = elements,model1 =group,Y1=X2)
head(tem2)

tem3 = nodes2 %>% 
  dplyr::select(elements,group,X1,X2) %>%
  dplyr::right_join(Edges,by =c("elements" = "OTU_2" ) ) %>%
  dplyr::rename(OTU_2 = elements,model2 =group,X2=X1,Y2=X2)
head(tem3)

tem4 = tem2 %>%inner_join(tem3)
head(tem4)

edge2 = tem4 %>% mutate(color = ifelse(model1 == model2,model1,"others"))

### 出图
pnet <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = color),
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
pnet
#export 6x7

save.image("net-try-1130.Rdata")






## igraph layout
# 计算 layout
sub_net_layout <- layout_with_fr(Oneyear, niter=999,grid = 'nogrid')
data = as.data.frame(sub_net_layout)
data$OTU = igraph::vertex_attr(Oneyear)$name
colnames(data) = c("X1","X2","elements")

tem =  V(Oneyear)$modularity
tem[!tem %in% modu_name] = "mini_model"
tem[tem %in% modu_name] = paste("model_",tem[tem %in% modu_name],sep = "")

row.names(data) = data$elements
dat = data.frame(orig_model = V(Oneyear)$modularity,
                 model = tem,
                 color = V(Oneyear)$color,
                 OTU = igraph::vertex_attr(Oneyear)$name,
                 X1 = data$X1,
                 X2 = data$X2)

result2<-list(data,dat,Oneyear)

node = result2[[1]]
head(node)

dat = result2[[2]]
head(dat)
tem = data.frame(mod = dat$model,col = dat$color) %>%  
  dplyr::distinct( mod, .keep_all = TRUE)  
col = tem$col
names(col) = tem$mod

#---node节点注释#-----------
# otu_table = as.data.frame(t(vegan_otu(ps)))
# tax_table = as.data.frame(vegan_tax(ps))
# nodes = nodeadd(plotcord =node,otu_table = otu_table,tax_table = tax_table)
# head(nodes)
# #-----计算边#--------
# edge = edgeBuild(cor = cor,node = node)
# colnames(edge)[8] = "cor"
# head(edge)

Edges<-edge_attr(Oneyear)
Edges<-data.frame(Edges)

Edges<-Edges %>% separate(shared.name, c('OTU_1', 'OTU_2'),sep="->",remove=F)

Edges<-Edges %>% dplyr::select(SUID, weight, OTU_1, OTU_2)

tem2 = dat %>% 
  dplyr::select(OTU,model,color,X1,X2) %>%
  dplyr::right_join(Edges,by =c("OTU" = "OTU_1" ) ) %>%
  dplyr::rename(OTU_1 = OTU,model1 = model,color1 = color,Y1=X2)
head(tem2)

tem3 = dat %>% 
  dplyr::select(OTU,model,color,X1,X2) %>%
  dplyr::right_join(Edges,by =c("OTU" = "OTU_2" ) ) %>%
  dplyr::rename(OTU_2 = OTU,model2 = model,color2 = color,X2=X1,Y2=X2)
head(tem3)

tem4 = tem2 %>%inner_join(tem3)
head(tem4)

edge2 = tem4 %>% mutate(color = ifelse(model1 == model2,as.character(model1),"across"),
                        manual = ifelse(model1 == model2,as.character(color1),"#C1C1C1")
)

head(edge2)
col_edge = edge2 %>% dplyr::distinct(color, .keep_all = TRUE)  %>% dplyr::select(color,manual)
col0 = col_edge$manual
names(col0) = col_edge$color

library(ggnewscale)

p1 <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = color),
                              data = edge2, linewidth = 1) +
  scale_colour_manual(values = col0) 

# ggsave("./cs1.pdf",p1,width = 16,height = 14)
p2 = p1 +
  new_scale_color() +
  geom_point(aes(X1, X2,color =model), data = dat,size = 4) +
  scale_colour_manual(values = col) +
  scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
  theme(panel.background = element_blank()) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.background = element_rect(colour = NA)) +
  theme(panel.background = element_rect(fill = "white",  colour = NA)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
p2














membership<- membership(fc)
modularity <- modularity(Oneyear, membership = fc)   #0.6382

# 创建一个数据框来存储节点和它们所属的社群信息
nodes_df <- data.frame(
  Node = names(V(Oneyear)),  # 假设 g.undir 是您的图对象，并且节点有名称
  Cluster = membership
)

#把模块节点信息赋值给Oneyear
nodes_df<- separate(nodes_df,col="Node",into=c("a","key"),sep = "-",
                    remove = TRUE, convert = FALSE)
rownames(nodes_df)<-nodes_df[,2]
V(Oneyear)$modul<-nodes_df[key,]$Cluster
#给模块节点信息赋值颜色
mem.col<-rainbow(length(unique(V(Oneyear)$modul)),alpha = 0.3)
V(Oneyear)$modcol<-mem.col[V(Oneyear)$modul]

sphere<- layout_on_sphere(Oneyear)
plot(Oneyear, layout=sphere, vertex.label="",edge.curved=.2)

plot(Oneyear, layout=sphere, vertex.label="",vertex.color=V(Oneyear)$modcol,edge.curved=.2) #其位置分布没有按照模块分类
plot(member, Oneyear,vertex.size=V(Oneyear)$size,vertex.label=NA)

#画出子网络图
table(nodes_df$Cluster)
#1  2  3  4  5  6  7  8  9  10  11 12 13 14 15 16 17 18 19 
#43 21 21 16 34 34 26 14 10 16  7 22  9  9  1  2  1  3  2 
selected_vertices <- which(V(Oneyear)$modul == 1)   #将.graph数据转为向量，再转为数据框
sub_graph <- subgraph(Oneyear, selected_vertices)
plot(sub_graph, layout=layout.sphere, vertex.label="", vertex.color=V(sub_graph)$color, edge.curved=.2)


#网络聚类系数——transitivity   值越大代表交互关系越大，说明网络越复杂。
transitivity(Oneyear)    #[1] 0.6265925
#网络密度——graph.density    越大，说明网络越复杂
edge_density(Oneyear)    #[1] 0.02948217
#graph.density(Oneyear)  #从中可以看到不同社群与整体之间的网络密度情况（关联程度）
#平均最短路径
mean_distance(Oneyear)   #[1] 0.001964403



#可信度检验
install.packages("survival")
install.packages("MASS")   #未安装成功
install.packages("fitdistrplus")
library(fitdistrplus)
#生成1000个和最大社团相同节点数和密度的随机网络
result1 <- c()
for(i in 1:1000){
  g<-erdos.renyi.game(43, 0.35) #节点数和密度
  julei <- transitivity(g)
  result1 <- c(result1,julei)}
fit <- fitdist(result1,"norm")
plot(fit)
#mean和sd值，输入fit即可获取。
left<-mean-1.96*sd
left<-fit$estimate-1.96*fit$sd






#加载sql相关模块____________废
install.packages("gsubfn")
install.packages("proto")
install.packages("RSQLite")
install.packages("sqldf")
library(gsubfn)    
library(proto)
library(RSQLite)
library(sqldf)  

g1 <- graph_from_data_frame(Oneyear)
att2<-as.data.frame(V(Oneyear))
#查询图中的点
nodes_df<-nodes_df[,-1]
att1<-merge(att,nodes_df,by="key")
table(att1$Cluster)
#1  2  3  4  5  6  7  8  9  10  11 12 13 14 15 16 17 18 19 
#43 21 21 16 34 34 26 14 10 16  7 22  9  9  1  2  1  3  2 
resultmc <- sqldf("select * from att1 where Cluster=1")   #查询最大子网络
gmc <- graph_from_data_frame(resultmc,directed = FALSE)

#标签传播最大社团图
a <- data.frame(lc[2])
resultlc <- sqldf("select * from data where x in a and y in a")
glc <-  graph.data.frame(resultlc,directed = FALSE)
plot(glc,vertex.size=5,vertex.label=NA)
#GN算法最大社团图
a <- data.frame(ec[3])
resultec <- sqldf("select * from data where x in a and y in a")
gec <-  graph.data.frame(resultec,directed = FALSE)
plot(gec,vertex.size=5,vertex.label=NA)

#绘制度的图
degreemc<- degree(Oneyear,mode="all",normalized=T)
plot(degreemc)
degreelc <- degree(sub_graph,mode="all",normalized=T)
plot(table(degreelc))
degreeec<- degree(gec,mode="all",normalized=T)
plot(degreeec)
#计算度的最大值、聚类系数、平均最短路径
max(degreemc)
transitivity(Oneyear)
mean_distance(Oneyear)
max(degreelc)
transitivity(sub_graph)
edge_density(sub_graph)   
mean_distance(sub_graph)

















