#2024.1.23 数据前处理
bin<-read.table("abund.txt",header=T,sep="\t",quote="")
row.names(bin)<-bin$ID
bin<-bin[,3:118]

min_nonzero_value <- min(bin[bin != 0], na.rm = TRUE)    #查找非零最小值
# 0.0762780          0.001
bin1 <- bin[,2:21]+0.001
rownames(bin1)<-bin[,1]

write.table(bin1,"cytoscape.binCon.txt",sep="\t",row.names=T,quote=F)
#然后在excel中分开各组数据，导入Cytoscape中 进行网络模型构建


#1.数据导入与整理——核心数据框Oneyear与att
install.packages("igraph")
install.packages("RColorBrewer")

library(igraph)
Oneyear<-read_graph("Network_60.graphml",format = "graphml")
Oneyear<-read_graph("Network_Con.graphml",format = "graphml")
Oneyear<-as.undirected(Oneyear,mode = "each")
plot(Oneyear)

bin.cyto<-read.table("bin.cyto_Con.txt",header=T,sep="\t",quote="")
species.anno<-read.table("species.anno314.txt",header=T,sep="\t",quote="")
node.fam<-merge(bin.cyto,species.anno,by="Genomic.bins")
write.table(node.fam,"node.fam_Con.txt",sep="\t",row.names=T,quote=F)  #数据框太大时使用
#再次整理数据_Excel，给node.fam添加key  Color 两列 
att<-read.table("att_Con.txt",header=T,sep="\t",quote="")

#整理数据索引列key
row.names(att)<-att$key 
ver<-data.frame(V(Oneyear))
ver$name<-row.names(ver)
library(tidyr)
ver<-separate(ver,name,into=c("o","ID"),sep=4)
ver<-ver[,c(1,3)]
ver$ID<-sub("-","_",ver$ID)     
key<-ver$ID

#2. 将domian  phylum  key;  abund 添加到Oneyear，只能一列一列加(注意行名对应起来)
V(Oneyear)$name     #为查看数据行名顺序，以便下一步按顺序合并。
d.name<-att[,c(5,13,15,20)]
#增加bin丰度信息
abund1<-read.table("bin_abundance0.001.txt",header = T,sep="\t",quote="") #运行一次即可
row.names(abund1)<-abund1$Genomic.bins

#abund<-abund[,2:21]  #在所有样本中取了中值，下面考虑分组别
abund<-abund1[,18:21]    #增雨60%
abund<-abund1[,2:5]    #减雨60%
abund<-abund1[,10:13]    #对照组
abund<-abund1[,14:16]    #增雨40%
abund<-abund1[,6:9]    #减雨40%

abund$mean<-rowMeans(abund)        
abund[,6]<-row.names(abund) #前面取的列数+2
colnames(abund)[6]<-"Genomic.bins"

d.name[,5]<-rownames(d.name)
colnames(d.name)[5]<-"Genomic.bins"
d.name<-merge(abund,d.name,by="Genomic.bins")
d.name$name<-factor(d.name$name,
                    levels = c("OTU-bin.239",  "OTU-bin.374",  "OTU-bin.633",  "OTU-bin.634",  "OTU-bin.247", 
                               "OTU-bin.248",  "OTU-bin.249",  "OTU-bin.644",  "OTU-bin.647",  "OTU-bin.648", 
                               "OTU-bin.251",  "OTU-bin.1000", "OTU-bin.1002", "OTU-bin.258",  "OTU-bin.259", 
                               "OTU-bin.1006", "OTU-bin.656",  "OTU-bin.657",  "OTU-bin.629",  "OTU-bin.261", 
                               "OTU-bin.641",  "OTU-bin.262",  "OTU-bin.1014", "OTU-bin.302",  "OTU-bin.662", 
                               "OTU-bin.1017", "OTU-bin.666",  "OTU-bin.667",  "OTU-bin.704",  "OTU-bin.705", 
                               "OTU-bin.870",  "OTU-bin.793",  "OTU-bin.270",  "OTU-bin.1020", "OTU-bin.1022",
                               "OTU-bin.277",  "OTU-bin.278",  "OTU-bin.1024", "OTU-bin.311",  "OTU-bin.671", 
                               "OTU-bin.312",  "OTU-bin.675",  "OTU-bin.1170", "OTU-bin.318",  "OTU-bin.712", 
                               "OTU-bin.694",  "OTU-bin.716",  "OTU-bin.717",  "OTU-bin.718",  "OTU-bin.719", 
                               "OTU-bin.282",  "OTU-bin.1031", "OTU-bin.321",  "OTU-bin.322",  "OTU-bin.1036",
                               "OTU-bin.323",  "OTU-bin.684",  "OTU-bin.685",  "OTU-bin.687",  "OTU-bin.264", 
                               "OTU-bin.293",  "OTU-bin.1040", "OTU-bin.1042", "OTU-bin.297",  "OTU-bin.692", 
                               "OTU-bin.333",  "OTU-bin.1047", "OTU-bin.1049", "OTU-bin.697",  "OTU-bin.698", 
                               "OTU-bin.731",  "OTU-bin.699",  "OTU-bin.474",  "OTU-bin.340",  "OTU-bin.1054",
                               "OTU-bin.1055", "OTU-bin.12",   "OTU-bin.347",  "OTU-bin.740",  "OTU-bin.17",  
                               "OTU-bin.741",  "OTU-bin.742",  "OTU-bin.19",   "OTU-bin.745",  "OTU-bin.747", 
                               "OTU-bin.484",  "OTU-bin.353",  "OTU-bin.1065", "OTU-bin.352",  "OTU-bin.1067",
                               "OTU-bin.354",  "OTU-bin.1100", "OTU-bin.355",  "OTU-bin.462",  "OTU-bin.27",  
                               "OTU-bin.751",  "OTU-bin.1106", "OTU-bin.1107", "OTU-bin.757",  "OTU-bin.254", 
                               "OTU-bin.759",  "OTU-bin.360",  "OTU-bin.1075", "OTU-bin.31",   "OTU-bin.366", 
                               "OTU-bin.400",  "OTU-bin.368",  "OTU-bin.36",   "OTU-bin.1115", "OTU-bin.1116",
                               "OTU-bin.406",  "OTU-bin.766",  "OTU-bin.407",  "OTU-bin.769",  "OTU-bin.801", 
                               "OTU-bin.803",  "OTU-bin.805",  "OTU-bin.806",  "OTU-bin.30",   "OTU-bin.72",  
                               "OTU-bin.370",  "OTU-bin.40",   "OTU-bin.373",  "OTU-bin.1087", "OTU-bin.376", 
                               "OTU-bin.44",   "OTU-bin.410",  "OTU-bin.1124", "OTU-bin.46",   "OTU-bin.379", 
                               "OTU-bin.1127", "OTU-bin.1129", "OTU-bin.775",  "OTU-bin.50",   "OTU-bin.812", 
                               "OTU-bin.1090", "OTU-bin.1093", "OTU-bin.1094", "OTU-bin.1131", "OTU-bin.53",  
                               "OTU-bin.54",   "OTU-bin.56",   "OTU-bin.780",  "OTU-bin.421",  "OTU-bin.1135",
                               "OTU-bin.58",   "OTU-bin.1",    "OTU-bin.788",  "OTU-bin.821",  "OTU-bin.1119",
                               "OTU-bin.822",  "OTU-bin.824",  "OTU-bin.1096", "OTU-bin.390",  "OTU-bin.392", 
                               "OTU-bin.63",   "OTU-bin.68",   "OTU-bin.1147", "OTU-bin.1035", "OTU-bin.438", 
                               "OTU-bin.70",   "OTU-bin.74",   "OTU-bin.1153", "OTU-bin.76",   "OTU-bin.1155",
                               "OTU-bin.77",   "OTU-bin.443",  "OTU-bin.1157", "OTU-bin.79",   "OTU-bin.1159",
                               "OTU-bin.848",  "OTU-bin.1162", "OTU-bin.1166", "OTU-bin.456",  "OTU-bin.631", 
                               "OTU-bin.1038", "OTU-bin.853",  "OTU-bin.856",  "OTU-bin.857",  "OTU-bin.102", 
                               "OTU-bin.98",   "OTU-bin.1178", "OTU-bin.465",  "OTU-bin.106",  "OTU-bin.109", 
                               "OTU-bin.469",  "OTU-bin.386",  "OTU-bin.502",  "OTU-bin.862",  "OTU-bin.790", 
                               "OTU-bin.505",  "OTU-bin.507",  "OTU-bin.900",  "OTU-bin.901",  "OTU-bin.903", 
                               "OTU-bin.905",  "OTU-bin.470",  "OTU-bin.119",  "OTU-bin.513",  "OTU-bin.342", 
                               "OTU-bin.515",  "OTU-bin.516",  "OTU-bin.517",  "OTU-bin.265",  "OTU-bin.912", 
                               "OTU-bin.919",  "OTU-bin.120",  "OTU-bin.735",  "OTU-bin.123",  "OTU-bin.125", 
                               "OTU-bin.486",  "OTU-bin.1077", "OTU-bin.880",  "OTU-bin.523",  "OTU-bin.882", 
                               "OTU-bin.524",  "OTU-bin.883",  "OTU-bin.526",  "OTU-bin.527",  "OTU-bin.920", 
                               "OTU-bin.889",  "OTU-bin.923",  "OTU-bin.924",  "OTU-bin.221",  "OTU-bin.926", 
                               "OTU-bin.135",  "OTU-bin.330",  "OTU-bin.499",  "OTU-bin.532",  "OTU-bin.894", 
                               "OTU-bin.122",  "OTU-bin.896",  "OTU-bin.537",  "OTU-bin.898",  "OTU-bin.539", 
                               "OTU-bin.932",  "OTU-bin.937",  "OTU-bin.939",  "OTU-bin.241",  "OTU-bin.140", 
                               "OTU-bin.145",  "OTU-bin.146",  "OTU-bin.545",  "OTU-bin.157",  "OTU-bin.158", 
                               "OTU-bin.550",  "OTU-bin.552",  "OTU-bin.554",  "OTU-bin.952",  "OTU-bin.1009",
                               "OTU-bin.160",  "OTU-bin.165",  "OTU-bin.206",  "OTU-bin.207",  "OTU-bin.600", 
                               "OTU-bin.601",  "OTU-bin.961",  "OTU-bin.962",  "OTU-bin.965",  "OTU-bin.968", 
                               "OTU-bin.969",  "OTU-bin.175",  "OTU-bin.217",  "OTU-bin.218",  "OTU-bin.579", 
                               "OTU-bin.594",  "OTU-bin.734",  "OTU-bin.808",  "OTU-bin.612",  "OTU-bin.616", 
                               "OTU-bin.996",  "OTU-bin.363",  "OTU-bin.183",  "OTU-bin.223",  "OTU-bin.224", 
                               "OTU-bin.225",  "OTU-bin.586",  "OTU-bin.229",  "OTU-bin.654",  "OTU-bin.509", 
                               "OTU-bin.624",  "OTU-bin.984",  "OTU-bin.985",  "OTU-bin.627",  "OTU-bin.987", 
                               "OTU-bin.988",  "OTU-bin.989",  "OTU-bin.190",  "OTU-bin.590",  "OTU-bin.596", 
                               "OTU-bin.597",  "OTU-bin.819"  ))     #Control
d.name$name<-factor(d.name$name,
                    levels = c("OTU-bin.1145", "OTU-bin.969",  "OTU-bin.354",  "OTU-bin.793",  "OTU-bin.586", 
                               "OTU-bin.1035", "OTU-bin.769",  "OTU-bin.880", "OTU-bin.269",  "OTU-bin.832", 
                               "OTU-bin.537",  "OTU-bin.321",  "OTU-bin.1100", "OTU-bin.248",  "OTU-bin.1147",
                               "OTU-bin.1067", "OTU-bin.848",  "OTU-bin.368",  "OTU-bin.905",  "OTU-bin.120", 
                               "OTU-bin.996",  "OTU-bin.644",  "OTU-bin.719",  "OTU-bin.657",  "OTU-bin.1116",
                               "OTU-bin.590",  "OTU-bin.374",  "OTU-bin.808",  "OTU-bin.924",  "OTU-bin.901", 
                               "OTU-bin.390",  "OTU-bin.1024", "OTU-bin.898",  "OTU-bin.1075", "OTU-bin.579",
                               "OTU-bin.766",  "OTU-bin.853",  "OTU-bin.539",  "OTU-bin.339",  "OTU-bin.206", 
                               "OTU-bin.1129", "OTU-bin.353",  "OTU-bin.1096", "OTU-bin.965",  "OTU-bin.76",  
                               "OTU-bin.1157", "OTU-bin.698",  "OTU-bin.817",  "OTU-bin.526",  "OTU-bin.550", 
                               "OTU-bin.12",  "OTU-bin.1159", "OTU-bin.379" , "OTU-bin.19",   "OTU-bin.1166",
                               "OTU-bin.254",  "OTU-bin.125",  "OTU-bin.932",  "OTU-bin.685",  "OTU-bin.523", 
                               "OTU-bin.516",  "OTU-bin.671",  "OTU-bin.499",  "OTU-bin.641",  "OTU-bin.247", 
                               "OTU-bin.282",  "OTU-bin.667",  "OTU-bin.340",  "OTU-bin.1006", "OTU-bin.207", 
                               "OTU-bin.699",  "OTU-bin.1038", "OTU-bin.1055", "OTU-bin.407",  "OTU-bin.1014",
                               "OTU-bin.856",  "OTU-bin.17",   "OTU-bin.465",  "OTU-bin.27",   "OTU-bin.775", 
                               "OTU-bin.406",  "OTU-bin.742",  "OTU-bin.54",   "OTU-bin.40",   "OTU-bin.1162",
                               "OTU-bin.30",   "OTU-bin.470",  "OTU-bin.363",  "OTU-bin.654",  "OTU-bin.355", 
                               "OTU-bin.352",  "OTU-bin.44",   "OTU-bin.612",  "OTU-bin.819",  "OTU-bin.63",  
                               "OTU-bin.1155", "OTU-bin.532",  "OTU-bin.862",  "OTU-bin.79",   "OTU-bin.502", 
                               "OTU-bin.469",  "OTU-bin.456",  "OTU-bin.517",  "OTU-bin.515",  "OTU-bin.278", 
                               "OTU-bin.347",  "OTU-bin.704",  "OTU-bin.1022", "OTU-bin.1054", "OTU-bin.912", 
                               "OTU-bin.894",  "OTU-bin.1",    "OTU-bin.157",  "OTU-bin.1124", "OTU-bin.31",  
                               "OTU-bin.687",  "OTU-bin.1093", "OTU-bin.165",  "OTU-bin.221",  "OTU-bin.624", 
                               "OTU-bin.527",  "OTU-bin.265",  "OTU-bin.801",  "OTU-bin.1115", "OTU-bin.962", 
                               "OTU-bin.513",  "OTU-bin.507",  "OTU-bin.1049", "OTU-bin.552",  "OTU-bin.183", 
                               "OTU-bin.718",  "OTU-bin.1131", "OTU-bin.249",  "OTU-bin.889",  "OTU-bin.1135", 
                               "OTU-bin.937",  "OTU-bin.740",  "OTU-bin.1090", "OTU-bin.596",  "OTU-bin.666",
                               "OTU-bin.1065", "OTU-bin.627",  "OTU-bin.1094" ))     #-40%
d.name$name<-factor(d.name$name,
                    levels = c("OTU-bin.293",  "OTU-bin.968",  "OTU-bin.926",  "OTU-bin.590",  "OTU-bin.1017",
                               "OTU-bin.354",  "OTU-bin.386",  "OTU-bin.594",  "OTU-bin.353",  "OTU-bin.832", 
                               "OTU-bin.269",  "OTU-bin.77",   "OTU-bin.363",  "OTU-bin.261",  "OTU-bin.677", 
                               "OTU-bin.654",  "OTU-bin.790",  "OTU-bin.870",  "OTU-bin.698",  "OTU-bin.901", 
                               "OTU-bin.444",  "OTU-bin.239",  "OTU-bin.984",  "OTU-bin.1106", "OTU-bin.112", 
                               "OTU-bin.996",  "OTU-bin.793",  "OTU-bin.160",  "OTU-bin.1019", "OTU-bin.988", 
                               "OTU-bin.819",  "OTU-bin.1009", "OTU-bin.30",   "OTU-bin.682",  "OTU-bin.241", 
                               "OTU-bin.321"))   #减雨60%

library(dplyr)
d.name<-arrange(d.name,name)
V(Oneyear)$domian<-d.name$domian
V(Oneyear)$phylum<-d.name$Color
V(Oneyear)$key<-d.name$key
V(Oneyear)$abundance<-d.name$mean
V(Oneyear)$abundance     #检查合并结果
V(Oneyear)$domian        #检查合并结果
V(Oneyear)$phylum

#2.1 颜色设置 phylum
library(RColorBrewer)
#colors<-brewer.pal(11,"Paired")    #该调色板最多12种颜色
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
shapes<-c("square","circle")
Xingz<-att[key,]$domian
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
#colors.ed        #检查对应关系

#2.4 节点大小及连线宽度设置
V(Oneyear)$size <- log(V(Oneyear)$abundance*100)
E(Oneyear)$width <- abs(log10(E(Oneyear)$weight+0.0001))/2

#2.4 针对Control组修改
V(Oneyear)$size    #不可以有负数，故缩放数值
library(scales)
V(Oneyear)$size<- rescale(V(Oneyear)$size, to = c(1, 5), from = c(-2.303, 8.204))


#3. 作图
circle<- layout_in_circle(Oneyear)   
plot(Oneyear, layout=circle, vertex.label="",edge.curved=.2)   #单独画一个图
plot(Oneyear, layout=circle, vertex.label="")   #单独画一个图
dev.off()

#4 补充作图
sphere<- layout_on_sphere(Oneyear)
circle<- layout_in_circle(Oneyear)
fr<- layout_with_fr(Oneyear)     #?
kk<- layout_with_kk(Oneyear)     #?
lgl<-layout_with_lgl(Oneyear)
gem<-layout_with_gem(Oneyear)
plot(Oneyear, layout=gem, vertex.label="",edge.curved=.2)

layouts <- grep("^layout_", ls("package:igraph"), value=TRUE)[-1] 
# Remove layouts that do not apply to our graph.
layouts <- layouts[!grepl("bipartite|merge|norm|sugiyama|tree", layouts)]
layouts<-layouts[c(1,3:15)]
par(mfrow=c(3,3), mar=c(1,1,1,1))

for (layout in layouts) {
  print(layout)
  l <- do.call(layout, list(Oneyear)) 
  plot(Oneyear, edge.arrow.mode=0, layout=l, main=layout,vertex.label="") }
dev.off()


#做一个bin的散点图，整一个图例
att$Color<-factor(att$Color,
                    levels = c("Acidobacteriota","Actinobacteriota","Armatimonadota","Bacteroidota","Bdellovibrionota",
                               "Chloroflexota", "Cyanobacteria","Desulfobacterota","Elusimicrobiota","Fibrobacterota",
                               "Firmicutes","Gemmatimonadota", "Halobacteriota","Margulisbacteria","Methanobacteriota",
                               "Methylomirabilota","Myxococcota","Nitrospirota","Patescibacteria","Planctomycetota",
                               "Proteobacteria","Thermoproteota" ))  

library(ggplot2)
ggplot(att,aes(x = abundance,y=degree))+
  geom_point(aes(x = abundance,y=degree,fill=Color),shape=21)+
  scale_fill_manual(values = c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#AB7519",
                               "#CAB2D6","#6A3D9A","#FFFF99","#B15928","#E64E00","#E6EB00","#65B48E","#3E5CC5","#FB1299",
                               "#F3C846","#D882AD","#DDDEDE","#AAAB20"))+
  theme_classic()+
  theme(axis.text = element_text(colour = "black"))






#2024.3.20 网络节点和链接可视化——百分比条形图
net.perc<-read.table("net.perc.txt",header=T,sep="\t",quote="")
library(reshape2)
note<-net.perc[,c(1:3)]
note<-melt(note,id.vars= c("Group"),variable.name = "mic",value.name = "percent")
line<-net.perc[,c(1,4,5)]
line<-melt(line,id.vars= c("Group"),variable.name = "line",value.name = "percent")
#排序
note$mic<-factor(note$mic,levels = c("Archaea","Bacteria"))


library(ggplot2)
note$Group<-factor(note$Group,levels = c("-60%","-40%","Control","40%","60%"))
ggplot(note,aes(fill=mic,x=Group,y=percent))+
  geom_bar(position = "stack",stat = "identity",width = 0.3)+
  theme_classic()+
  theme(axis.text = element_text(colour = "black"))+
  ylab("node percentage (%)")+labs(fill="domian")+
  scale_fill_manual(values=c(Archaea = "#FB9A99",Bacteria = "#3CCEE3"),
                    labels=c("Archaea","Bacteria"))

line$Group<-factor(line$Group,levels = c("-60%","-40%","Control","40%","60%"))
ggplot(line,aes(fill=line,x=Group,y=percent))+
  geom_bar(position = "stack",stat = "identity",width = 0.3)+
  theme_classic()+
  theme(axis.text = element_text(colour = "black"))+
  ylab("line percentage (%)")+labs(fill="interaction")+
  scale_fill_manual(values=c(mutualExclusion = "#E25A1C",copresence = "#33A02C"),
                    labels=c("mutualExclusion","copresence"))



###对每个处理中phylum进行计数
att<-read.table("att_60.txt",header=T,sep="\t",quote="")
phylum_60 <- aggregate(att$Color, by = list(att$Color), FUN = length)   #计数公式
colnames(phylum_60)<-c("phylum","phylum_60")

att_40<-read.table("att_40.txt",header=T,sep="\t",quote="")
phylum_40 <- aggregate(att_40$Color, by = list(att_40$Color), FUN = length)
colnames(phylum_40)<-c("phylum","phylum_40")

att40<-read.table("att40.txt",header=T,sep="\t",quote="")
phylum40 <- aggregate(att40$Color, by = list(att40$Color), FUN = length)
colnames(phylum40)<-c("phylum","phylum40")

att60<-read.table("att60.txt",header=T,sep="\t",quote="")
phylum60 <- aggregate(att60$Color, by = list(att60$Color), FUN = length)
colnames(phylum60)<-c("phylum","phylum60")

att_Con<-read.table("att_Con.txt",header=T,sep="\t",quote="")
phylum_Con <- aggregate(att_Con$Color, by = list(att_Con$Color), FUN = length)
colnames(phylum_Con)<-c("phylum","phylum_Con")

library(dplyr)
df1 <- full_join(phylum_60, phylum_40,by = "phylum")    #合并
df1 <- full_join(df1, phylum_Con,by = "phylum")
df1 <- full_join(df1, phylum40,by = "phylum")
df1 <- full_join(df1, phylum60,by = "phylum")

df1[is.na(df1)]=0
df1$sum <- rowSums(df1[, 2:6])
df1<-arrange(df1,sum)
#计数完成，见df1

df<-df1[,c(1:6)]
df<-melt(df,id.vars= c("phylum"),variable.name = "Group",value.name = "count")
#条形图
df$phylum<-factor(df$phylum,
                  levels = c("Acidobacteriota","Actinobacteriota","Armatimonadota","Bacteroidota","Bdellovibrionota",
                             "Chloroflexota", "Cyanobacteria","Desulfobacterota","Elusimicrobiota","Fibrobacterota",
                             "Firmicutes","Gemmatimonadota", "Halobacteriota","Margulisbacteria","Methanobacteriota",
                             "Methylomirabilota","Myxococcota","Nitrospirota","Patescibacteria","Planctomycetota",
                             "Proteobacteria","Thermoproteota" ))  

ggplot(df, aes(x = Group, y = count,fill=phylum)) + 
  geom_bar(position="stack", stat="identity",width = 0.3)+
  scale_fill_manual(values = c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#AB7519",
                               "#CAB2D6","#6A3D9A","#FFFF99","#B15928","#E64E00","#E6EB00","#65B48E","#3E5CC5","#FB1299",
                               "#F3C846","#D882AD","#DDDEDE","#AAAB20"))+
  theme_classic()+ theme(axis.text = element_text(colour = "black"))+
  labs( y = "count")

#百分比条形图
df1<-mutate(df1,p_60=phylum_60/36*100)
df1<-mutate(df1,p_40=phylum_40/143*100)
df1<-mutate(df1,p_Con=phylum_Con/292*100)
df1<-mutate(df1,p40=phylum40/157*100)
df1<-mutate(df1,p60=phylum60/58*100)
df2<-df1[,c(1,8:12)]

#条形图
df2<-melt(df2,id.vars= c("phylum"),variable.name = "Group",value.name = "percentage")
df2$phylum<-factor(df2$phylum,
                  levels = c("Acidobacteriota","Actinobacteriota","Armatimonadota","Bacteroidota","Bdellovibrionota",
                             "Chloroflexota", "Cyanobacteria","Desulfobacterota","Elusimicrobiota","Fibrobacterota",
                             "Firmicutes","Gemmatimonadota", "Halobacteriota","Margulisbacteria","Methanobacteriota",
                             "Methylomirabilota","Myxococcota","Nitrospirota","Patescibacteria","Planctomycetota",
                             "Proteobacteria","Thermoproteota" ))  

ggplot(df2, aes(x = Group, y = percentage,fill=phylum)) + 
  geom_bar(position="stack", stat="identity",width = 0.3)+
  scale_fill_manual(values = c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#AB7519",
                               "#CAB2D6","#6A3D9A","#FFFF99","#B15928","#E64E00","#E6EB00","#65B48E","#3E5CC5","#FB1299",
                               "#F3C846","#D882AD","#DDDEDE","#AAAB20"))+
  theme_classic()+ theme(axis.text = element_text(colour = "black"))+
  labs( y = "percentage (%)")

#分细菌 古菌
df3<-df1[-c(5,9,13),c(1,8:12)]
df4<-df1[c(5,9,13),c(1,8:12)]

df4<-melt(df4,id.vars= c("phylum"),variable.name = "Group",value.name = "percentage")
ggplot(df3, aes(x = Group, y = percentage,fill=phylum)) + 
  geom_bar(position="stack", stat="identity",width = 0.3)+
  scale_fill_manual(values = c(Acidobacteriota="#A6CEE3",Actinobacteriota="#1F78B4",Armatimonadota="#B2DF8A",
                               Bacteroidota="#33A02C",Bdellovibrionota="#FB9A99",Chloroflexota="#E31A1C",Cyanobacteria="#FDBF6F",
                               Desulfobacterota="#FF7F00",Elusimicrobiota="#AB7519",Fibrobacterota="#CAB2D6",
                               Firmicutes="#6A3D9A",Gemmatimonadota="#FFFF99",Margulisbacteria="#E64E00",
                               Methylomirabilota="#65B48E",Myxococcota="#3E5CC5",Nitrospirota="#FB1299",
                               Patescibacteria="#F3C846",Planctomycetota="#D882AD",Proteobacteria="#DDDEDE"))+
  theme_classic()+ theme(axis.text = element_text(colour = "black"))+
  labs( y = "percentage (%)")
ggplot(df4, aes(x = Group, y = percentage,fill=phylum)) + 
  geom_bar(position="stack", stat="identity",width = 0.3)+
  scale_fill_manual(values = c(Methanobacteriota="#E6EB00",Thermoproteota="#AAAB20",Halobacteriota="#B15928"))+
  theme_classic()+ theme(axis.text = element_text(colour = "black"))+
  scale_y_continuous(limits = c(0,10))+
  labs(y = "percentage (%)")
"#6A3D11"
#参考
scale_fill_manual(values = c(Acidobacteriota="#A6CEE3",Actinobacteriota="#1F78B4",Armatimonadota="#B2DF8A",
                             Bacteroidota="#33A02C",Bdellovibrionota="#FB9A99",Chloroflexota="#E31A1C",Cyanobacteria="#FDBF6F",
                             Desulfobacterota="#FF7F00",Elusimicrobiota="#AB7519",Fibrobacterota="#CAB2D6",
                             Firmicutes="#6A3D9A",Gemmatimonadota="#FFFF99",Halobacteriota="#B15928",Margulisbacteria="#E64E00",
                             Methanobacteriota="#E6EB00",Methylomirabilota="#65B48E",Myxococcota="#3E5CC5",Nitrospirota="#FB1299",
                             Patescibacteria="#F3C846",Planctomycetota="#D882AD",Proteobacteria="#DDDEDE",Thermoproteota="#AAAB20"))



#2024.3.26  子网络作图
graph<-Oneyear
dev.off()

library(igraph)
subgraph <- induced_subgraph(graph, "OTU-bin.17")

# 绘制子网络图
plot(subgraph, edge.label=E(Oneyear)$weight)





#---------------------------------------------------
#2024.11.5 (一) 数据前处理
bin<-read.table("abund.txt",header = T,sep="\t") 
row.names(bin)<-bin$ID

bin<-bin[,3:118]
min_nonzero_value <- min(bin[bin != 0], na.rm = TRUE)    #查找非零最小值
# 0.0064055          0.00001
bin1 <- bin+0.00001

#分生态系统
bin1.t<-as.data.frame(t(bin1))
bin1.t$Category<-row.names(bin1.t)
library(tidyr)
bin1.t<- separate(bin1.t,col="Category",into=c("Category","id"),sep = "_",
                remove = F, convert = FALSE)

bin.MF.t<-bin1.t[bin1.t$Category=="MF",1:1954]
bin.MF<-as.data.frame(t(bin.MF.t))
bin.SM.t<-bin1.t[bin1.t$Category=="SM",1:1954]
bin.SM<-as.data.frame(t(bin.SM.t))
bin.SG.t<-bin1.t[bin1.t$Category=="SG",1:1954]
bin.SG<-as.data.frame(t(bin.SG.t))

#所有样本加和，选取每个生态系统丰度最高的前300个，构建网络图
library(dplyr)
bin.MF<-mutate(bin.MF,sum=rowSums(bin.MF))
bin.MF_desc <- bin.MF %>%
  arrange(desc(sum))
bin.MF<-bin.MF_desc[-c(301:1954),-35]

bin.SM<-mutate(bin.SM,sum=rowSums(bin.SM))
bin.SM_desc <- bin.SM %>%
  arrange(desc(sum))
bin.SM<-bin.SM_desc[-c(301:1954),-58]

bin.SG<-mutate(bin.SG,sum=rowSums(bin.SG))
bin.SG_desc <- bin.SG %>%
  arrange(desc(sum))
bin.SG<-bin.SG_desc[-c(301:1954),-26]

write.table(bin.MF,"cytoscape.binMF.txt",sep="\t",row.names=T,quote=F)
write.table(bin.SM,"cytoscape.binSM.txt",sep="\t",row.names=T,quote=F)
write.table(bin.SG,"cytoscape.binSG.txt",sep="\t",row.names=T,quote=F)
#导出txt文件，导入Cytoscape中 进行网络模型构建



#（二）MF 2024.11.7
#1.Cytoscape运行结果数据导入与整理——核心数据框Oneyear与att
install.packages("igraph")
install.packages("RColorBrewer")

library(igraph)
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



###SM------------------------------------------
Oneyear<-read_graph("SM.graphml",format = "graphml")
Oneyear<-as.undirected(Oneyear,mode = "each")
plot(Oneyear)

bin.cyto<-read.table("SM.node.txt",header=T,sep="\t",quote="")
anno<-read.table("anno.txt",header = T,sep="\t") 
node.fam<-merge(bin.cyto,anno,by="ID")
#整理数据，给node.fam添加key(ID)  Color (Phylum)两列 
node.fam$Color<-node.fam$Phylum
node.fam$key<-node.fam$ID
#增加bin丰度均值信息到node.fam数据框
abund<-bin.SM
abund$mean<-rowMeans(abund)        
abund[,59]<-row.names(abund) 
colnames(abund)[59]<-"ID"
node.fam<-merge(node.fam,abund,by="ID")
node.fam<-node.fam[,c(1:23,81)]

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
unique(att$Phylum) #检查门的数量，MF共12个,,SM 14个
library(RColorBrewer)
colors<-brewer.pal(12,"Paired")    #该调色板最多12种颜色

colors.point<-c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C",
                "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928",
                "#E6EB00", "#65B48E", "#AAAB20", "#D882AD") 
point.type<-c("Acidobacteriota",   "Actinomycetota",    "Bacillota" ,        "Bacteroidota",  "Bdellovibrionota","Campylobacterota",    
              "Chloroflexota" ,    "Desulfobacterota",  "Gemmatimonadota", "Patescibacteria", "Pseudomonadota", "Verrucomicrobiota",  
              "Methylomirabilota", "Myxococcota",   "Planctomycetota",  "Thermoproteota")
point.type<-factor(point.type,
                      levels = c("Acidobacteriota",   "Actinomycetota",    "Bacillota" ,        "Bacteroidota",  "Bdellovibrionota","Campylobacterota",    
                                "Chloroflexota" ,    "Desulfobacterota",  "Gemmatimonadota", "Patescibacteria", "Pseudomonadota", "Verrucomicrobiota",  
                                "Methylomirabilota", "Myxococcota",   "Planctomycetota",  "Thermoproteota"))
names(colors.point)<-levels(point.type)
colors.point        #检查对应关系   色板设置，仅运行一次即可
V(Oneyear)$color <- colors.point[V(Oneyear)$phylum]   
V(Oneyear)$color      


#2.2形状 domian
shapes<-c("square","circle")
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
#8*8出图



###SG------------------------------------------
Oneyear<-read_graph("SG.graphml",format = "graphml")
Oneyear<-as.undirected(Oneyear,mode = "each")
plot(Oneyear)

bin.cyto<-read.table("SG.node.txt",header=T,sep="\t",quote="")
anno<-read.table("anno.txt",header = T,sep="\t") 
node.fam<-merge(bin.cyto,anno,by="ID")
#整理数据，给node.fam添加key(ID)  Color (Phylum)两列 
node.fam$Color<-node.fam$Phylum
node.fam$key<-node.fam$ID
#增加bin丰度均值信息到node.fam数据框
abund<-bin.SG
abund$mean<-rowMeans(abund)        
abund[,27]<-row.names(abund) 
colnames(abund)[27]<-"ID"
node.fam<-merge(node.fam,abund,by="ID")
node.fam<-node.fam[,c(1:23,49)]

att<-node.fam
row.names(att)<-att$key 
#整理数据索引列key
ver<-data.frame(V(Oneyear))
ver$name<-row.names(ver)
ver<-separate(ver,name,into=c("o","ID"),sep=4)
ver<-ver[,c(1,3)]
key<-ver$ID   

#2. 将domian  phylum  key;  abund 添加到Oneyear，只能一列一列加(注意行名对应起来)
###特别注意：att需要索引的全部列必须在索引key创建前整理好，否则后来添加的列无效；
V(Oneyear)$Yanse<-att[key,]$Color
V(Oneyear)$domian<-att[key,]$Domain
V(Oneyear)$phylum<-att[key,]$Color
V(Oneyear)$abundance<-att[key,]$mean
V(Oneyear)$key<-att[key,]$key

#2.1 颜色设置 phylum
unique(att$Phylum) #检查门的数量，MF共12个,,SM 14个;  SG 12个
library(RColorBrewer)
colors<-brewer.pal(12,"Paired")    #该调色板最多12种颜色

colors.point<-c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C",
                "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928",
                "#E6EB00", "#65B48E", "#AAAB20", "#D882AD") 
point.type<-c("Acidobacteriota",   "Actinomycetota",    "Bacillota" ,        "Bacteroidota",  "Bdellovibrionota","Campylobacterota",    
              "Chloroflexota" ,    "Desulfobacterota",  "Gemmatimonadota", "Patescibacteria", "Pseudomonadota", "Verrucomicrobiota",  
              "Methylomirabilota", "Myxococcota",   "Planctomycetota",  "Thermoproteota")
point.type<-factor(point.type,
                   levels = c("Acidobacteriota",   "Actinomycetota",    "Bacillota" ,        "Bacteroidota",  "Bdellovibrionota","Campylobacterota",    
                              "Chloroflexota" ,    "Desulfobacterota",  "Gemmatimonadota", "Patescibacteria", "Pseudomonadota", "Verrucomicrobiota",  
                              "Methylomirabilota", "Myxococcota",   "Planctomycetota",  "Thermoproteota"))
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


#---------------做图例
#直接使用上面的画板
library(ggplot2)
att.1<-att
att.1$size <- log(att.1$mean*100)
ggplot(att,aes(x = abundance,y=degree))+
  geom_point(aes(x = abundance,y=degree,fill=Color),shape=21)+
  scale_fill_manual(values = colors.point)+
  theme_classic()+
  theme(axis.text = element_text(colour = "black"))

Acidobacteriota    Actinomycetota         Bacillota      Bacteroidota  Bdellovibrionota 
"#A6CEE3"         "#1F78B4"         "#B2DF8A"         "#33A02C"         "#FB9A99" 
Campylobacterota     Chloroflexota  Desulfobacterota   Gemmatimonadota   Patescibacteria 
"#E31A1C"         "#FDBF6F"         "#FF7F00"         "#CAB2D6"         "#6A3D9A" 
Pseudomonadota Verrucomicrobiota Methylomirabilota       Myxococcota   Planctomycetota 
"#FFFF99"         "#B15928"         "#E6EB00"         "#65B48E"         "#AAAB20" 
Thermoproteota 
"#D882AD" 













#---------------------#网络模块化2024.11.25
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
if (min_weight <= 0) {
  E(Oneyear)$weight[E(Oneyear)$weight <= 0] <- 2.2e-16
}

member<-edge.betweenness.community(Oneyear,weight=E(Oneyear)$weight,directed=F)   #directed=T时就代表有向线。
print(member)  #groups: 19, mod: 0.61  模块度modularity
#警告：在执行基于边介数中心性的社区划分时，算法计算了多个可能的社区划分，
#并从中选择了模块度（modularity）最高的一个作为最终结果。
#模块度是衡量社区划分质量的一个指标，它反映了社区内部节点连接的紧密程度相对于社区间节点连接的紧密程度。
#membership <- cluster_edge_betweenness(Oneyear, directed = FALSE, weights = E(Oneyear)$weight)
#cluster_fast_greedy 马尔科夫
#模块化指标Q——modularity
membership<- membership(member)
modularity <- modularity(Oneyear, membership = membership)   #0.6382

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

g1 <- graph.data.frame(Oneyear)
att2$<-as.data.frame(V(Oneyear)$)
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