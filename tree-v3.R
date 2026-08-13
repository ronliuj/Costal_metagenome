library(ggtreeExtra) # 设置叠加的包
library(ggstar) # 提供几何图形
library(ggplot2) # 
library(ggtree) # 绘制进化树
library(treeio)
library(ggnewscale) # 创建新的scale，多个fill或者color
library(TDbook)
library(RColorBrewer)
library(reshape2)
library(picante)
library(dplyr)

arc.tree.raw <- phytools::read.newick("blue.ar53.user_msa.fasta.gz.treefile")
bac.tree.raw <- phytools::read.newick("blue.bac120.user_msa.fasta.gz.treefile")

tree.raw <- ape::bind.tree(bac.tree.raw, arc.tree.raw)

bins<-read.table("bins-anno.txt",header = T,sep="\t")

row.names(bins)<-bins$genome

##remove unmattched
phylo = match.phylo.data(tree.raw, bins)

tree = phylo$phy

treephy<-as.phylo(tree)

x <- as_tibble(treephy)

tree2<-as.treedata(x)

meta<-read.table("para-0528.txt",header = T,sep="\t")
row.names(meta)<-meta$ID

abund<-bins[,3:118]


prg<-data.frame(matrix(nrow=13671,ncol=3))
colnames(prg)<-c("bin","Cat","Pgr")
prg$bin<-rep(x$label[1:1953],7)
prg$Cat<-rep(c("MF","SM","SG","HHS","BHS","SCS","ECS"),each=1953)
prg$Pgr<-"DL"


set.seed(123)
prg <- prg %>%
  mutate(
    Pgr = if (
      is.factor(Pgr) && !"HoS" %in% levels(Pgr)
    ) {
      factor(v, levels = c(levels(Pgr), "HoS"))
    } else {
      Pgr
    },
    Pgr = replace(
      Pgr,
      sample(n(), max(1, round(0.1 * n()))),
      "HoS"
    )
  )


bin.cha<-bins[,122:127]

aa<-data.frame(table(bin.cha$Phylum))
aa<-arrange(aa,desc(aa$Freq))

bin.cha$Phy2<-ifelse(is.element(bin.cha$Phylum,aa[1:15,]$Var1),bin.cha$Phylum,"Other")
##只有genmome个数大于15个的才进行着色，其他的均为Others

unique(bin.cha$Phy2)



colors<-brewer.pal(8,"Dark2")
co2<-brewer.pal(9,"Set1")
co3<-brewer.pal(12,"Paired")
colors<-c(colors,co2,co3)
colors<-unique(colors)[1:16]
colors<-c(colors[c(1:7,9:16)],"lightgrey")
names(colors)<-c(unique(bin.cha$Phy2)[c(1:8,10:16)],"Other")

Pseudomonadota      Bacteroidota     Chloroflexota    Actinomycetota  Campylobacterota 
"#1B9E77"         "#D95F02"         "#7570B3"         "#E7298A"         "#66A61E" 
Myxococcota         Bacillota   Gemmatimonadota   Planctomycetota    Thermoproteota 
"#E6AB02"         "#A6761D"         "#E41A1C"         "#377EB8"         "#4DAF4A" 
Desulfobacterota   Patescibacteria   Acidobacteriota  Bdellovibrionota Verrucomicrobiota 
"#984EA3"         "#FF7F00"         "#FFFF33"         "#A65628"         "#F781BF" 
Other 
"lightgrey" 

#write.table(colors,"phylum-color.code",sep="\t",quote=T)

x<-data.frame(x)[1:1953,]
row.names(x)<-x$label

bin.cha$node<-x[row.names(bin.cha),]$node


# The circular layout tree.

bin.cha$ID<-row.names(bin.cha)


ggtree(tree2, layout="fan", size=0.15, open.angle=5) +
  geom_tiplab()+ #添加tiplabe
  geom_hilight(data=bin.cha,mapping=aes(node=node,fill=Phy2,color=Phy2),
               extendto=1.5, alpha=1, 
               size=0.05) + xlim(-3, NA)+
  geom_text(aes(label=node))+ #添加node号码，以便整理着色表格
  scale_color_manual(values = colors)+
  scale_fill_manual(values = colors)#添加门水平注释

#geom_hilight是一个node,一个node进行着色的，因此如果想要一个phylum一起着色，应当适用phylum最顶上的node
#根据上图整理phylum着色表
node.ph<-read.table("Node-phy.txt",header = T,sep="\t")

ggtree(tree2, layout="fan", size=0.15,
          branch.length = "none", open.angle=180, color="grey",alpha=0.3) + ##此处设置线的颜色
  #geom_tiplab()+ #添加tiplabe
  #geom_text(aes(label=node))+ #添加node号码，以便整理着色表格
  geom_hilight(data=node.ph,mapping=aes(node=node,fill=Phy2,color=Phy2),
               extendto=-Inf, alpha=1, 
               size=0.3,to.bottom = T) + xlim(-3, NA)+
  scale_color_manual(values = colors)+
  scale_fill_manual(values = colors)#添加门水平注释

##branch.length="none" 端点平齐，内向延长；不选择本选项，extendto可以增大使得着色框平齐

group1 <- split(row.names(bin.cha),bin.cha$Phy2)
tree2<-groupOTU(tree2,group1)

p <- ggtree(tree2, layout="fan", size=0.25, open.angle=10)
p2<-p + 
  geom_fruit(data=bin.cha, geom=geom_tile,
             mapping=aes(y=ID, fill=Phy2), width = 0.2,
             offset = 0.02) + scale_fill_manual(values = colors) + new_scale_fill()+
  geom_fruit(
    data=prg,
    geom=geom_tile,
    mapping=aes(y=bin, x=Cat, fill=Pgr),
    offset=0.08,   # The distance between external layers, default is 0.03 times of x range of tree.
    pwidth=0.25 # width of the external layer, default is 0.2 times of x range of tree.
  ) + 
  scale_fill_manual(
    values=c("#4DBBD5FF", "#F39B7FFF"),
    guide=guide_legend(keywidth=0.5, keyheight=0.5, order=3)
  ) 





selected<-read.table("固碳/selected.fix.genome",header = T,sep="\t")

setdiff(selected$ID,bins$genome)

colnames(selected)[1]<-"genome"

row.names(selected)<-selected$genome

selected$ID<-bins[row.names(selected),]$ID
selected$Phylum<-bins[row.names(selected),]$Phylum

library(ggtree)
library(ape)


group<-bin.cha[,c("ID","Phy2")]

group$group<-"A"

group[row.names(selected),]$group<-selected$Phylum

groupInfo <- split(row.names(group), group$group)


tree2 <- groupOTU(tree2, groupInfo)

colors<-c("#FF7F00","#6A3D9A","lightgrey")
names(colors)<-c("Desulfobacterota","Chloroflexota","A")


group1<-selected[selected$Phylum=="Desulfobacterota",]$genome

group1.id<-selected[selected$Phylum=="Desulfobacterota",]$ID

group2<-selected[selected$Phylum=="Chloroflexota",]$genome

group2.id<-selected[selected$Phylum=="Chloroflexota",]$ID


label_data1 <- data.frame(
  label = group1,
  new_label = group1.id
)

label_data2 <- data.frame(
  label = group2,
  new_label = group2.id
)


p1<-ggtree(tree2, layout="circular",aes(color=group),size=0.5) +
  geom_tiplab2(aes(color=group))+
  geom_tiplab2(aes(color=group),size=2)+
  scale_color_manual(values = colors)+
  theme(legend.position = "right")


p2<-ggtree(tree2, layout="circular",aes(color=group),size=0.5) +
  #geom_tiplab2(aes(color=group))+
  #geom_tiplab2(aes(color=group),size=2)+
  scale_color_manual(values = colors)+
  theme(legend.position = "right")

p3<-ggpubr::ggarrange(p1,p2,ncol=1)

ggsave("hight.fixgenome.pdf",p3,device = "pdf",width = 40,height = 30)

ggsave("hight.fixgenome2.pdf",p2,device = "pdf",width = 40,height =20)

save.image("tree_progress.Rdata")
