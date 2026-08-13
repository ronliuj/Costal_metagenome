library(igraph)
library(ggClusterNet)
library(ggplot2)
library(ggraph)


bin.MF<-read.table("bin.MF.txt",header=T,sep="\t",quote="")
anno<-read.table("anno.txt",header = T,sep="\t") 

cor<-cor(t(bin.MF),method="spearman")
cor[abs(cor)<0.5]<-0

result2 = model_maptree2(cor = cor,
                         method = "cluster_fast_greedy"
)
node = result2[[1]]
aa<- result2[[2]]
bb<- result2[[3]]

nodes2 = node %>% inner_join(aa,by = c("elements" = "ID"))
edge = edgeBuild(cor = cor,node = node)

pnet <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor)),
                                data = edge, size = 0.5) +
  geom_point(aes(X1, X2,fill = group,size = 2),pch = 21, data = nodes2) +
  scale_colour_brewer(palette = "Set1") +
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



result2 <- model_igraph(cor = cor,
                        method = "cluster_fast_greedy",
                        seed = 12
)
node = result2[[1]]
head(node)

dat = result2[[2]]
head(dat)
tem = data.frame(mod = dat$model,col = dat$color) %>%  
  dplyr::distinct( mod, .keep_all = TRUE)  
col = tem$col
names(col) = tem$mod


edge = edgeBuild(cor = cor,node = node)
colnames(edge)[8] = "cor"
head(edge)

tem2 = dat %>% 
  dplyr::select(OTU,model,color) %>%
  dplyr::right_join(edge,by =c("OTU" = "OTU_1" ) ) %>%
  dplyr::rename(OTU_1 = OTU,model1 = model,color1 = color)
head(tem2)

tem3 = dat %>% 
  dplyr::select(OTU,model,color) %>%
  dplyr::right_join(edge,by =c("OTU" = "OTU_2" ) ) %>%
  dplyr::rename(OTU_2 = OTU,model2 = model,color2 = color)
head(tem3)

tem4 = tem2 %>%inner_join(tem3)
head(tem4)

edge2 = tem4 %>% mutate(color = ifelse(model1 == model2,as.character(model1),"across"),
                        manual = ifelse(model1 == model2,as.character(color1),"#C1C1C1")
)
head(edge2)
col_edge = edge2 %>% dplyr::distinct(color, .keep_all = TRUE)  %>% 
  dplyr::select(color,manual)
col0 = col_edge$manual
names(col0) = col_edge$color

library(ggnewscale)

p1 <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = color),
                              data = edge2, size = 1) +
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




















