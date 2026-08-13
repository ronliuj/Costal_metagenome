##生态位
library(spaa)
abund<-read.table("abund.txt",header = T,sep="\t")
row.names(abund)<-abund$ID

comm.tab<-data.frame(t(abund[,3:118]))
comm.tab<-ceiling(comm.tab*10^4)##注意这里不能太大，会超过32位int限制

##permatswap不接受强制转整数，改为presence absence数据

#Lev.index2=as.data.frame(niche.width(comm.tab, method ="levins"))

#aa<-rbind(Lev.index,Lev.index2)#均乘以一个整数不影响结果

##adapt from https://github.com/GuillemSalazar/EcolUtils/blob/master/R/EcolUtils_functions.R sepc.gen
#n=999
library(spaa)
library(vegan)
levin.index.real<-as.numeric(niche.width(comm.tab,method="levins"))
names(levin.index.real)<-colnames(comm.tab)
#levin.index.simul<-matrix(NA,ncol=dim(comm.tab)[2],nrow=999)

library(foreach)
library(doParallel)
detectCores()#40
registerDoParallel(cores = 30)#使用服务器才可以

#ps -eLf -u liujian | grep "liujian"
for (n in 1:10) {
  
simindex<-foreach(i=1:100,.combine = rbind) %dopar% {
  # 置换指定变量
  library(spaa)
  library(vegan)
  index<-as.numeric(niche.width(permatswap(comm.tab,"quasiswap",times=1)$perm,
                              method="levins"))
  
}
  write.table(simindex,paste("sim/indexsim_",n,sep=""),row.names = F,sep="\t",quote=F)
}
##以上服务器运行nohup Rscript levin-generalist.R

##以下服务器运行nohup Rscript levin-generalist2.R

abund<-read.table("abund.txt",header = T,sep="\t")
row.names(abund)<-abund$ID

comm.tab<-data.frame(t(abund[,3:118]))

comm.tab<-ceiling(comm.tab*10^4)##注意这里不能太大，会超过32位int限制

library(spaa)
library(vegan)

levin.index.real<-as.numeric(niche.width(comm.tab,method="levins"))
names(levin.index.real)<-colnames(comm.tab)

simfiles <- list.files("sim", full.names = TRUE)

levin.index.simul<-read.table(simfiles[1],header = T,sep="\t")

for(i in 2:10){
  aa<-read.table(simfiles[i],header = T,sep="\t")
  levin.index.simul<-rbind(levin.index.simul,aa)
}


colnames(levin.index.simul)<-colnames(comm.tab)
levin.index.simul<-as.data.frame(levin.index.simul)

write.table(levin.index.simul,"simul.lev.txt",row.names = T,sep="\t",quote = F)

media<-apply(levin.index.simul,2,mean)
ci<-apply(levin.index.simul,2,quantile,probs=c(0.025,0.975))
resultats<-data.frame(observed=levin.index.real,mean.simulated=media,lowCI=ci[1,],uppCI=ci[2,],sign=NA)

for (j in 1:dim(resultats)[1]){
  if (resultats$observed[j]>resultats$uppCI[j]) resultats$sign[j]<-"GENERALIST"
  if (resultats$observed[j]<resultats$lowCI[j]) resultats$sign[j]<-"SPECIALIST"
  if (resultats$observed[j]>=resultats$lowCI[j] & resultats$observed[j]<=resultats$uppCI[j]) resultats$sign[j]<-"NON SIGNIFICANT"
}

resultats$sign<-as.factor(resultats$sign)
write.table(resultats,"levin-stats.txt",row.names = T,sep="\t",quote=F)


## on pc own
abund<-read.table("abund.txt",header = T,sep="\t")
row.names(abund)<-abund$ID

comm.tab<-data.frame(t(abund[,3:118]))
comm.tab<-ceiling(comm.tab*10^4)##注意这里不能太大，会超过32位int限制

library(spaa)
library(vegan)
levin.index.real<-as.numeric(niche.width(comm.tab,method="levins"))
names(levin.index.real)<-colnames(comm.tab)

comm.tab<-comm.tab>0

Prevalence<-colSums(comm.tab)/116



resultats<-read.table("levin-stats.txt",header = T,sep="\t")

resultats$Prevalence<-Prevalence

resultats$abund<-log10(rowMeans(abund[,3:118]))

resultats$sign[resultats$sign=="NON SIGNIFICANT"]<-"NON_SIGNIFICANT"

write.table(resultats,"levin-stats2.txt",sep="\t",row.names=F,quote = F)

library(RColorBrewer)
brewer.pal(11,"Set2")

[1] "#66C2A5" "#FC8D62" "#8DA0CB" "#E78AC3" "#A6D854" "#FFD92F" "#E5C494" "#B3B3B3"

library(ggplot2)
p1<-ggplot(data=resultats[1:1953,])+
  geom_point(aes(x=Prevalence,y=abund,color=sign),size=2.5)+
  annotate("text",x=0.5,y=-3,
           label=paste0("GENERALIST = 12","\n",
                        "SPECIALIST = 1737","\n",
                        "NON_SIGNIFICANT = 204","\n"))+
  scale_color_manual(values=c(SPECIALIST = "#FC8D62", GENERALIST = "#66C2A5", NON_SIGNIFICANT= "lightgrey" ))+
  theme_bw()+ylab("log10(abundance)")+xlab("Prevalence")+ylim(c(-4.3,2))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))

genome.stat<-read.table("statgenome.txt",header = T,sep="\t")
row.names(genome.stat)<-genome.stat$file

resultats$genome<-abund[row.names(resultats),]$Genomic.bins

resultats$size<-genome.stat[resultats$genome,]$sum_len/1000000

resultats2<-resultats[1:1953,]

resultats2<-resultats2[resultats2$sign=="SPECIALIST"|resultats2$sign=="GENERALIST",]

library(ggbeeswarm)
p2<-ggplot(data=resultats2)+
  geom_boxplot(aes(x=sign,y=size,fill=sign),width=0.2,alpha=0.2,
               outlier.shape = NA)+
  geom_violin(aes(x=sign,y=size,fill=sign),width=0.5,alpha=0.2)+
  geom_quasirandom(aes(x=sign,y=size,color=sign),width=0.22,
                   size=0.2,alpha=0.5)+
  scale_color_manual(values=c(SPECIALIST = "#FC8D62", GENERALIST = "#66C2A5" ))+
  scale_fill_manual(values=c(SPECIALIST = "#FC8D62", GENERALIST = "#66C2A5" ))+
  theme_bw()+
  ylab("Genome size (Mb)")+xlab(NULL)+
  #ggtitle("CRAM-abund-permol-box")+
  annotate("text",x=1,y=8,label="Ge=3.26, Sp=3.19, p = 0.793")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black"))

shapiro.test(resultats2$size)

wilcox.test(size~sign,resultats2)
#W = 11303, p-value = 0.6135
t.test(size~sign,resultats2)
#t = 0.26855, df = 11.306, p-value = 0.7931

p2
#export 5x7
 
p1
#export 5x8
save.image("Levin.Genral.Rdata")

anno<-read.table("anno.txt",header = T,sep="\t")
row.names(anno)<-anno$genome

anno$size<-genome.stat[row.names(anno),]$sum_len


t.test(anno[anno$Familiy=="f__Nannocystaceae",]$size,anno[anno$Familiy!="f__Nannocystaceae",]$size)

sel.qg<-read.table("magsname.txt",sep="\t",header = F)

sel.qg.re<-resultats[sel.qg$V1,]
table(sel.qg.re$sign)

cor.qg<-read.table("Q10_bins_significant.tsv",sep="\t",header = T)
cor.qg$ID<-anno[cor.qg$bin,]$ID
row.names(cor.qg)<-cor.qg$ID

sel.qg.re$cor<-cor.qg[row.names(sel.qg.re),]$spearman_rho
library(dplyr)  

sel.qg.re<- sel.qg.re %>% mutate(cor.s=ifelse(cor>0,"pos","neg"))
  
  
  
  
