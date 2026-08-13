##生态位
library(spaa)
##KO 
KO<-read.table("KO.sum.TPM",header = T,sep="\t")
row.names(KO)<-KO$Name

comm.tab<-data.frame(t(KO[,2:117]))
comm.tab<-ceiling(comm.tab*10)##注意这里不能太大，会超过32位int限制

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

abund<-read.table("../KO.sum.TPM",header = T,sep="\t")
row.names(abund)<-abund$Name

comm.tab<-data.frame(t(abund[,2:117]))

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
write.table(resultats,"levin-stats_KO.txt",row.names = T,sep="\t",quote=F)


## on pc own
abund<-read.table("../KO.sum.TPM",header = T,sep="\t")
row.names(abund)<-abund$ID

comm.tab<-data.frame(t(abund[,2:117]))
comm.tab<-ceiling(comm.tab*10)##注意这里不能太大，会超过32位int限制

library(spaa)
library(vegan)
levin.index.real<-as.numeric(niche.width(comm.tab,method="levins"))
names(levin.index.real)<-colnames(comm.tab)

comm.tab<-comm.tab>0

Prevalence<-colSums(comm.tab)/116



resultats<-read.table("levin-stats_KO.txt",header = T,sep="\t")

resultats$Prevalence<-Prevalence

write.table(resultats,"levin-stats_KO2.txt",sep="\t",quote = F)

resultats$abund<-log10(rowMeans(abund[,2:117]))

resultats$sign[resultats$sign=="NON SIGNIFICANT"]<-"NON_SIGNIFICANT"

resultats$sign<-factor(resultats$sign,levels =c("GENERALIST","SPECIALIST","NON_SIGNIFICANT"))

library(RColorBrewer)
brewer.pal(11,"Set2")

[1] "#66C2A5" "#FC8D62" "#8DA0CB" "#E78AC3" "#A6D854" "#FFD92F" "#E5C494" "#B3B3B3"

library(ggplot2)
p1<-ggplot()+
  #geom_point(aes(x=Prevalence,y=abund,color=sign,alpha=sign),size=2.5)+
  geom_jitter(data=resultats[resultats$sign=="NON_SIGNIFICANT",],
              aes(x=Prevalence,y=abund,color=sign,alpha=sign),
              size = 2,width=0.003,shape=16)+
  geom_jitter(data=resultats[resultats$sign=="SPECIALIST",],
              aes(x=Prevalence,y=abund,color=sign,alpha=sign),
              size = 2,width=0.003,shape=16)+
  geom_jitter(data=resultats[resultats$sign=="GENERALIST",],
              aes(x=Prevalence,y=abund,color=sign,alpha=sign),
              size = 2,width=0.003,shape=16)+
  annotate("text",x=0.5,y=-2,
           label=paste0("GENERALIST = 1865","\n",
                        "SPECIALIST = 8670","\n",
                        "NON_SIGNIFICANT = 2757","\n"))+
  scale_alpha_manual(values=c(SPECIALIST = 0.5, GENERALIST = 1, NON_SIGNIFICANT= 0.5 ))+
  scale_color_manual(values=c(SPECIALIST = "#FC8D62", GENERALIST = "#66C2A5", NON_SIGNIFICANT= "lightgrey" ))+
  theme_bw()+ylab("log10(abundance)")+xlab("Prevalence")+#ylim(c(-4.3,2))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))
#export 6x10

KO<-read.table("KO.sum.TPM",header = T,sep="\t")
row.names(KO)<-KO$Name

comm.tab<-data.frame(t(KO[,2:117]))
comm.tab<-ceiling(comm.tab*10)

comm.tab<-data.frame(t(comm.tab))

library(vegan)
comm.tab2 <-decostand(comm.tab, MARGIN=2, method="total")


KO.sim <-simper(t(comm.tab2))
KO.sim.res <-data.frame(summary(KO.sim)$total)
KO.sim.res <-KO.sim.res[order(rownames(KO.sim.res)),]
summary(rownames(KO.sim.res)==rownames(resultats))

KO.sim.res$spc.gen <-resultats$sign
KO.sim.res2 <-KO.sim.res[KO.sim.res$average>0,]

KO.sim.res2.sum<-aggregate(average~spc.gen,data=KO.sim.res2,sum)

KO.sim.res2.mean<-aggregate(average~spc.gen,data=KO.sim.res2,mean)

KO.sim.res2.sd<-aggregate(average~spc.gen,data=KO.sim.res2,sd)

KO.sim.res.sum<-aggregate(average~spc.gen,data=KO.sim.res,sum)

boxplot(KO.sim.res2$average~KO.sim.res2$spc.gen, xlab="", log="y", ylab="Contribution to Bray-Curtis dissimilarity", col="white")


write.table(KO.sim.res2,"contribution.txt",sep="\t",quote=F)

kruskal.test(average~spc.gen, data=KO.sim.res2)#< 2.2e-16
PMCMRplus::kwAllPairsDunnTest(average~spc.gen, data=KO.sim.res2, 
                                  p.adjust.method ="bonferroni")






resultats2<-resultats[resultats$Prevalence>0.5,]
table(resultats2$sign)

save.image("Levin.Genral_KO.Rdata")
