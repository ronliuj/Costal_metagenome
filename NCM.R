##Fits the neutral model from Sloan et al. 2006 to an OTU table and returns several fitting statistics. Alternatively, will return predicted occurrence frequencies for each OTU based on their abundance in the metacommunity
#Install the following packages if they haven't been availabled in your computer yet 
library(Hmisc)
library(minpack.lm)
library(stats4)

#using Non-linear least squares (NLS) to calculate R2:
#spp: A community table with taxa as rows and samples as columns

abund<-  read.table("data/abund.txt",sep="\t",header = T)
row.names(abund)<-abund$Genomic.bins

anno<-  read.table("data/anno.txt",sep="\t",header = T)
row.names(anno)<-anno$genome

meta<-  read.table("data/para-0528.txt",sep="\t",header = T)
row.names(meta)<-meta$ID

spp<-t(abund[,3:118])

# #get the mean of abundance of each sample
N <- mean(apply(spp, 1, sum))
# #get the mean of species relative abundance in metacommmunity
 p.m <- apply(spp, 2, mean)
 p.m <- p.m[p.m != 0]
# #get the percentage of each species in metacommmunity
 p <- p.m/N


#get the binary data of community abundance matrix
spp.bi <- 1*(spp>0)
#get the frequncy of species occurrence in metacommunity
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]
#get a table record species percentage and occurrence frequency in metacommunity
C <- merge(p, freq, by=0)
#sort the table according to occurence frquency of each species
C <- C[order(C[,2]),]
C <- as.data.frame(C)
#delete rows containning zero
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]
d = 1/N
##Fit model parameter m (or Nm) using Non-linear least squares (NLS)
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
m.fit #get the m value
m.ci <- confint(m.fit, 'm', level=0.95)
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
# get the R2 value
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
Rsqr

#data visulization
#Drawing the figure using grid package:
#p is the mean relative abundance
#freq is occurrence frequency
#freq.pred is predicted occurrence frequency
bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])

bacnlsALL$ID<-anno[row.names(bacnlsALL),]$ID
write.table(bacnlsALL,"ncmresult.txt",sep="\t",row.names = T,quote=F)

inter.col<-rep('black',nrow(bacnlsALL))
inter.col[bacnlsALL$freq <= bacnlsALL$Lower]<-'#A52A2A'#define the color of below points
inter.col[bacnlsALL$freq >= bacnlsALL$Upper]<-'#29A6A6'#define the color of up points
library(grid)
grid.newpage()
pushViewport(viewport(h=0.6,w=0.6))
pushViewport(dataViewport(xData=range(log10(bacnlsALL$p)), yData=c(0,1.02),extension=c(0.02,0)))
grid.rect()
grid.points(log10(bacnlsALL$p), bacnlsALL$freq,pch=20,gp=gpar(col=inter.col,cex=0.7))
grid.yaxis()
grid.xaxis()
grid.lines(log10(bacnlsALL$p),bacnlsALL$freq.pred,gp=gpar(col='blue',lwd=2),default='native')

grid.lines(log10(bacnlsALL$p),bacnlsALL$Lower ,gp=gpar(col='blue',lwd=2,lty=2),default='native') 
grid.lines(log10(bacnlsALL$p),bacnlsALL$Upper,gp=gpar(col='blue',lwd=2,lty=2),default='native')  
grid.text(y=unit(0,'npc')-unit(2.5,'lines'),label='Mean Relative Abundance (log10)', gp=gpar(fontface=2)) 
grid.text(x=unit(0,'npc')-unit(3,'lines'),label='Frequency of Occurance',gp=gpar(fontface=2),rot=90) 
#grid.text(x=unit(0,'npc')-unit(-1,'lines'), y=unit(0,'npc')-unit(-15,'lines'),label='Mean Relative Abundance (log)', gp=gpar(fontface=2)) 
#grid.text(round(coef(m.fit)*N),x=unit(0,'npc')-unit(-5,'lines'), y=unit(0,'npc')-unit(-15,'lines'),gp=gpar(fontface=2)) 
#grid.text(label = "Nm=",x=unit(0,'npc')-unit(-3,'lines'), y=unit(0,'npc')-unit(-15,'lines'),gp=gpar(fontface=2))
#grid.text(round(Rsqr,2),x=unit(0,'npc')-unit(-5,'lines'), y=unit(0,'npc')-unit(-16,'lines'),gp=gpar(fontface=2))
#grid.text(label = "Rsqr=",x=unit(0,'npc')-unit(-3,'lines'), y=unit(0,'npc')-unit(-16,'lines'),gp=gpar(fontface=2))
draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=",round(Rsqr,3),"\n","Nm=",round(coef(m.fit)*N)), x=x[j], y=y[i], just=just)
  #grid.text(deparse(substitute(just)), x=x[j], y=y[i] + unit(2, "lines"),
  #          gp=gpar(col="grey", fontsize=8))
}
x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")
draw.text(c("centre", "bottom"), 4, 1)

##MF

sppMF<-t(abund[,row.names(meta[meta$Category=="MF",])])

# #get the mean of abundance of each sample
N.MF <- mean(apply(sppMF, 1, sum))
# #get the mean of species relative abundance in metacommmunity
p.mMF <- apply(sppMF, 2, mean)
p.mMF <- p.mMF[p.mMF != 0]
# #get the percentage of each species in metacommmunity
p.MF <- p.mMF/N.MF


#get the binary data of community abundance matrix
sppMF.bi <- 1*(sppMF>0)
#get the frequncy of species occurrence in metacommunity
freqMF <- apply(sppMF.bi, 2, mean)
freqMF <- freqMF[freqMF != 0]
#get a table record species percentage and occurrence frequency in metacommunity
CMF <- merge(p.MF, freqMF, by=0)
#sort the table according to occurence frquency of each species
CMF <- CMF[order(CMF[,2]),]
CMF <- as.data.frame(CMF)
#delete rows containning zero
C.0MF <- CMF[!(apply(CMF, 1, function(y) any(y == 0))),]
pMF <- C.0MF[,2]
freqMF <- C.0MF[,3]
names(pMF) <- C.0MF[,1]
names(freqMF) <- C.0MF[,1]
dMF = 1/N.MF
##Fit model parameter m (or Nm) using Non-linear least squares (NLS)
m.fitMF <- nlsLM(freqMF ~ pbeta(dMF, N.MF*m*pMF, N.MF*m*(1 -pMF), lower.tail=FALSE),start=list(m=0.1))
m.fitMF #get the m value
m.ciMF <- confint(m.fitMF, 'm', level=0.95)
freq.predMF <- pbeta(dMF, N.MF*coef(m.fitMF)*pMF, N.MF*coef(m.fitMF)*(1 -pMF), lower.tail=FALSE)
pred.ciMF <- binconf(freq.predMF*nrow(sppMF), nrow(sppMF), alpha=0.05, method="wilson", return.df=TRUE)
# get the R2 value
RsqrMF <- 1 - (sum((freqMF - freq.predMF)^2))/(sum((freqMF - mean(freqMF))^2))
RsqrMF

#data visulization
#Drawing the figure using grid package:
#p is the mean relative abundance
#freq is occurrence frequency
#freq.pred is predicted occurrence frequency
bacnlsALLMF <-data.frame(pMF,freqMF,freq.predMF,pred.ciMF[,2:3])
inter.colMF<-rep('black',nrow(bacnlsALLMF))
inter.colMF[bacnlsALLMF$freqMF <= bacnlsALLMF$Lower]<-'#A52A2A'#define the color of below points
inter.colMF[bacnlsALLMF$freqMF >= bacnlsALLMF$Upper]<-'#29A6A6'#define the color of up points
library(grid)
grid.newpage()
pushViewport(viewport(h=0.6,w=0.6))
pushViewport(dataViewport(xData=range(log10(bacnlsALLMF$pMF)), yData=c(0,1.02),extension=c(0.02,0)))
grid.rect()
grid.points(log10(bacnlsALLMF$pMF), bacnlsALLMF$freqMF,pch=20,gp=gpar(col=inter.colMF,cex=0.7))
grid.yaxis()
grid.xaxis()
grid.lines(log10(bacnlsALLMF$pMF),bacnlsALLMF$freq.predMF,gp=gpar(col='blue',lwd=2),default='native')

grid.lines(log10(bacnlsALLMF$pMF),bacnlsALLMF$Lower ,gp=gpar(col='blue',lwd=2,lty=2),default='native') 
grid.lines(log10(bacnlsALLMF$pMF),bacnlsALLMF$Upper,gp=gpar(col='blue',lwd=2,lty=2),default='native')  
grid.text(y=unit(0,'npc')-unit(2.5,'lines'),label='Mean Relative Abundance (log10)', gp=gpar(fontface=2)) 
grid.text(x=unit(0,'npc')-unit(3,'lines'),label='Frequency of Occurance',gp=gpar(fontface=2),rot=90) 
#grid.text(x=unit(0,'npc')-unit(-1,'lines'), y=unit(0,'npc')-unit(-15,'lines'),label='Mean Relative Abundance (log)', gp=gpar(fontface=2)) 
#grid.text(round(coef(m.fit)*N),x=unit(0,'npc')-unit(-5,'lines'), y=unit(0,'npc')-unit(-15,'lines'),gp=gpar(fontface=2)) 
#grid.text(label = "Nm=",x=unit(0,'npc')-unit(-3,'lines'), y=unit(0,'npc')-unit(-15,'lines'),gp=gpar(fontface=2))
#grid.text(round(Rsqr,2),x=unit(0,'npc')-unit(-5,'lines'), y=unit(0,'npc')-unit(-16,'lines'),gp=gpar(fontface=2))
#grid.text(label = "Rsqr=",x=unit(0,'npc')-unit(-3,'lines'), y=unit(0,'npc')-unit(-16,'lines'),gp=gpar(fontface=2))
draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=",round(RsqrMF,3),"\n","Nm=",round(coef(m.fitMF)*N.MF)), x=x[j], y=y[i], just=just)
  #grid.text(deparse(substitute(just)), x=x[j], y=y[i] + unit(2, "lines"),
  #          gp=gpar(col="grey", fontsize=8))
}
x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")
draw.text(c("centre", "bottom"), 4, 1)

##SM

sppSM<-t(abund[,row.names(meta[meta$Category=="SM",])])

# #get the mean of abundance of each sample
N.SM <- mean(apply(sppSM, 1, sum))
# #get the mean of species relative abundance in metacommmunity
p.mSM <- apply(sppSM, 2, mean)
p.mSM <- p.mSM[p.mSM != 0]
# #get the percentage of each species in metacommmunity
p.SM <- p.mSM/N.SM


#get the binary data of community abundance matrix
sppSM.bi <- 1*(sppSM>0)
#get the frequncy of species occurrence in metacommunity
freqSM <- apply(sppSM.bi, 2, mean)
freqSM <- freqSM[freqSM != 0]
#get a table record species percentage and occurrence frequency in metacommunity
CSM <- merge(p.SM, freqSM, by=0)
#sort the table according to occurence frquency of each species
CSM <- CSM[order(CSM[,2]),]
CSM <- as.data.frame(CSM)
#delete rows containning zero
C.0SM <- CSM[!(apply(CSM, 1, function(y) any(y == 0))),]
pSM <- C.0SM[,2]
freqSM <- C.0SM[,3]
names(pSM) <- C.0SM[,1]
names(freqSM) <- C.0SM[,1]
dSM = 1/N.SM
##Fit model parameter m (or Nm) using Non-linear least squares (NLS)
m.fitSM <- nlsLM(freqSM ~ pbeta(dSM, N.SM*m*pSM, N.SM*m*(1 -pSM), lower.tail=FALSE),start=list(m=0.1))
m.fitSM #get the m value
m.ciSM <- confint(m.fitSM, 'm', level=0.95)
freq.predSM <- pbeta(dSM, N.SM*coef(m.fitSM)*pSM, N.SM*coef(m.fitSM)*(1 -pSM), lower.tail=FALSE)
pred.ciSM <- binconf(freq.predSM*nrow(sppSM), nrow(sppSM), alpha=0.05, method="wilson", return.df=TRUE)
# get the R2 value
RsqrSM <- 1 - (sum((freqSM - freq.predSM)^2))/(sum((freqSM - mean(freqSM))^2))
RsqrSM

#data visulization
#Drawing the figure using grid package:
#p is the mean relative abundance
#freq is occurrence frequency
#freq.pred is predicted occurrence frequency
bacnlsALLSM <-data.frame(pSM,freqSM,freq.predSM,pred.ciSM[,2:3])
inter.colSM<-rep('black',nrow(bacnlsALLSM))
inter.colSM[bacnlsALLSM$freqSM <= bacnlsALLSM$Lower]<-'#A52A2A'#define the color of below points
inter.colSM[bacnlsALLSM$freqSM >= bacnlsALLSM$Upper]<-'#29A6A6'#define the color of up points
library(grid)
grid.newpage()
pushViewport(viewport(h=0.6,w=0.6))
pushViewport(dataViewport(xData=range(log10(bacnlsALLSM$pSM)), yData=c(0,1.02),extension=c(0.02,0)))
grid.rect()
grid.points(log10(bacnlsALLSM$pSM), bacnlsALLSM$freqSM,pch=20,gp=gpar(col=inter.colSM,cex=0.7))
grid.yaxis()
grid.xaxis()
grid.lines(log10(bacnlsALLSM$pSM),bacnlsALLSM$freq.predSM,gp=gpar(col='blue',lwd=2),default='native')

grid.lines(log10(bacnlsALLSM$pSM),bacnlsALLSM$Lower ,gp=gpar(col='blue',lwd=2,lty=2),default='native') 
grid.lines(log10(bacnlsALLSM$pSM),bacnlsALLSM$Upper,gp=gpar(col='blue',lwd=2,lty=2),default='native')  
grid.text(y=unit(0,'npc')-unit(2.5,'lines'),label='Mean Relative Abundance (log10)', gp=gpar(fontface=2)) 
grid.text(x=unit(0,'npc')-unit(3,'lines'),label='Frequency of Occurance',gp=gpar(fontface=2),rot=90) 
#grid.text(x=unit(0,'npc')-unit(-1,'lines'), y=unit(0,'npc')-unit(-15,'lines'),label='Mean Relative Abundance (log)', gp=gpar(fontface=2)) 
#grid.text(round(coef(m.fit)*N),x=unit(0,'npc')-unit(-5,'lines'), y=unit(0,'npc')-unit(-15,'lines'),gp=gpar(fontface=2)) 
#grid.text(label = "Nm=",x=unit(0,'npc')-unit(-3,'lines'), y=unit(0,'npc')-unit(-15,'lines'),gp=gpar(fontface=2))
#grid.text(round(Rsqr,2),x=unit(0,'npc')-unit(-5,'lines'), y=unit(0,'npc')-unit(-16,'lines'),gp=gpar(fontface=2))
#grid.text(label = "Rsqr=",x=unit(0,'npc')-unit(-3,'lines'), y=unit(0,'npc')-unit(-16,'lines'),gp=gpar(fontface=2))
draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=",round(RsqrSM,3),"\n","Nm=",round(coef(m.fitSM)*N.SM)), x=x[j], y=y[i], just=just)
  #grid.text(deparse(substitute(just)), x=x[j], y=y[i] + unit(2, "lines"),
  #          gp=gpar(col="grey", fontsize=8))
}
x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")
draw.text(c("centre", "bottom"), 4, 1)


##SG

sppSG<-t(abund[,row.names(meta[meta$Category=="SG",])])

# #get the mean of abundance of each sample
N.SG <- mean(apply(sppSG, 1, sum))
# #get the mean of species relative abundance in metacommmunity
p.mSG <- apply(sppSG, 2, mean)
p.mSG <- p.mSG[p.mSG != 0]
# #get the percentage of each species in metacommmunity
p.SG <- p.mSG/N.SG


#get the binary data of community abundance matrix
sppSG.bi <- 1*(sppSG>0)
#get the frequncy of species occurrence in metacommunity
freqSG <- apply(sppSG.bi, 2, mean)
freqSG <- freqSG[freqSG != 0]
#get a table record species percentage and occurrence frequency in metacommunity
CSG <- merge(p.SG, freqSG, by=0)
#sort the table according to occurence frquency of each species
CSG <- CSG[order(CSG[,2]),]
CSG <- as.data.frame(CSG)
#delete rows containning zero
C.0SG <- CSG[!(apply(CSG, 1, function(y) any(y == 0))),]
pSG <- C.0SG[,2]
freqSG <- C.0SG[,3]
names(pSG) <- C.0SG[,1]
names(freqSG) <- C.0SG[,1]
dSG = 1/N.SG
##Fit model parameter m (or Nm) using Non-linear least squares (NLS)
m.fitSG <- nlsLM(freqSG ~ pbeta(dSG, N.SG*m*pSG, N.SG*m*(1 -pSG), lower.tail=FALSE),start=list(m=0.1))
m.fitSG #get the m value
m.ciSG <- confint(m.fitSG, 'm', level=0.95)
freq.predSG <- pbeta(dSG, N.SG*coef(m.fitSG)*pSG, N.SG*coef(m.fitSG)*(1 -pSG), lower.tail=FALSE)
pred.ciSG <- binconf(freq.predSG*nrow(sppSG), nrow(sppSG), alpha=0.05, method="wilson", return.df=TRUE)
# get the R2 value
RsqrSG <- 1 - (sum((freqSG - freq.predSG)^2))/(sum((freqSG - mean(freqSG))^2))
RsqrSG

#data visulization
#Drawing the figure using grid package:
#p is the mean relative abundance
#freq is occurrence frequency
#freq.pred is predicted occurrence frequency
bacnlsALLSG <-data.frame(pSG,freqSG,freq.predSG,pred.ciSG[,2:3])
inter.colSG<-rep('black',nrow(bacnlsALLSG))
inter.colSG[bacnlsALLSG$freqSG <= bacnlsALLSG$Lower]<-'#A52A2A'#define the color of below points
inter.colSG[bacnlsALLSG$freqSG >= bacnlsALLSG$Upper]<-'#29A6A6'#define the color of up points
library(grid)
grid.newpage()
pushViewport(viewport(h=0.6,w=0.6))
pushViewport(dataViewport(xData=range(log10(bacnlsALLSG$pSG)), yData=c(0,1.02),extension=c(0.02,0)))
grid.rect()
grid.points(log10(bacnlsALLSG$pSG), bacnlsALLSG$freqSG,pch=20,gp=gpar(col=inter.colSG,cex=0.7))
grid.yaxis()
grid.xaxis()
grid.lines(log10(bacnlsALLSG$pSG),bacnlsALLSG$freq.predSG,gp=gpar(col='blue',lwd=2),default='native')

grid.lines(log10(bacnlsALLSG$pSG),bacnlsALLSG$Lower ,gp=gpar(col='blue',lwd=2,lty=2),default='native') 
grid.lines(log10(bacnlsALLSG$pSG),bacnlsALLSG$Upper,gp=gpar(col='blue',lwd=2,lty=2),default='native')  
grid.text(y=unit(0,'npc')-unit(2.5,'lines'),label='Mean Relative Abundance (log10)', gp=gpar(fontface=2)) 
grid.text(x=unit(0,'npc')-unit(3,'lines'),label='Frequency of Occurance',gp=gpar(fontface=2),rot=90) 
#grid.text(x=unit(0,'npc')-unit(-1,'lines'), y=unit(0,'npc')-unit(-15,'lines'),label='Mean Relative Abundance (log)', gp=gpar(fontface=2)) 
#grid.text(round(coef(m.fit)*N),x=unit(0,'npc')-unit(-5,'lines'), y=unit(0,'npc')-unit(-15,'lines'),gp=gpar(fontface=2)) 
#grid.text(label = "Nm=",x=unit(0,'npc')-unit(-3,'lines'), y=unit(0,'npc')-unit(-15,'lines'),gp=gpar(fontface=2))
#grid.text(round(Rsqr,2),x=unit(0,'npc')-unit(-5,'lines'), y=unit(0,'npc')-unit(-16,'lines'),gp=gpar(fontface=2))
#grid.text(label = "Rsqr=",x=unit(0,'npc')-unit(-3,'lines'), y=unit(0,'npc')-unit(-16,'lines'),gp=gpar(fontface=2))
draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=",round(RsqrSG,3),"\n","Nm=",round(coef(m.fitSG)*N.SG)), x=x[j], y=y[i], just=just)
  #grid.text(deparse(substitute(just)), x=x[j], y=y[i] + unit(2, "lines"),
  #          gp=gpar(col="grey", fontsize=8))
}
x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")
draw.text(c("centre", "bottom"), 4, 1)


all.pie<-data.frame(table(inter.col))

all.pie$part<-c("above","below","neutral")
all.pie$prot<-all.pie$Freq/1913

"#29A6A6" "#A52A2A"   black 
436     539     938 

pie(all.pie$Freq,border="white",col=c("#29A6A6" ,"#A52A2A",   "black"),
                                      lwd=0.5,labels=all.pie$part)


MF.pie<-data.frame(table(inter.colMF))

MF.pie$part<-c("above","below","neutral")
MF.pie$prot<-MF.pie$Freq/1308

pie(MF.pie$Freq,border="white",col=c("#29A6A6" ,"#A52A2A",   "black"),
    lwd=0.5,labels=MF.pie$part)


SM.pie<-data.frame(table(inter.colSM))

SM.pie$part<-c("above","below","neutral")
SM.pie$prot<-SM.pie$Freq/1482

pie(SM.pie$Freq,border="white",col=c("#29A6A6" ,"#A52A2A",   "black"),
    lwd=0.5,labels=SM.pie$part)


SG.pie<-data.frame(table(inter.colSG))

SG.pie$part<-c("above","below","neutral")
SG.pie$prot<-SG.pie$Freq/1172

pie(SG.pie$Freq,border="white",col=c("#29A6A6" ,"#A52A2A",   "black"),
    lwd=0.5,labels=SG.pie$part)

ncmSM<-bacnlsALLSM
colnames(ncmSM)<-c("p","freq","freq.pred","Lower","Upper")

ncmSM$ID<-anno[row.names(ncmSM),]$ID
ncmSM$Category<-"SM"


ncmSG<-bacnlsALLSG
colnames(ncmSG)<-c("p","freq","freq.pred","Lower","Upper")

ncmSG$ID<-anno[row.names(ncmSG),]$ID
ncmSG$Category<-"SG"

ncmMF<-bacnlsALLMF
colnames(ncmMF)<-c("p","freq","freq.pred","Lower","Upper")

ncmMF$ID<-anno[row.names(ncmMF),]$ID
ncmMF$Category<-"MF"

ncm2<-rbind(ncmSM,ncmSG,ncmMF)

write.table(ncm2,"ncmresult_eco.txt",sep="\t",row.names = T,quote=F)






##another method
# 清理工作环境中的所有对象
rm(list = ls())

# 安装并加载所需的包
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
if (!requireNamespace("ggimage", quietly = TRUE)) {
  install.packages("ggimage")
}
if (!requireNamespace("plyr", quietly = TRUE)) {
  install.packages("plyr")
}

library(devtools)
library(ggplot2)
library(ggimage)
library(plyr)


# 加载 remotes 包（用于 install_github 函数）
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
# 检查并安装 MicEco 包
if (!requireNamespace("MicEco", quietly = TRUE)) {
  remotes::install_github("Russel88/MicEco")
}
library(MicEco)

# 读取数据
read_table <- read.csv(file = "read_table.csv", row.names = 1)

# 运行中性模型分析
res <- neutral.fit(t(read_table))

m <- res[[1]][1]
N <- res[[1]][4]
Nm<- N*m
r2 <- res[[1]][3]
out <- res[[2]]

# 处理数据
out$group <- with(out, ifelse(freq < Lower, "#509579",
                              ifelse(freq > Upper, "#cf9198", "#485970")))

# 绘制模型结果图
p1 <- ggplot(data = out) +
  geom_line(aes(x = log(p), y = freq.pred), size = 1.2, linetype = 1) +
  geom_line(aes(x = log(p), y = Lower), size = 1.2, linetype = 2) +
  geom_line(aes(x = log(p), y = Upper), size = 1.2, linetype = 2) +
  geom_point(aes(x = log(p), y = freq, color = group), size = 2) +
  xlab("log10(mean relative abundance)") +
  ylab("Occurrence frequency") +
  scale_colour_manual(values = c("#485970", "#cf9198", "#509579")) +  # 应用新的颜色方案
  annotate("text", x = -4, y = 0.10, label = paste("Nm = ", round(Nm, 3)), size = 7) +
  annotate("text", x = -4, y = 0.20, label = paste("R² = ", round(r2, 3)), size = 7) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line.x = element_line(size = 0.5, colour = "black"),
    axis.line.y = element_line(size = 0.5, colour = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black", size = 18),
    legend.position = "none",
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.text = element_text(size = 18),
    text = element_text(family = "sans", size = 18),
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # 添加图表边框
  )


p1
###########要是需要计算各分类占比可以使用在图中加入饼图，代码如下#########

# 计算饼图数据(下面不清晰的话可以参考)
#low = nrow(out[out[,6]== "#509579",])
#med = nrow(out[out[,6]== "#485970",])
#high = nrow(out[out[,6]== "#cf9198",])
#type <- c('med','high','low')
#nums <- c(med,high,low)
#df <- data.frame(type = type, nums = nums)

data_summary <- as.data.frame(table(out$group))
colnames(data_summary) <- c("group", "nums")
data_summary$type<-(c("med","low","high"))
data_summary$percentage <- round(data_summary$nums / sum(data_summary$nums) * 100, 1)
data_summary$label <- paste(data_summary$type, paste(data_summary$percentage, "%", sep = ''))

# 绘制饼图
p2 <- ggplot(data_summary, aes(x = "", y = nums, fill = group)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(
    values = c("#485970","#509579","#cf9198"),  # 应用新的颜色方案
    labels = data_summary$label
  ) +
  theme_void() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.text = element_text(size = 18)
  )
p2

# 嵌入饼图到主图中
p_final <- p1 + geom_subview(subview = p2, x = -11.5, y = 0.95, w = 4.5, h = 4.5)

# 显示最终图形
print(p_final)





















