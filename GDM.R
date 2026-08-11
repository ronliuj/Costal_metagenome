##data from mcm

#install.packages("gdm")
library(gdm)
library(dplyr)
##MAG_fun
# library(foreach)
# library(doParallel)
# detectCores()#40
# registerDoParallel(cores = 30)

##grep search mapids
##collect mapids in each sample

aa<-read.table(paste("data/mag_fun/",row.names(abund)[1],".txt",sep=""),header = F)
aa[,2:117]<-matrix(rep(as.numeric(abund[1,3:118]),nrow(aa)),ncol=116,nrow=nrow(aa),byrow = T)
aa<-aggregate(.~V1,data=aa,sum,na.rm=T,na.action = NULL)

# foreach (i=2:1954) %dopar%{
  
for(i in 2:1954 ){
  bb<-read.table(paste("data/mag_fun/",row.names(abund)[i],".txt",sep=""),header = F)
  bb[,2:117]<-matrix(rep(as.numeric(abund[row.names(abund)[i],3:118]),nrow(bb)),ncol=116,nrow=nrow(bb),byrow = T)
  aa<-rbind(aa,bb)
  aa<-aggregate(.~V1,data=aa,sum,na.rm=T,na.action = NULL)
}


MAG.fun<-aa
rm(aa)

row.names(MAG.fun)<-MAG.fun$V1

colnames(MAG.fun)<-c("MAPID",colnames(abund)[3:118])

MAG.fun$MAPID<-gsub("map","ko",MAG.fun$MAPID)

row.names(MAG.fun)<-MAG.fun$MAPID

##read human factors

Human<-read.table("data/Human_factor.txt",header=T,sep="\t")
row.names(Human)<-Human$ID

meta[,86:93]<-Human[row.names(meta),6:13]


meta<-meta %>% mutate(Vanillyl_real=Vanillyl/0.333,
                      Syringyl_real=Syringyl/0.9,
                      Lignin.phenol_real=Cinnamyl+Vanillyl_real+Syringyl_real)


climate<-read.table("data/climate.txt",header = T,sep="\t")
row.names(climate)<-climate$ID
  
meta[,97:98]<-climate[row.names(meta),5:6]

write.table(meta,"meta_116-250124.txt",sep="\t",quote=F)

selected<-c("TOC", "POC",  "MOC", "DOC", "Ca.OC", "Fe.OC","Lignin.phenol_real",
"Cinnamyl",  "Syringyl_real",  "Vanillyl_real" ,"Microbial.necromass",
"Bacterial.necromass" ,"Fungal.necromass" ,   
"pH" , "EC" , "TN", "Fed",   "Feo", "Fep" , "MeanSize" ,"X.4", "X4.63" ,"X.63" ,"C_N",
"Ald", "Alo" , "Alp",  "C.P",  "P" ,"MAP" ,"MAT",
"Lon", "Lat",
"Subarea","Population", "Density", "GDP","Agriculture", "Public_in","Public_out","Interprise")

#is.element(selected,colnames(meta))

meta_sel<-meta[,selected]  

##species should be in site-species struc, the row should be the same as meta
##species and meta must have the same colume of sites
abundT<-data.frame(t(abund[,3:118]))
abundT<-abundT[row.names(meta),]

# library(vegan)
# abundT <- decostand(abundT, method = 'hellinger')

abundT$sites<-row.names(abundT)


meta_sel_scale<-data.frame(scale(meta_sel[,1:41]))

library(mice)
##修改不使用scaled数据
# 使用 mice 进行多重插补
#imp <- mice(scale(rf.d[,3:26]), method = 'pmm', m = 5)
imp <- mice(meta_sel, method = 'rf', m = 5)

meta_filled <- complete(imp, 1)

meta_filled$sites<-row.names(meta_filled)

# abundT.F<-data.frame(abundT[,1:1954]!=0)
# abundT.F<-data.frame(apply(abundT.F,1,as.numeric))
# abundT.F$sites<-abundT$sites
# 
# ##abundance = T, 丰度数据
gdmTab <- formatsitepair(bioData=abundT,
                         bioFormat=1,
                         predData=meta_filled,
                         XColumn="Lon", 
                         YColumn="Lat",
                         siteColumn="sites",abundance = T)

gdm.1 <- gdm(data=gdmTab, geo=TRUE,)
summary(gdm.1) ## 45.69906

## default abundance = F, presence 数据
gdmTab2 <- formatsitepair(bioData=abundT,
                         bioFormat=1,
                         predData=meta_filled,
                         XColumn="Lon", 
                         YColumn="Lat",
                         siteColumn="sites")

gdm.2 <- gdm(data=gdmTab2, geo=TRUE)
summary(gdm.2)#39.1927

##pred数据必须包含sites和坐标数据
all<-c(selected,"sites")
carbon<-c(selected[1:13],"sites","Lon","Lat")
carbon<-carbon[c(1,4,2,3,6,5,7,11,14,15,16)]

soil<-c(selected[14:29],"sites","Lon","Lat")
soil<-soil[c(1,2,3,11,16,15,4,5,6,12,13,14,8,9,10,17,18,19)]

climate<-c(selected[30:31],"sites","Lon","Lat")

human<-c(selected[36:40],"sites","Lon","Lat")
human<-human[c(1,2,3,6,7,8)]

list.par<-NULL
list.par[[1]]<-all
list.par[[2]]<-carbon
list.par[[3]]<-soil
list.par[[4]]<-climate
list.par[[5]]<-human
pars<-c("all","carbon","soil","climate","human")

##物种丰度数据
list.gdm.abund.mag<-NULL
for(i in 1:5){
  gdmTab0 <- formatsitepair(bioData=abundT,
                           bioFormat=1,
                           predData=meta_filled[,list.par[[i]]],
                           XColumn="Lon", 
                           YColumn="Lat",
                           siteColumn="sites",abundance = T)
  
  list.gdm.abund.mag[[i]] <- gdm(data=gdmTab0, geo=TRUE)
  list.gdm.abund.mag[[i]]$dataname<-paste("MAG.abund.",pars[i],sep="")
 }
  
##物种presence数据
list.gdm.AbPre.mag<-NULL
for(i in 1:5){
  gdmTab0 <- formatsitepair(bioData=abundT,
                            bioFormat=1,
                            predData=meta_filled[,list.par[[i]]],
                            XColumn="Lon", 
                            YColumn="Lat",
                            siteColumn="sites")
  
  list.gdm.AbPre.mag[[i]] <- gdm(data=gdmTab0, geo=TRUE)
  list.gdm.AbPre.mag[[i]]$dataname<-paste("MAG.absense.",pars[i],sep="")
}

## MAG function distance
MAG.fun<-data.frame(t(MAG.fun[,2:117]))
MAG.fun<-MAG.fun[row.names(meta),]
MAG.fun$sites<-row.names(MAG.fun)

##功能丰度数据

list.gdm.abund.FUN<-NULL
for(i in 1:5){
  gdmTab0 <- formatsitepair(bioData=MAG.fun,
                            bioFormat=1,
                            predData=meta_filled[,list.par[[i]]],
                            XColumn="Lon", 
                            YColumn="Lat",
                            siteColumn="sites",abundance = T)
  
  list.gdm.abund.FUN[[i]] <- gdm(data=gdmTab0, geo=TRUE)
  list.gdm.abund.FUN[[i]]$dataname<-paste("MAG.abund.",pars[i],sep="")
}

##功能presence数据
# list.gdm.AbPre.FUN<-NULL
# for(i in 2:5){
#   gdmTab0 <- formatsitepair(bioData=MAG.fun,
#                             bioFormat=1,
#                             predData=meta_filled[,list.par[[2]]],
#                             XColumn="Lon", 
#                             YColumn="Lat",
#                             siteColumn="sites")
#   
#   list.gdm.AbPre.FUN[[2]] <- gdm(data=gdmTab0, geo=TRUE)
#   list.gdm.AbPre.FUN[[i]]$dataname<-paste("MAG.absense.",pars[i],sep="")
# }

##功能presence情况比较一致，无结果


## read KEGG levels
kegglvls<-read.table("data/KEGG_levels.txt",header = T,sep="\t")

metabolisms<-kegglvls[kegglvls$Level.1=="Metabolism",]

metalevel<-unique(metabolisms$Level.2)

list.gdb.metabo<-NULL
devia<-data.frame(matrix(nrow=64,ncol=3))

n=1
for(i in 1:12){
  aa<-c("sites",metabolisms[metabolisms$Level.2==metalevel[i],]$KID)
  aa<-aa[is.element(aa,colnames(MAG.fun))]
  for(j in 2:5){
  if(length(aa)>2){
  MAG.fun0<-MAG.fun[,aa]
  gdmTab0 <- formatsitepair(bioData=MAG.fun0,
                            bioFormat=1,
                            predData=meta_filled[,list.par[[j]]],
                            XColumn="Lon", 
                            YColumn="Lat",
                            siteColumn="sites",abundance = T)
  
  list.gdb.metabo[[n]] <- gdm(data=gdmTab0, geo=TRUE)
  list.gdb.metabo[[n]]$dataname<-paste(metalevel[i],pars[i],sep="")
  devia[n,1]<-metalevel[i]
  devia[n,2]<-pars[j]
  devia[n,3]<-list.gdb.metabo[[n]]$explained
    } else {
    list.gdb.metabo[[n]]<-paste(metalevel[i],pars[i],sep="")
    devia[n,1]<-metalevel[i]
    devia[n,2]<-pars[j]
    devia[n,3]<-NA
  }
    n=n+1  
    }
  
}

colnames(devia)<-c("Micro","Predict","Variance.explained")


list.gdb.metabo<-c(list.gdb.metabo,list.gdm.abund.mag[2:5])
for(i in 2:5){
  devia[47+i,1]<-"MAG-abund"
  devia[47+i,2]<-pars[i]
  devia[47+i,3]<-list.gdm.abund.mag[[i]]$explained
}
  
list.gdb.metabo<-c(list.gdb.metabo,list.gdm.AbPre.mag[2:5])
for(i in 2:5){
  devia[51+i,1]<-"MAG-presence"
  devia[51+i,2]<-pars[i]
  devia[51+i,3]<-list.gdm.AbPre.mag[[i]]$explained
}

list.gdb.metabo<-c(list.gdb.metabo,list.gdm.abund.FUN[2:5])
for(i in 2:5){
  devia[55+i,1]<-"Fun-abund"
  devia[55+i,2]<-pars[i]
  devia[55+i,3]<-list.gdm.abund.FUN[[i]]$explained
}


devia$Index<-1:64

devia$Cate<-"all"



##SM
SM.sams<-row.names(meta[meta$Category=="SM",])
##物种丰度数据
for(i in 2:5){
  gdmTab0 <- formatsitepair(bioData=abundT[SM.sams,],
                            bioFormat=1,
                            predData=meta_filled[SM.sams,list.par[[i]]],
                            XColumn="Lon", 
                            YColumn="Lat",
                            siteColumn="sites",abundance = T)
  
  list.gdb.metabo[[i+59]] <- gdm(data=gdmTab0, geo=TRUE)
  list.gdb.metabo[[i+59]]$dataname<-paste("MAG.abund.",pars[i],sep="")
  devia[59+i,1]<-"MAG-abund"
  devia[59+i,2]<-pars[i]
  devia[59+i,3]<-list.gdb.metabo[[i+59]]$explained
  devia[59+i,4]<-59+i
  devia[59+i,5]<-"SM"
}


for(i in 2:5){
  gdmTab0 <- formatsitepair(bioData=abundT[SM.sams,],
                            bioFormat=1,
                            predData=meta_filled[SM.sams,list.par[[i]]],
                            XColumn="Lon", 
                            YColumn="Lat",
                            siteColumn="sites")
  
  list.gdb.metabo[[i+63]] <- gdm(data=gdmTab0, geo=TRUE)
  list.gdb.metabo[[i+63]]$dataname<-paste("MAG-presence.",pars[i],sep="")
  devia[63+i,1]<-"MAG-presence"
  devia[63+i,2]<-pars[i]
  devia[63+i,3]<-list.gdb.metabo[[i+63]]$explained
  devia[63+i,4]<-63+i
  devia[63+i,5]<-"SM"
}



##功能丰度数据

for(i in 2:5){
  gdmTab0 <- formatsitepair(bioData=MAG.fun[SM.sams,],
                            bioFormat=1,
                            predData=meta_filled[SM.sams,list.par[[i]]],
                            XColumn="Lon", 
                            YColumn="Lat",
                            siteColumn="sites",abundance = T)
  
  list.gdb.metabo[[i+67]] <- gdm(data=gdmTab0, geo=TRUE)
  list.gdb.metabo[[i+67]]$dataname<-paste("Fun-abund.",pars[i],sep="")
  devia[67+i,1]<-"Fun-abund"
  devia[67+i,2]<-pars[i]
  devia[67+i,3]<-list.gdb.metabo[[i+67]]$explained
  devia[67+i,4]<-67+i
  devia[67+i,5]<-"SM"
}

n=73
for(i in 1:12){
  aa<-c("sites",metabolisms[metabolisms$Level.2==metalevel[i],]$KID)
  aa<-aa[is.element(aa,colnames(MAG.fun))]
  for(j in 2:5){
    if(length(aa)>2){
      MAG.fun0<-MAG.fun[SM.sams,aa]
      gdmTab0 <- formatsitepair(bioData=MAG.fun0,
                                bioFormat=1,
                                predData=meta_filled[SM.sams,list.par[[j]]],
                                XColumn="Lon", 
                                YColumn="Lat",
                                siteColumn="sites",abundance = T)
      
      list.gdb.metabo[[n]] <- gdm(data=gdmTab0, geo=TRUE)
      list.gdb.metabo[[n]]$dataname<-paste(metalevel[i],pars[i],sep="")
      devia[n,1]<-metalevel[i]
      devia[n,2]<-pars[j]
      devia[n,3]<-list.gdb.metabo[[n]]$explained
      devia[n,4]<-n
      devia[n,5]<-"SM"
    } else {
      list.gdb.metabo[[n]]<-paste(metalevel[i],pars[i],sep="")
      devia[n,1]<-metalevel[i]
      devia[n,2]<-pars[j]
      devia[n,3]<-NA
      devia[n,4]<-n
      devia[n,5]<-"SM"
    }
    n=n+1  
  }
  
}


##MF

MF.sams<-row.names(meta[meta$Category=="MF",])
##物种丰度数据
for(i in 2:5){
  gdmTab0 <- formatsitepair(bioData=abundT[MF.sams,],
                            bioFormat=1,
                            predData=meta_filled[MF.sams,list.par[[i]]],
                            XColumn="Lon", 
                            YColumn="Lat",
                            siteColumn="sites",abundance = T)
  
  list.gdb.metabo[[i+119]] <- gdm(data=gdmTab0, geo=TRUE)
  list.gdb.metabo[[i+119]]$dataname<-paste("MAG.abund.",pars[i],sep="")
  devia[119+i,1]<-"MAG-abund"
  devia[119+i,2]<-pars[i]
  devia[119+i,3]<-list.gdb.metabo[[i+119]]$explained
  devia[119+i,4]<-119+i
  devia[119+i,5]<-"MF"
}


for(i in 2:5){
  gdmTab0 <- formatsitepair(bioData=abundT[MF.sams,],
                            bioFormat=1,
                            predData=meta_filled[MF.sams,list.par[[i]]],
                            XColumn="Lon", 
                            YColumn="Lat",
                            siteColumn="sites")
  
  list.gdb.metabo[[i+123]] <- gdm(data=gdmTab0, geo=TRUE)
  list.gdb.metabo[[i+123]]$dataname<-paste("MAG-presence.",pars[i],sep="")
  devia[123+i,1]<-"MAG-presence"
  devia[123+i,2]<-pars[i]
  devia[123+i,3]<-list.gdb.metabo[[i+123]]$explained
  devia[123+i,4]<-123+i
  devia[123+i,5]<-"MF"
}



##功能丰度数据

for(i in 2:5){
  gdmTab0 <- formatsitepair(bioData=MAG.fun[MF.sams,],
                            bioFormat=1,
                            predData=meta_filled[MF.sams,list.par[[i]]],
                            XColumn="Lon", 
                            YColumn="Lat",
                            siteColumn="sites",abundance = T)
  
  list.gdb.metabo[[i+127]] <- gdm(data=gdmTab0, geo=TRUE)
  list.gdb.metabo[[i+127]]$dataname<-paste("Fun-abund.",pars[i],sep="")
  devia[127+i,1]<-"Fun-abund"
  devia[127+i,2]<-pars[i]
  devia[127+i,3]<-list.gdb.metabo[[i+127]]$explained
  devia[127+i,4]<-127+i
  devia[127+i,5]<-"MF"
}

n=133
for(i in 1:12){
  aa<-c("sites",metabolisms[metabolisms$Level.2==metalevel[i],]$KID)
  aa<-aa[is.element(aa,colnames(MAG.fun))]
  for(j in 2:5){
    if(length(aa)>2){
      MAG.fun0<-MAG.fun[MF.sams,aa]
      gdmTab0 <- formatsitepair(bioData=MAG.fun0,
                                bioFormat=1,
                                predData=meta_filled[MF.sams,list.par[[j]]],
                                XColumn="Lon", 
                                YColumn="Lat",
                                siteColumn="sites",abundance = T)
      
      list.gdb.metabo[[n]] <- gdm(data=gdmTab0, geo=TRUE)
      list.gdb.metabo[[n]]$dataname<-paste(metalevel[i],pars[i],sep="")
      devia[n,1]<-metalevel[i]
      devia[n,2]<-pars[j]
      devia[n,3]<-list.gdb.metabo[[n]]$explained
      devia[n,4]<-n
      devia[n,5]<-"MF"
    } else {
      list.gdb.metabo[[n]]<-paste(metalevel[i],pars[i],sep="")
      devia[n,1]<-metalevel[i]
      devia[n,2]<-pars[j]
      devia[n,3]<-NA
      devia[n,4]<-n
      devia[n,5]<-"MF"
    }
    n=n+1  
  }
  
}

##SG

SG.sams<-row.names(meta[meta$Category=="SG",])
##物种丰度数据
for(i in 2:5){
  gdmTab0 <- formatsitepair(bioData=abundT[SG.sams,],
                            bioFormat=1,
                            predData=meta_filled[SG.sams,list.par[[i]]],
                            XColumn="Lon", 
                            YColumn="Lat",
                            siteColumn="sites",abundance = T)
  
  list.gdb.metabo[[i+179]] <- gdm(data=gdmTab0, geo=TRUE)
  list.gdb.metabo[[i+179]]$dataname<-paste("MAG.abund.",pars[i],sep="")
  devia[179+i,1]<-"MAG-abund"
  devia[179+i,2]<-pars[i]
  devia[179+i,3]<-list.gdb.metabo[[i+179]]$explained
  devia[179+i,4]<-179+i
  devia[179+i,5]<-"SG"
}


for(i in 2:5){
  gdmTab0 <- formatsitepair(bioData=abundT[SG.sams,],
                            bioFormat=1,
                            predData=meta_filled[SG.sams,list.par[[i]]],
                            XColumn="Lon", 
                            YColumn="Lat",
                            siteColumn="sites")
  
  list.gdb.metabo[[i+183]] <- gdm(data=gdmTab0, geo=TRUE)
  list.gdb.metabo[[i+183]]$dataname<-paste("MAG-presence.",pars[i],sep="")
  devia[183+i,1]<-"MAG-presence"
  devia[183+i,2]<-pars[i]
  devia[183+i,3]<-list.gdb.metabo[[i+183]]$explained
  devia[183+i,4]<-183+i
  devia[183+i,5]<-"SG"
}



##功能丰度数据

for(i in 2:5){
  gdmTab0 <- formatsitepair(bioData=MAG.fun[SG.sams,],
                            bioFormat=1,
                            predData=meta_filled[SG.sams,list.par[[i]]],
                            XColumn="Lon", 
                            YColumn="Lat",
                            siteColumn="sites",abundance = T)
  
  list.gdb.metabo[[i+187]] <- gdm(data=gdmTab0, geo=TRUE)
  list.gdb.metabo[[i+187]]$dataname<-paste("Fun-abund.",pars[i],sep="")
  devia[187+i,1]<-"Fun-abund"
  devia[187+i,2]<-pars[i]
  devia[187+i,3]<-list.gdb.metabo[[i+187]]$explained
  devia[187+i,4]<-187+i
  devia[187+i,5]<-"SG"
}

n=193
for(i in 1:12){
  aa<-c("sites",metabolisms[metabolisms$Level.2==metalevel[i],]$KID)
  aa<-aa[is.element(aa,colnames(MAG.fun))]
  for(j in 2:5){
    if(length(aa)>2){
      MAG.fun0<-MAG.fun[SG.sams,aa]
      gdmTab0 <- formatsitepair(bioData=MAG.fun0,
                                bioFormat=1,
                                predData=meta_filled[SG.sams,list.par[[j]]],
                                XColumn="Lon", 
                                YColumn="Lat",
                                siteColumn="sites",abundance = T)
      
      list.gdb.metabo[[n]] <- gdm(data=gdmTab0, geo=TRUE)
      list.gdb.metabo[[n]]$dataname<-paste(metalevel[i],pars[i],sep="")
      devia[n,1]<-metalevel[i]
      devia[n,2]<-pars[j]
      devia[n,3]<-list.gdb.metabo[[n]]$explained
      devia[n,4]<-n
      devia[n,5]<-"SG"
    } else {
      list.gdb.metabo[[n]]<-paste(metalevel[i],pars[i],sep="")
      devia[n,1]<-metalevel[i]
      devia[n,2]<-pars[j]
      devia[n,3]<-NA
      devia[n,4]<-n
      devia[n,5]<-"SG"
    }
    n=n+1  
  }
  
}

lev<-unique(devia$Micro)

lev<-lev[c(13,14,15,1:12)]

library(RColorBrewer)
colors<-brewer.pal(8,"Dark2")
co2<-brewer.pal(9,"Set1")
co3<-brewer.pal(12,"Paired")
colors<-c(colors,co2,co3)
colors<-unique(colors)[1:15]

names(colors)<-lev

MAG-abund                                MAG-presence 
"#1B9E77"                                   "#D95F02" 
Fun-abund Biosynthesis of other secondary metabolites 
"#7570B3"                                   "#E7298A" 
Energy metabolism          Glycan biosynthesis and metabolism 
"#66A61E"                                   "#E6AB02" 
Lipid metabolism             Metabolism of other amino acids 
"#A6761D"                                   "#666666" 
Amino acid metabolism                     Carbohydrate metabolism 
"#E41A1C"                                   "#377EB8" 
Global and overview maps        Metabolism of cofactors and vitamins 
"#4DAF4A"                                   "#984EA3" 
Metabolism of terpenoids and polyketides                       Nucleotide metabolism 
"#FF7F00"                                   "#FFFF33" 
Xenobiotics biodegradation and metabolism 
"#A65628" 

devia$Micro<-factor(devia$Micro,levels = rev(lev))
devia$Cate<-factor(devia$Cate,levels = c("all","MF","SM","SG"))
devia$Predict<-factor(devia$Predict,levels = c("soil","carbon","climate","human"))

library(ggplot2)

ggplot(devia,aes(fill=Micro,x=Micro,y=Variance.explained))+
  geom_bar(position = "stack",stat = "identity",width = 0.8)+
  theme_bw()+
  theme( panel.grid.minor = element_blank(),panel.grid.major.y = element_blank(),
         axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black"))+
  ylab("Variance explained (%)")+labs(fill="Category")+
  scale_fill_manual(values=colors)+coord_flip()+
  facet_wrap(~ Cate+Predict, nrow = 4)
#export 20X20

summary(list.gdb.metabo[[50]])

row.names(devia)<-paste("gdm_",devia$Index,sep="")

gdm.1.splineDat <- isplineExtract(gdm.1)

relimpor<-NULL

for (i in 1:60){
  CC<-NULL
  for(j in c(4*i-2,4*i-3,4*i-1,4*i)){
  aa <- isplineExtract(list.gdb.metabo[[j]])
  bb <- data.frame(apply(aa$y,2,max,na.action=NULL))
  bb$Factor<-row.names(bb)
  colnames(bb)[1]<-"Importance"
  CC<-rbind(CC,bb)
  CC<-CC[!duplicated(CC$Factor),]
  }
  CC$Micro<-devia[4*i,]$Micro
  CC$Cate<-devia[4*i,]$Cate
  relimpor<-rbind(relimpor,CC)
}


level.fa<-c("Geographic",soil,carbon,climate,human)
level.fa<-unique(level.fa)

relimpor$Factor<-factor(relimpor$Factor,levels = level.fa)
  
MAG.aa<-relimpor[relimpor$Micro=="MAG-abund"|relimpor$Micro=="MAG-presence",]

MAG.aa[MAG.aa$Micro=="MAG-presence",]$Importance<-MAG.aa[MAG.aa$Micro=="MAG-presence",]$Importance*-1

p1<-ggplot(MAG.aa,aes(fill=Micro,x=Factor,y=Importance))+
  geom_bar(position = "stack",stat = "identity",width = 0.8)+
  theme_bw()+
  theme( panel.grid.minor = element_blank(),panel.grid.major.y = element_blank(),
        axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("relative importance")+labs(fill="Category")+
  scale_fill_manual(values=c('MAG-abund'="#E52925",'MAG-presence'="#FF6A72"))+
  facet_wrap(~ Cate, nrow = 4)

MAG.bb<-relimpor[relimpor$Micro=="MAG-abund"|relimpor$Micro=="Fun-abund",]

MAG.bb[MAG.bb$Micro=="Fun-abund",]$Importance<-MAG.bb[MAG.bb$Micro=="Fun-abund",]$Importance*-1

p2<-ggplot(MAG.bb,aes(fill=Micro,x=Factor,y=Importance))+
  geom_bar(position = "stack",stat = "identity",width = 0.8)+
  theme_bw()+
  theme( panel.grid.minor = element_blank(),panel.grid.major.y = element_blank(),
         axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("relative importance")+labs(fill="Category")+
  scale_fill_manual(values=c('MAG-abund'="#005BA3",'Fun-abund'="#0070A3"))+
  facet_wrap(~ Cate, nrow = 4)

ggpubr::ggarrange(p1,p2,ncol=2)
#export 10x20

relimpor$Micro<-factor(relimpor$Micro,levels = rev(levels(relimpor$Micro)))

relimpor2<-relimpor[is.element(relimpor$Micro,levels(relimpor$Micro)[4:15]),]

ggplot(relimpor2,aes(x=Factor,y=Importance))+
  geom_bar(position = "stack",stat = "identity",fill="#E52925",width = 0.8)+
  theme_bw()+
  theme( panel.grid.minor = element_blank(),panel.grid.major.y = element_blank(),
         axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("relative importance")+labs(fill="Category")+
  facet_wrap(~ Micro+Cate, ncol = 4)

##export 20X20

# PH,EC,CN,P of soil,MAG-abund,FUN-abund,MAG-presence

magAl<-c("MAG-abund","Fun-abund","MAG-presence")

ph.ec<-devia[devia$Predict=="soil"&is.element(devia$Micro,magAl),]

ph.da<-NULL

for(i in 1:12){
  bb<-data.frame(matrix(nrow=200,ncol=5))
  colnames(bb)<-c("x","y","Cate","Micro","Factor")
  aa<-isplineExtract(list.gdb.metabo[[ph.ec[i,]$Index]])
    bb$x<-aa$x[,"pH"]
    bb$y<-aa$y[,"pH"]
    bb$Cate<-ph.ec[i,]$Cate
    bb$Micro<-ph.ec[i,]$Micro
    bb$Factor<-"pH"
    ph.da<-rbind(ph.da,bb)
  }
  

ec.da<-NULL

for(i in 1:12){
  bb<-data.frame(matrix(nrow=200,ncol=5))
  colnames(bb)<-c("x","y","Cate","Micro","Factor")
  aa<-isplineExtract(list.gdb.metabo[[ph.ec[i,]$Index]])
  bb$x<-aa$x[,"EC"]
  bb$y<-aa$y[,"EC"]
  bb$Cate<-ph.ec[i,]$Cate
  bb$Micro<-ph.ec[i,]$Micro
  bb$Factor<-"EC"
  ec.da<-rbind(ec.da,bb)
}

cn.da<-NULL

for(i in 1:12){
  bb<-data.frame(matrix(nrow=200,ncol=5))
  colnames(bb)<-c("x","y","Cate","Micro","Factor")
  aa<-isplineExtract(list.gdb.metabo[[ph.ec[i,]$Index]])
  bb$x<-aa$x[,"C_N"]
  bb$y<-aa$y[,"C_N"]
  bb$Cate<-ph.ec[i,]$Cate
  bb$Micro<-ph.ec[i,]$Micro
  bb$Factor<-"C_N"
  cn.da<-rbind(cn.da,bb)
}


p.da<-NULL

for(i in 1:12){
  bb<-data.frame(matrix(nrow=200,ncol=5))
  colnames(bb)<-c("x","y","Cate","Micro","Factor")
  aa<-isplineExtract(list.gdb.metabo[[ph.ec[i,]$Index]])
  bb$x<-aa$x[,"P"]
  bb$y<-aa$y[,"P"]
  bb$Cate<-ph.ec[i,]$Cate
  bb$Micro<-ph.ec[i,]$Micro
  bb$Factor<-"P"
  p.da<-rbind(p.da,bb)
}

carbon.mod<-devia[devia$Predict=="carbon"&is.element(devia$Micro,magAl),]

caoc.da<-NULL

for(i in 1:12){
  bb<-data.frame(matrix(nrow=200,ncol=5))
  colnames(bb)<-c("x","y","Cate","Micro","Factor")
  aa<-isplineExtract(list.gdb.metabo[[carbon.mod[i,]$Index]])
  bb$x<-aa$x[,"Ca.OC"]
  bb$y<-aa$y[,"Ca.OC"]
  bb$Cate<-carbon.mod[i,]$Cate
  bb$Micro<-carbon.mod[i,]$Micro
  bb$Factor<-"Ca.OC"
  caoc.da<-rbind(caoc.da,bb)
}

doc.da<-NULL

for(i in 1:12){
  bb<-data.frame(matrix(nrow=200,ncol=5))
  colnames(bb)<-c("x","y","Cate","Micro","Factor")
  aa<-isplineExtract(list.gdb.metabo[[carbon.mod[i,]$Index]])
  bb$x<-aa$x[,"DOC"]
  bb$y<-aa$y[,"DOC"]
  bb$Cate<-carbon.mod[i,]$Cate
  bb$Micro<-carbon.mod[i,]$Micro
  bb$Factor<-"DOC"
  doc.da<-rbind(doc.da,bb)
}


climate.mod<-devia[devia$Predict=="climate"&is.element(devia$Micro,magAl),]

mat.da<-NULL

for(i in 1:12){
  bb<-data.frame(matrix(nrow=200,ncol=5))
  colnames(bb)<-c("x","y","Cate","Micro","Factor")
  aa<-isplineExtract(list.gdb.metabo[[climate.mod[i,]$Index]])
  bb$x<-aa$x[,"MAT"]
  bb$y<-aa$y[,"MAT"]
  bb$Cate<-climate.mod[i,]$Cate
  bb$Micro<-climate.mod[i,]$Micro
  bb$Factor<-"MAT"
  mat.da<-rbind(mat.da,bb)
}


map.da<-NULL

for(i in 1:12){
  bb<-data.frame(matrix(nrow=200,ncol=5))
  colnames(bb)<-c("x","y","Cate","Micro","Factor")
  aa<-isplineExtract(list.gdb.metabo[[climate.mod[i,]$Index]])
  bb$x<-aa$x[,"MAP"]
  bb$y<-aa$y[,"MAP"]
  bb$Cate<-climate.mod[i,]$Cate
  bb$Micro<-climate.mod[i,]$Micro
  bb$Factor<-"MAP"
  map.da<-rbind(map.da,bb)
}


human.mod<-devia[devia$Predict=="human"&is.element(devia$Micro,magAl),]

agri.da<-NULL

for(i in 1:12){
  bb<-data.frame(matrix(nrow=200,ncol=5))
  colnames(bb)<-c("x","y","Cate","Micro","Factor")
  aa<-isplineExtract(list.gdb.metabo[[human.mod[i,]$Index]])
  bb$x<-aa$x[,"Agriculture"]
  bb$y<-aa$y[,"Agriculture"]
  bb$Cate<-human.mod[i,]$Cate
  bb$Micro<-human.mod[i,]$Micro
  bb$Factor<-"Agriculture"
  agri.da<-rbind(agri.da,bb)
}


densi.da<-NULL

for(i in 1:12){
  bb<-data.frame(matrix(nrow=200,ncol=5))
  colnames(bb)<-c("x","y","Cate","Micro","Factor")
  aa<-isplineExtract(list.gdb.metabo[[human.mod[i,]$Index]])
  bb$x<-aa$x[,"Density"]
  bb$y<-aa$y[,"Density"]
  bb$Cate<-human.mod[i,]$Cate
  bb$Micro<-human.mod[i,]$Micro
  bb$Factor<-"Density"
  densi.da<-rbind(densi.da,bb)
}

max(ph.da$x)#9.13
min(ph.da$x)#4.04
median(ph.da$x)#7.013291
max(ec.da$x)#12.43
min(ec.da$x)#0.32
median(ec.da$x)#4.745025

a=(max(ph.da$x)-min(ph.da$x))/(max(ec.da$x)-min(ec.da$x))
b=min(ph.da$x)-a*min(ec.da$x)


p3<-ggplot()+
  #scale_color_manual(values =c("EC"="#005BA3","pH"="#E52925"))+
  geom_line(data=ph.da,aes(x=x,y=y,linetype = Micro),color="#E52925",size=1)+
  geom_line(data=ec.da,aes(x=(x*a+b),y=y,linetype = Micro),color="#005BA3",size=1)+
  scale_x_continuous(name = "pH",#y1特征
                     sec.axis = sec_axis(trans=~((. - b)/a), name="EC"))+
  theme_bw()+
  theme( panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
         axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black"))+
  #theme(axis.text.x = element_text(angle = 45, hjust = 2)) +
  ylab("Turnover")+labs(fill="Category")+guides(color = guide_legend(order=2))+
  ggtitle("red PH, blue EC")+
  facet_wrap(~ Cate, ncol = 4)

max(cn.da$x)#70
min(cn.da$x)#5
max(p.da$x)#1498
min(p.da$x)#90

c=(max(p.da$x)-min(p.da$x))/(max(cn.da$x)-min(cn.da$x))
d=min(p.da$x)-c*min(cn.da$x)


p4<-ggplot()+
  #scale_color_manual(values =c("EC"="#005BA3","pH"="#E52925"))+
  geom_line(data=p.da,aes(x=x,y=y,linetype = Micro),color="#E52925",size=1)+
  geom_line(data=cn.da,aes(x=(x*c+d),y=y,linetype = Micro),color="#005BA3",size=1)+
  scale_x_continuous(name = "P",#y1特征
                     sec.axis = sec_axis(trans=~((. - d)/c), name="C_N"))+
  theme_bw()+
  theme( panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
         axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black"))+
  #theme(axis.text.x = element_text(angle = 45, hjust = 2)) +
  ylab("Turnover")+labs(fill="Category")+guides(color = guide_legend(order=2))+
  ggtitle("red P, blue C_N")+
  facet_wrap(~ Cate, ncol = 4)



max(caoc.da$x)#1.19
min(caoc.da$x)#0.04
max(doc.da$x)#137.
min(doc.da$x)#3.5

e=(max(doc.da$x)-min(doc.da$x))/(max(caoc.da$x)-min(caoc.da$x))
f=min(doc.da$x)-c*min(caoc.da$x)


p5<-ggplot()+
  #scale_color_manual(values =c("EC"="#005BA3","pH"="#E52925"))+
  geom_line(data=doc.da,aes(x=x,y=y,linetype = Micro),color="#E52925",size=1)+
  geom_line(data=caoc.da,aes(x=(x*e+f),y=y,linetype = Micro),color="#005BA3",size=1)+
  scale_x_continuous(name = "DOC",#y1特征
                     sec.axis = sec_axis(trans=~((. - f)/e), name="Ca.OC"))+
  theme_bw()+
  theme( panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
         axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black"))+
  #theme(axis.text.x = element_text(angle = 45, hjust = 2)) +
  ylab("Turnover")+labs(fill="Category")+guides(color = guide_legend(order=2))+
  ggtitle("red DOC, blue Ca.OC")+
  facet_wrap(~ Cate, ncol = 4)

max(map.da$x)#2180.1
min(map.da$x)#560.65
max(mat.da$x)#24.34
min(mat.da$x)#8.89


g=(max(map.da$x)-min(map.da$x))/(max(mat.da$x)-min(mat.da$x))
h=min(map.da$x)-c*min(mat.da$x)


p6<-ggplot()+
  #scale_color_manual(values =c("EC"="#005BA3","pH"="#E52925"))+
  geom_line(data=map.da,aes(x=x,y=y,linetype = Micro),color="#E52925",size=1)+
  geom_line(data=mat.da,aes(x=(x*e+f),y=y,linetype = Micro),color="#005BA3",size=1)+
  scale_x_continuous(name = "MAP",#y1特征
                     sec.axis = sec_axis(trans=~((. - h)/g), name="MAT"))+
  theme_bw()+
  theme( panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
         axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black"))+
  #theme(axis.text.x = element_text(angle = 45, hjust = 2)) +
  ylab("Turnover")+labs(fill="Category")+guides(color = guide_legend(order=2))+
  ggtitle("red MAP, blue MAT")+
  facet_wrap(~ Cate, ncol = 4)



max(agri.da$x)#74746
min(agri.da$x)#13
max(densi.da$x)#9833
min(densi.da$x)#0


m=(max(agri.da$x)-min(agri.da$x))/(max(densi.da$x)-min(densi.da$x))
n=min(agri.da$x)-c*min(densi.da$x)


p7<-ggplot()+
  #scale_color_manual(values =c("EC"="#005BA3","pH"="#E52925"))+
  geom_line(data=agri.da,aes(x=x,y=y,linetype = Micro),color="#E52925",size=1)+
  geom_line(data=densi.da,aes(x=(x*e+f),y=y,linetype = Micro),color="#005BA3",size=1)+
  scale_x_continuous(name = "Agricultrue",#y1特征
                     sec.axis = sec_axis(trans=~((. - h)/g), name="Density"))+
  theme_bw()+
  theme( panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
         axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black"))+
  #theme(axis.text.x = element_text(angle = 45, hjust = 2)) +
  ylab("Turnover")+labs(fill="Category")+guides(color = guide_legend(order=2))+
  ggtitle("red Agriculture, blue Density")+
  facet_wrap(~ Cate, ncol = 4)

ggpubr::ggarrange(p3,p4,p5,p6,p7,ncol=1)
#export 20*20

relimpor3<-relimpor[(relimpor$Micro=="MAG-abund"|relimpor$Micro=="Fun-abund")&relimpor$Cate=="all",]
write.table(relimpor3,"Gdm-importance.txt",row.names = F,sep="\t",quote=F)


## spline selected


# CN,P of soil,MAG-abund,FUN-abund

magA2<-c("MAG-abund","Fun-abund")

cn.p<-devia[devia$Predict=="soil"&is.element(devia$Micro,magA2)&devia$Cate=="all",]

cn.da2<-NULL

for(i in 1:2){
  bb<-data.frame(matrix(nrow=200,ncol=5))
  colnames(bb)<-c("x","y","Cate","Micro","Factor")
  aa<-isplineExtract(list.gdb.metabo[[cn.p[i,]$Index]])
  bb$x<-aa$x[,"C_N"]
  bb$y<-aa$y[,"C_N"]
  bb$Cate<-cn.p[i,]$Cate
  bb$Micro<-cn.p[i,]$Micro
  bb$Factor<-"C_N"
  cn.da2<-rbind(cn.da2,bb)
}
cn.da2$Micro<-factor(cn.da2$Micro,levels = c("MAG-abund","Fun-abund"))

p.da2<-NULL

for(i in 1:2){
  bb<-data.frame(matrix(nrow=200,ncol=5))
  colnames(bb)<-c("x","y","Cate","Micro","Factor")
  aa<-isplineExtract(list.gdb.metabo[[cn.p[i,]$Index]])
  bb$x<-aa$x[,"P"]
  bb$y<-aa$y[,"P"]
  bb$Cate<-cn.p[i,]$Cate
  bb$Micro<-cn.p[i,]$Micro
  bb$Factor<-"P"
  p.da2<-rbind(p.da2,bb)
}

p.da2$Micro<-factor(p.da2$Micro,levels = c("MAG-abund","Fun-abund"))
##DOC MAOC
carbon.mod2<-devia[devia$Cate=="all"&devia$Predict=="carbon"&is.element(devia$Micro,magA2),]

moc.da<-NULL

for(i in 1:2){
  bb<-data.frame(matrix(nrow=200,ncol=5))
  colnames(bb)<-c("x","y","Cate","Micro","Factor")
  aa<-isplineExtract(list.gdb.metabo[[carbon.mod2[i,]$Index]])
  bb$x<-aa$x[,"MOC"]
  bb$y<-aa$y[,"MOC"]
  bb$Cate<-carbon.mod2[i,]$Cate
  bb$Micro<-carbon.mod2[i,]$Micro
  bb$Factor<-"MOC"
  moc.da<-rbind(moc.da,bb)
}
moc.da$Micro<-factor(moc.da$Micro,levels = c("MAG-abund","Fun-abund"))

doc.da2<-NULL

for(i in 1:2){
  bb<-data.frame(matrix(nrow=200,ncol=5))
  colnames(bb)<-c("x","y","Cate","Micro","Factor")
  aa<-isplineExtract(list.gdb.metabo[[carbon.mod2[i,]$Index]])
  bb$x<-aa$x[,"DOC"]
  bb$y<-aa$y[,"DOC"]
  bb$Cate<-carbon.mod2[i,]$Cate
  bb$Micro<-carbon.mod2[i,]$Micro
  bb$Factor<-"DOC"
  doc.da2<-rbind(doc.da2,bb)
}
doc.da2$Micro<-factor(doc.da2$Micro,levels = c("MAG-abund","Fun-abund"))

## map mat
climate.mod2<-devia[devia$Cate=="all"&devia$Predict=="climate"&is.element(devia$Micro,magA2),]

mat.da2<-NULL

for(i in 1:2){
  bb<-data.frame(matrix(nrow=200,ncol=5))
  colnames(bb)<-c("x","y","Cate","Micro","Factor")
  aa<-isplineExtract(list.gdb.metabo[[climate.mod2[i,]$Index]])
  bb$x<-aa$x[,"MAT"]
  bb$y<-aa$y[,"MAT"]
  bb$Cate<-climate.mod2[i,]$Cate
  bb$Micro<-climate.mod2[i,]$Micro
  bb$Factor<-"MAT"
  mat.da2<-rbind(mat.da2,bb)
}
mat.da2$Micro<-factor(mat.da2$Micro,levels = c("MAG-abund","Fun-abund"))

map.da2<-NULL

for(i in 1:2){
  bb<-data.frame(matrix(nrow=200,ncol=5))
  colnames(bb)<-c("x","y","Cate","Micro","Factor")
  aa<-isplineExtract(list.gdb.metabo[[climate.mod2[i,]$Index]])
  bb$x<-aa$x[,"MAP"]
  bb$y<-aa$y[,"MAP"]
  bb$Cate<-climate.mod2[i,]$Cate
  bb$Micro<-climate.mod2[i,]$Micro
  bb$Factor<-"MAP"
  map.da2<-rbind(map.da2,bb)
}
map.da2$Micro<-factor(map.da2$Micro,levels = c("MAG-abund","Fun-abund"))

human.mod2<-devia[devia$Cate=="all"&devia$Predict=="human"&is.element(devia$Micro,magA2),]

agri.da2<-NULL

for(i in 1:2){
  bb<-data.frame(matrix(nrow=200,ncol=5))
  colnames(bb)<-c("x","y","Cate","Micro","Factor")
  aa<-isplineExtract(list.gdb.metabo[[human.mod2[i,]$Index]])
  bb$x<-aa$x[,"Agriculture"]
  bb$y<-aa$y[,"Agriculture"]
  bb$Cate<-human.mod2[i,]$Cate
  bb$Micro<-human.mod2[i,]$Micro
  bb$Factor<-"Agriculture"
  agri.da2<-rbind(agri.da2,bb)
}
agri.da2$Micro<-factor(agri.da2$Micro,levels = c("MAG-abund","Fun-abund"))

densi.da2<-NULL

for(i in 1:2){
  bb<-data.frame(matrix(nrow=200,ncol=5))
  colnames(bb)<-c("x","y","Cate","Micro","Factor")
  aa<-isplineExtract(list.gdb.metabo[[human.mod2[i,]$Index]])
  bb$x<-aa$x[,"Density"]
  bb$y<-aa$y[,"Density"]
  bb$Cate<-human.mod2[i,]$Cate
  bb$Micro<-human.mod2[i,]$Micro
  bb$Factor<-"Density"
  densi.da2<-rbind(densi.da2,bb)
}
densi.da2$Micro<-factor(densi.da2$Micro,levels = c("MAG-abund","Fun-abund"))


max(cn.da2$x)#70
min(cn.da2$x)#5
max(p.da2$x)#1498
min(p.da2$x)#90

c2=(max(p.da2$x)-min(p.da2$x))/(max(cn.da2$x)-min(cn.da2$x))
d2=min(p.da2$x)-c2*min(cn.da2$x)


p8<-ggplot()+
  #scale_color_manual(values =c("EC"="#005BA3","pH"="#E52925"))+
  geom_line(data=p.da2,aes(x=x,y=y,linetype = Micro),color="#E52925",size=1)+
  geom_line(data=cn.da2,aes(x=(x*c2+d2),y=y,linetype = Micro),color="#005BA3",size=1)+
  scale_x_continuous(name = "P",#y1特征
                     sec.axis = sec_axis(trans=~((. - d2)/c2), name="C_N"))+
  theme_bw()+
  theme( panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
         axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black"))+
  #theme(axis.text.x = element_text(angle = 45, hjust = 2)) +
  ylab("Turnover")+labs(fill="Category")+guides(color = guide_legend(order=2))+
  ggtitle("red P, blue C_N")



max(moc.da$x)#1.19
min(moc.da$x)#0.04
max(doc.da2$x)#137.
min(doc.da2$x)#3.5

e2=(max(doc.da2$x)-min(doc.da2$x))/(max(moc.da$x)-min(moc.da$x))
f2=min(doc.da2$x)-e2*min(moc.da$x)


p9<-ggplot()+
  #scale_color_manual(values =c("EC"="#005BA3","pH"="#E52925"))+
  geom_line(data=doc.da2,aes(x=x,y=y,linetype = Micro),color="#E52925",size=1)+
  geom_line(data=moc.da,aes(x=(x*e2+f2),y=y,linetype = Micro),color="#005BA3",size=1)+
  scale_x_continuous(name = "DOC",#y1特征
                     sec.axis = sec_axis(trans=~((. - f2)/e2), name="MOC"))+
  theme_bw()+
  theme( panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
         axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black"))+
  #theme(axis.text.x = element_text(angle = 45, hjust = 2)) +
  ylab("Turnover")+labs(fill="Category")+guides(color = guide_legend(order=2))+
  ggtitle("red DOC, blue MOC")

max(map.da2$x)#2180.1
min(map.da2$x)#560.65
max(mat.da2$x)#24.34
min(mat.da2$x)#8.89


g2=(max(map.da2$x)-min(map.da2$x))/(max(mat.da2$x)-min(mat.da2$x))
h2=min(map.da2$x)-g2*min(mat.da2$x)


p10<-ggplot()+
  #scale_color_manual(values =c("EC"="#005BA3","pH"="#E52925"))+
  geom_line(data=map.da2,aes(x=x,y=y,linetype = Micro),color="#E52925",size=1)+
  geom_line(data=mat.da2,aes(x=(x*g2+h2),y=y,linetype = Micro),color="#005BA3",size=1)+
  scale_x_continuous(name = "MAP",#y1特征
                     sec.axis = sec_axis(trans=~((. - h2)/g2), name="MAT"))+
  theme_bw()+
  theme( panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
         axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black"))+
  #theme(axis.text.x = element_text(angle = 45, hjust = 2)) +
  ylab("Turnover")+labs(fill="Category")+guides(color = guide_legend(order=2))+
  ggtitle("red MAP, blue MAT")



max(agri.da2$x)#74746
min(agri.da2$x)#13
max(densi.da2$x)#9833
min(densi.da2$x)#0


m2=(max(agri.da2$x)-min(agri.da2$x))/(max(densi.da2$x)-min(densi.da2$x))
n2=min(agri.da2$x)-m2*min(densi.da2$x)


p11<-ggplot()+
  #scale_color_manual(values =c("EC"="#005BA3","pH"="#E52925"))+
  geom_line(data=agri.da2,aes(x=x,y=y,linetype = Micro),color="#E52925",size=1)+
  geom_line(data=densi.da2,aes(x=(x*m2+n2),y=y,linetype = Micro),color="#005BA3",size=1)+
  scale_x_continuous(name = "Agricultrue",#y1特征
                     sec.axis = sec_axis(trans=~((. - n2)/m2), name="Density"))+
  theme_bw()+
  theme( panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
         axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black"))+
  #theme(axis.text.x = element_text(angle = 45, hjust = 2)) +
  ylab("Turnover")+labs(fill="Category")+guides(color = guide_legend(order=2))+
  ggtitle("red Agriculture, blue Density")

ggpubr::ggarrange(p8,p9,p10,p11,ncol=1)
#export 8*15


devia3<-devia[devia$Cate=="all"&is.element(devia$Micro,magA2),]

geo.mag.da <- formatsitepair(bioData=abundT,
                          bioFormat=1,
                          predData=meta_filled[,c("sites","Lat","Lon")],
                          XColumn="Lon", 
                          YColumn="Lat",
                          siteColumn="sites",abundance = T)

geo.mag.da.mod <- gdm(data=geo.mag.da, geo=TRUE)
geo.mag.da.mod$dataname<-"MAG.abund.geo"
geo.mag.da.mod$explained

geo.fun.da <- formatsitepair(bioData=MAG.fun,
                             bioFormat=1,
                             predData=meta_filled[,c("sites","Lat","Lon")],
                             XColumn="Lon", 
                             YColumn="Lat",
                             siteColumn="sites",abundance = T)

geo.fun.da.mod <- gdm(data=geo.fun.da, geo=TRUE)
geo.fun.da.mod$dataname<-"MAG.fun.geo"
geo.fun.da.mod$explained

devia3[9,]$Variance.explained<-geo.mag.da.mod$explained
devia3[10,]$Variance.explained<-geo.fun.da.mod$explained

devia3[9,c(1,2,4,5)]<-c("MAG-abund","Geographic","0","all")
devia3[10,c(1,2,4,5)]<-c("Fun-abund","Geographic","0","all")

devia3$Predict<-as.character(devia3$Predict)
devia3[9:10,]$Predict<-"Geograph"

devia3$Predict<-factor(devia3$Predict,levels=c("Geograph","human","climate","carbon","soil"))

ggplot(devia3,aes(fill=Micro,x=Predict,y=Variance.explained))+
  geom_bar(position = "dodge",stat = "identity",width = 0.8)+
  theme_bw()+
  theme( panel.grid.minor = element_blank(),panel.grid.major.y = element_blank(),
         axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black"))+
  ylab("Variance explained (%)")+labs(fill="Category")+
  scale_fill_manual(values=c("MAG-abund"="#E52925","Fun-abund"="#005BA3"))+coord_flip()

write.table(devia3,"variance2.txt",sep="\t",quote = F)

##access p value
# devtools::install_github('skiptoniam/bbgdm')
# library(bbgdm)
# 
# bbgdm.wald.test(gdm.1,gdm=T)


source("gdm.test.R")























