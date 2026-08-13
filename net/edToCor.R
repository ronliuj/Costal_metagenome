#' # The algorithm creat the corelation matrix from edges exported from CoNet for the analysis of ggClusterNet.
#'
#' @title edges to cor
#' @description Enter edges exported from CoNet igraph, produce cor matrix for ggClusterNet
#' @param edeges edges dataframe exported from CoNet igraph
#' @details
#' @return cor matrix
#' @author Contact: Jian Liu \email{liujian@@dzu.edu.cn} 
#' @export
library(dplyr)
library(reshape2)
library(tidyr)
library(tidyverse)

edToCor<-function(edges){
  eds<-edges 
  eds$shared.name<-gsub("OTU-","",eds$shared.name)
  eds<-eds %>% separate(shared.name, c('OTU_1', 'OTU_2'),sep="->",remove=F)
  OTUs<-unique(c(eds$OTU_1,eds$OTU_2))
  cc<-data.frame(matrix(nrow=length(OTUs),ncol=length(OTUs)))
  row.names(cc)<-OTUs
  colnames(cc)<-OTUs
  for(i in 1:length(OTUs)){
    for(j in 1:length(OTUs)){
      cc[i,j]<-ifelse(i==j,1,eds[(eds$OTU_1==OTUs[i]&eds$OTU_2==OTUs[j])|(eds$OTU_2==OTUs[i]&eds$OTU_1==OTUs[j]),]$weight)
    }
  }
  cc[is.na(cc)]<-0
  cc<-as.matrix(cc)
  return(cc)
}
