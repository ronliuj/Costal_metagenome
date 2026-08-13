install.packages("tidyr")
install.packages("fmsb")
install.packages("reshape2")
install.packages("ggplot2")


library(ggplot2)
library(reshape2)
library(fmsb)
library(dplyr)
library(tidyr)


RADdata <- read.table("Blue carbon RAD116.txt", header=TRUE, sep="\t", quote="", fileEncoding="GBK")
RADdata <- read.table("Blue carbon116.all.txt", header=TRUE, sep="\t", quote="", fileEncoding="GBK")
RADdata<-RADdata[,-1]
Geop.Chem<-RADdata[c("Category","EC","TN", "pH", "C.N", "P", "C.P", "Fed", "Feo", "Fep","Ald","Alo","Alp")]
Org.C<-RADdata[c("Category","TOC","DOC", "POC", "MOC", "Fe.OC", "Ca.OC", "Lignin.phenol", "Microbial.necromass")]
Others<-RADdata[c("Category","Lon","Lat", "MAT", "MAP", "人口密度", "总种植面积")]
#Geop.Chem
RADdata_median <- as.data.frame(Geop.Chem %>% group_by(Category) %>% summarize_if(is.numeric, list("q25" = ~quantile(., 0.25,type=2),"q75" = ~quantile(., 0.75,type=2),"q50" = ~quantile(., 0.50,type=2))))
median_only <- as.data.frame(Geop.Chem %>% group_by(Category) %>% summarize_if(is.numeric, list("q50" = ~quantile(., 0.50,type=2))))
#Org.C
RADdata_median <- as.data.frame(Org.C %>% group_by(Category) %>% summarize_if(is.numeric, list("q25" = ~quantile(., 0.25,type=2,na.rm = TRUE),"q75" = ~quantile(., 0.75,type=2,na.rm = TRUE),"q50" = ~quantile(., 0.50,type=2,na.rm = TRUE))))
median_only <- as.data.frame(Org.C %>% group_by(Category) %>% summarize_if(is.numeric, list("q50" = ~quantile(., 0.50,type=2,na.rm = TRUE))))
#Others
RADdata_median <- as.data.frame(Others %>% group_by(Category) %>% summarize_if(is.numeric, list("q25" = ~quantile(., 0.25,type=2,na.rm = TRUE),"q75" = ~quantile(., 0.75,type=2,na.rm = TRUE),"q50" = ~quantile(., 0.50,type=2,na.rm = TRUE))))
median_only <- as.data.frame(Others %>% group_by(Category) %>% summarize_if(is.numeric, list("q50" = ~quantile(., 0.50,type=2,na.rm = TRUE))))

#Geop.Chem 13列
process_category_data <- function(data, category_name) {
  # 选择一个区域
  metadata_site <- subset(data, Category == category_name)
  melt_site <- melt(metadata_site)
  # 筛选q25, q50和q75数据
  q25_data <- subset(melt_site, grepl("_q25$", variable))
  q50_data <- subset(melt_site, grepl("_q50$", variable))
  q75_data <- subset(melt_site, grepl("_q75$", variable))
  
  # 数据透视
  q25_pivot <- pivot_wider(q25_data, names_from = variable, values_from = value)
  q50_pivot <- pivot_wider(q50_data, names_from = variable, values_from = value)
  q75_pivot <- pivot_wider(q75_data, names_from = variable, values_from = value)
  
  # 重命名列
  colnames(q25_pivot) <- c("group","EC","TN", "pH", "C.N", "P", "C.P", "Fed", "Feo", "Fep","Ald","Alo","Alp")
  colnames(q50_pivot) <- c("group","EC","TN", "pH", "C.N", "P", "C.P", "Fed", "Feo", "Fep","Ald","Alo","Alp")
  colnames(q75_pivot) <- c("group","EC","TN", "pH", "C.N", "P", "C.P", "Fed", "Feo", "Fep","Ald","Alo","Alp")
  
  # 合并q25, q50和q75数据
  min <- data.frame(matrix(rep(c(0), 13), nrow = 1))
  colnames(min)<-c("group","EC","TN", "pH", "C.N", "P", "C.P", "Fed", "Feo", "Fep","Ald","Alo","Alp")
  combined_data <- rbind(min, q75_pivot, q50_pivot)
  
  # 添加新列"group1"
  combined_data$group1 <- c("min", "max", "Q50")
  combined_data$group <- NULL
  
  combined_data <- combined_data %>%
    select(group1, everything())
  
  return(combined_data)
}
#Org.C 9列
process_category_data <- function(data, category_name) {
  # 选择一个区域
  metadata_site <- subset(data, Category == category_name)
  melt_site <- melt(metadata_site)
  # 筛选q25, q50和q75数据
  q25_data <- subset(melt_site, grepl("_q25$", variable))
  q50_data <- subset(melt_site, grepl("_q50$", variable))
  q75_data <- subset(melt_site, grepl("_q75$", variable))
  
  # 数据透视
  q25_pivot <- pivot_wider(q25_data, names_from = variable, values_from = value)
  q50_pivot <- pivot_wider(q50_data, names_from = variable, values_from = value)
  q75_pivot <- pivot_wider(q75_data, names_from = variable, values_from = value)
  
  # 重命名列
  colnames(q25_pivot) <- c("group","TOC","DOC", "POC", "MOC", "Fe.OC", "Ca.OC", "Lignin.phenol", "Microbial.necromass")
  colnames(q50_pivot) <- c("group","TOC","DOC", "POC", "MOC", "Fe.OC", "Ca.OC", "Lignin.phenol", "Microbial.necromass")
  colnames(q75_pivot) <- c("group","TOC","DOC", "POC", "MOC", "Fe.OC", "Ca.OC", "Lignin.phenol", "Microbial.necromass")
  
  # 合并q25, q50和q75数据
  min <- data.frame(matrix(rep(c(0), 9), nrow = 1))
  colnames(min)<-c("group","TOC","DOC", "POC", "MOC", "Fe.OC", "Ca.OC", "Lignin.phenol", "Microbial.necromass")
  combined_data <- rbind(min, q75_pivot, q50_pivot)
  
  # 添加新列"group1"
  combined_data$group1 <- c("min", "max", "Q50")
  combined_data$group <- NULL
  
  combined_data <- combined_data %>%
    select(group1, everything())
  
  return(combined_data)
}
#Others 7列
process_category_data <- function(data, category_name) {
  # 选择一个区域
  metadata_site <- subset(data, Category == category_name)
  melt_site <- melt(metadata_site)
  # 筛选q25, q50和q75数据
  q25_data <- subset(melt_site, grepl("_q25$", variable))
  q50_data <- subset(melt_site, grepl("_q50$", variable))
  q75_data <- subset(melt_site, grepl("_q75$", variable))
  
  # 数据透视
  q25_pivot <- pivot_wider(q25_data, names_from = variable, values_from = value)
  q50_pivot <- pivot_wider(q50_data, names_from = variable, values_from = value)
  q75_pivot <- pivot_wider(q75_data, names_from = variable, values_from = value)
  
  # 重命名列
  colnames(q25_pivot) <- c("group","Lon","Lat", "MAT", "MAP", "人口密度", "总种植面积")
  colnames(q50_pivot) <- c("group","Lon","Lat", "MAT", "MAP", "人口密度", "总种植面积")
  colnames(q75_pivot) <- c("group","Lon","Lat", "MAT", "MAP", "人口密度", "总种植面积")
  
  # 合并q25, q50和q75数据
  min <- data.frame(matrix(rep(c(0), 7), nrow = 1))
  colnames(min)<-c("group","Lon","Lat", "MAT", "MAP", "人口密度", "总种植面积")
  combined_data <- rbind(min, q75_pivot, q50_pivot)
  
  # 添加新列"group1"
  combined_data$group1 <- c("min", "max", "Q50")
  combined_data$group <- NULL
  
  combined_data <- combined_data %>%
    select(group1, everything())
  
  return(combined_data)
}

# 定义ggspider函数
ggspider <- function(data, axis_name_offset = 0.1, background_color = "white", fill_opacity = 0.5) {
  # 函数内容如上
}

## 生成和保存蜘蛛图的函数
generate_and_save_spider_plot <- function(data, category_name) {
  gg <- ggspider(data, axis_name_offset = 0.1, background_color = "white", fill_opacity = 0.4) +
    scale_fill_manual(values = c("Q25" = "black", "Q50" = "#3F459BFF", "Q75" = "white")) +
    scale_color_manual(values = c("Q25" = "black", "Q50" = "#3F459BFF", "Q75" = "#3F459BFF"))
  # 创建PDF设备
  pdf(paste0(category_name, "_spider_plot.pdf"), width = 8, height = 8)
  # 打印图表
  print(gg)
  # 关闭PDF设备
  dev.off()
}

## 计算各站点数据框
MF_data <- process_category_data(RADdata_median, "MF")
SM_data <- process_category_data(RADdata_median, "SM")
SG_data <- process_category_data(RADdata_median, "SG")

## Generate and save the plots
generate_and_save_spider_plot(MF_data, "MF")
generate_and_save_spider_plot(SM_data, "SM")
generate_and_save_spider_plot(SG_data, "SG")

## 将数据缩放到0-100范围的函数
scale_column_to_0_100 <- function(x) {
  # 直接计算，不使用 scale()
  scaled <- (x / max(x, na.rm = TRUE)) * 100
  # 确保返回向量
  return(as.vector(scaled))
}

## 对每个数值列应用缩放
for (col_name in colnames(median_only)[-1]) {
  median_only[, col_name] <- scale_column_to_0_100(median_only[, col_name])
}

median_only_rescaled <-median_only


#Geop.Chem
## 添加最大值和最小值行
min <- data.frame( 
  Category = "Min",
  EC_q50 = 0,
  pH_q50 = 0,
  TN_q50 = 0,
  P_q50 = 0,
  C.P_q50 = 0,
  C.N_q50 = 0,
  Feo_q50 = 0,
  Fep_q50=0,
  Fed_q50 = 0,
  Ald_q50 = 0,
  Alo_q50=0,
  Alp_q50=0
)
#Org.C
min <- data.frame( 
  Category = "Min",
  TOC_q50 = 0,
  DOC_q50 = 0,
  POC_q50 = 0,
  MOC_q50 = 0,
  Fe.OC_q50 = 0,
  Ca.OC_q50 = 0,
  Lignin.phenol_q50= 0,
  Microbial.necromass_q50 = 0
)
#Others
min <- data.frame( 
  Category = "Min",
  Lon_q50 = 0,
  Lat_q50 = 0,
  MAT_q50 = 0,
  MAP_q50 = 0,
  人口密度_q50 = 0,
  总种植面积_q50=0
)

#Geop.Chem
max <- data.frame(
  Category = "max",
  EC_q50 = 100,
  pH_q50 = 100,
  TN_q50 = 100,
  P_q50 = 100,
  C.P_q50 = 100,
  C.N_q50 = 100,
  Feo_q50 = 100,
  Fep_q50= 100,
  Fed_q50 = 100,
  Ald_q50 = 100,
  Alo_q50= 100,
  Alp_q50= 100
)
#Org.C
max <- data.frame( 
  Category = "Max",
  TOC_q50 = 100,
  DOC_q50 = 100,
  POC_q50 = 100,
  MOC_q50 = 100,
  Fe.OC_q50 = 100,
  Ca.OC_q50 = 100,
  Lignin.phenol_q50= 100,
  Microbial.necromass_q50 = 100
)
#Others
max <- data.frame( 
  Category = "Max",
  Lon_q50 = 100,
  Lat_q50 = 100,
  MAT_q50 = 100,
  MAP_q50 = 100,
  人口密度_q50 = 100,
  总种植面积_q50= 100
)

## Replace spr into SRP
colnames(median_only_rescaled)[colnames(median_only_rescaled) == "SRP_q50"] <- "srp_q50"

## Now bind the different columns together
bind_columns <- rbind(max,min,median_only_rescaled)
rownames(bind_columns) <- bind_columns$Category
bind_columns$Category <- NULL

## 定义颜色和标题
colors <- c("#E3F3CD","#81BAD8", "#FDAE61")
titles <- c("MF","SG","SM")

## 调整边距
op <- par(mar = c(1, 1, 1, 1))
par(mfrow = c(1, 1))  # 创建1行2列的布局

## 创建雷达图的函数
create_beautiful_radarchart <- function(data, color = "#00AFBB", vlabels = colnames(data), vlcex = 0.7,
   caxislabels = NULL, title = NULL, ...){
   radarchart(
    data, axistype = 1,
    # 自定义多边形
    pcol = color, pfcol = scales::alpha(color, 0.5), plwd = 2, plty = 1,
    # 自定义网格
    cglcol = "grey", cglty = 1, cglwd = 0.8,
    # 自定义坐标轴
    axislabcol = "grey", 
    # 变量标签
    vlcex = vlcex, vlabels = vlabels,
    caxislabels = caxislabels, title = title, ...
  )
}  

## 绘制雷达图
for(i in 1:3){
  create_beautiful_radarchart(
    data = bind_columns[c(1, 2, i+2), ], 
    caxislabels = c(0, 25, 50, 75, 100),
    color = colors[i], 
    title = titles[i]
  )
}

print(bind_columns[c(1, 2, 3), ])  # max, min, MF
print(bind_columns[c(1, 2, 4), ])  # max, min, SG  
print(bind_columns[c(1, 2, 5), ])  # max, min, SM
  