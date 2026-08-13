# 加载必要的包
library(ggplot2)
library(dplyr)
library(tidyr)
library(openxlsx)

#设置工作目录
setwd("C:/Users/dang/Documents")

#数据很全
#Pmax Q10
data<-read.table("Q10(116).txt",header=T,sep="\t",quote="")
#木质素值已更新
data<-mutate(data,Lignin.phenol=Lignin.phenol/1000)#列数值处理


#添加DOM
DOM<-read.table("para-1026.txt",header=T,sep="\t",quote="")
DOM6<-DOM[,c(1,4,5,115,121,124,127)]
data <- merge(data,DOM6,by="ID",all=T)
data<-data[-c(68:70),]


#Gompertz
Gompertz<-read.table("Gompertz.txt",header=T,sep="\t",quote="")
Gompertz<-Gompertz %>% rename( ID = Sample)
gompertz_subset <- Gompertz[, c(1,3:5)]

data <- merge(gompertz_subset, data, by = "ID")#匹配样本数量的对话框
#Stack1
Stack1<-read.table("Stack1.txt",header=T,sep="\t",quote="")
data <- merge(Stack1, data, by = "ID")
Stack2<-read.table("Stack2.txt",header=T,sep="\t",quote="")
data <- merge(Stack2, data, by = "ID")
row.names(data)<-data$ID


#write.xlsx(data, "Stack3.xlsx")

#Stack
data<-read.table("Stack.txt",header=T,sep="\t",quote="")
#Humanity <- read.table("Humanity.txt", header=TRUE, sep="\t", quote="", fileEncoding="GBK")
Humanity <- read.table("Humanity2025.7.13.txt", header=TRUE, sep="\t", quote="", fileEncoding="GBK")
#千公顷
Humanity<-mutate(Humanity,总种植面积=总种植面积/1000)
data <- merge(Humanity, data, by = "ID")
row.names(data)<-data$ID
data<-data[,-1]


write.xlsx(data, "Blue carbon116.xlsx")
write.xlsx(data, "Blue carbon116.all.xlsx")
#write.xlsx(data, "βNTI.xlsx")

# 转换为长格式
data_long <- data %>%
  pivot_longer(cols = -Category, names_to = "Variable", values_to = "Value")

data_long <- data_long %>%
  filter(!if_any(everything(), is.na))

library(RColorBrewer)  # 如果需要使用ColorBrewer调色板

# 定义一个颜色映射，将Category映射到颜色
category_colors <- c(
  "SM" = "#F58383",
  "MF" = "#FAEE85",
  "SG" = "#A8CAE8"
)

category_colors <- c(
  "SM" = "#c82423",
  "MF" = "#ffbe7a",
  "SG" = "#3480b8"
)

category_colors <- c(
  "SM" = "#fa8878",
  "MF" = "#ffbe7a",
  "SG" = "#82afda"
)

#老师选择
category_colors <- c(
  "MF" = "#E3F3CD",
  "SM" = "#FDAE61",
  "SG" = "#81BAD8"
)

#Lignin.phenol
filtered_data <- data_long %>% filter(Variable == "Lignin.phenol")
filtered_data$Category <- factor(filtered_data$Category, levels = c("MF", "SM", "SG"))
summary(filtered_data)
sum(is.na(filtered_data$Value))

#00000000000000000000000000
# 计算每个类别的中位数
medians <- filtered_data %>%
  group_by(Category) %>%
  summarise(MedianValue = median(Value, na.rm = TRUE), .groups = 'drop')

# 绘制堆叠直方图，并去除网格线
ggplot(filtered_data, aes(x = Value, fill = Category)) +
  geom_histogram(binwidth = 0.18,boundary = 0, position = "stack", alpha = 1) +
  scale_fill_manual(values = category_colors) +
  coord_cartesian(xlim = c(0, 4.3))+
  scale_x_continuous(breaks = seq(0, 4.3, by = 0.5))+
  theme_minimal() +
  labs(x = "Lignin.phenol", y = "No. Systerms") +
  theme(
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 为每个类别添加中位线，颜色随Category变化
  lapply(unique(filtered_data$Category), function(cat) {
    med_val <- medians$MedianValue[medians$Category == cat]
    color <- category_colors[cat] # 获取对应类别的颜色
    geom_segment(aes(x = med_val, xend = med_val, y = -Inf, yend = Inf), 
                 linetype = "dashed", colour = color,linewidth = 0.5) # 使用对应的颜色
  })


#Cinnamyl
filtered_data <- data_long %>% filter(Variable == "Cinnamyl")
filtered_data$Category <- factor(filtered_data$Category, levels = c("MF", "SM", "SG"))
summary(filtered_data)
sum(is.na(filtered_data$Value))

#00000000000000000000000000
# 计算每个类别的中位数
medians <- filtered_data %>%
  group_by(Category) %>%
  summarise(MedianValue = median(Value, na.rm = TRUE), .groups = 'drop')

# 绘制堆叠直方图，并去除网格线
ggplot(filtered_data, aes(x = Value, fill = Category)) +
  geom_histogram(binwidth = 10,boundary = 0, position = "stack", alpha = 1) +
  scale_fill_manual(values = category_colors) +
  coord_cartesian(xlim = c(0, 327))+
  scale_x_continuous(breaks = seq(0, 326, by = 30))+
  theme_minimal() +
  labs(x = "Cinnamyl", y = "No. Systerms") +
  theme(
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 为每个类别添加中位线，颜色随Category变化
  lapply(unique(filtered_data$Category), function(cat) {
    med_val <- medians$MedianValue[medians$Category == cat]
    color <- category_colors[cat] # 获取对应类别的颜色
    geom_segment(aes(x = med_val, xend = med_val, y = -Inf, yend = Inf), 
                 linetype = "dashed", colour = color,linewidth = 0.5) # 使用对应的颜色
  })


#Syringyl
filtered_data <- data_long %>% filter(Variable == "Syringyl")
filtered_data$Category <- factor(filtered_data$Category, levels = c("MF", "SM", "SG"))
summary(filtered_data)
sum(is.na(filtered_data$Value))

#00000000000000000000000000
# 计算每个类别的中位数
medians <- filtered_data %>%
  group_by(Category) %>%
  summarise(MedianValue = median(Value, na.rm = TRUE), .groups = 'drop')

# 绘制堆叠直方图，并去除网格线
ggplot(filtered_data, aes(x = Value, fill = Category)) +
  geom_histogram(binwidth = 35,boundary = 0, position = "stack", alpha = 1) +
  scale_fill_manual(values = category_colors) +
  coord_cartesian(xlim = c(0, 1050))+
  scale_x_continuous(breaks = seq(0, 1050, by = 100))+
  theme_minimal() +
  labs(x = "Syringyl", y = "No. Systerms") +
  theme(
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 为每个类别添加中位线，颜色随Category变化
  lapply(unique(filtered_data$Category), function(cat) {
    med_val <- medians$MedianValue[medians$Category == cat]
    color <- category_colors[cat] # 获取对应类别的颜色
    geom_segment(aes(x = med_val, xend = med_val, y = -Inf, yend = Inf), 
                 linetype = "dashed", colour = color,linewidth = 0.5) # 使用对应的颜色
  })


#Syringyl
filtered_data <- data_long %>% filter(Variable == "Syringyl")
filtered_data$Category <- factor(filtered_data$Category, levels = c("MF", "SM", "SG"))
summary(filtered_data)
sum(is.na(filtered_data$Value))

#00000000000000000000000000
# 计算每个类别的中位数
medians <- filtered_data %>%
  group_by(Category) %>%
  summarise(MedianValue = median(Value, na.rm = TRUE), .groups = 'drop')

# 绘制堆叠直方图，并去除网格线
ggplot(filtered_data, aes(x = Value, fill = Category)) +
  geom_histogram(binwidth = 35,boundary = 0, position = "stack", alpha = 1) +
  scale_fill_manual(values = category_colors) +
  coord_cartesian(xlim = c(0, 1050))+
  scale_x_continuous(breaks = seq(0, 1050, by = 100))+
  theme_minimal() +
  labs(x = "Syringyl", y = "No. Systerms") +
  theme(
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 为每个类别添加中位线，颜色随Category变化
  lapply(unique(filtered_data$Category), function(cat) {
    med_val <- medians$MedianValue[medians$Category == cat]
    color <- category_colors[cat] # 获取对应类别的颜色
    geom_segment(aes(x = med_val, xend = med_val, y = -Inf, yend = Inf), 
                 linetype = "dashed", colour = color,linewidth = 0.5) # 使用对应的颜色
  })


#Vanillyl
filtered_data <- data_long %>% filter(Variable == "Vanillyl")
filtered_data$Category <- factor(filtered_data$Category, levels = c("MF", "SM", "SG"))
summary(filtered_data)
sum(is.na(filtered_data$Value))

#00000000000000000000000000
# 计算每个类别的中位数
medians <- filtered_data %>%
  group_by(Category) %>%
  summarise(MedianValue = median(Value, na.rm = TRUE), .groups = 'drop')

# 绘制堆叠直方图，并去除网格线
ggplot(filtered_data, aes(x = Value, fill = Category)) +
  geom_histogram(binwidth = 35,boundary = 0, position = "stack", alpha = 1) +
  scale_fill_manual(values = category_colors) +
  coord_cartesian(xlim = c(0, 1120))+
  scale_x_continuous(breaks = seq(0, 1120, by = 100))+
  theme_minimal() +
  labs(x = "Vanillyl", y = "No. Systerms") +
  theme(
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 为每个类别添加中位线，颜色随Category变化
  lapply(unique(filtered_data$Category), function(cat) {
    med_val <- medians$MedianValue[medians$Category == cat]
    color <- category_colors[cat] # 获取对应类别的颜色
    geom_segment(aes(x = med_val, xend = med_val, y = -Inf, yend = Inf), 
                 linetype = "dashed", colour = color,linewidth = 0.5) # 使用对应的颜色
  })


#Microbial.necromass
filtered_data <- data_long %>% filter(Variable == "Microbial.necromass")
filtered_data$Category <- factor(filtered_data$Category, levels = c("MF", "SM", "SG"))
summary(filtered_data)
sum(is.na(filtered_data$Value))

#00000000000000000000000000
# 计算每个类别的中位数
medians <- filtered_data %>%
  group_by(Category) %>%
  summarise(MedianValue = median(Value, na.rm = TRUE), .groups = 'drop')

# 绘制堆叠直方图，并去除网格线
ggplot(filtered_data, aes(x = Value, fill = Category)) +
  geom_histogram(binwidth = 0.18,boundary = 0, position = "stack", alpha = 1) +
  scale_fill_manual(values = category_colors) +
  coord_cartesian(xlim = c(0, 5.5))+
  scale_x_continuous(breaks = seq(0, 5.5, by = 1))+
  theme_minimal() +
  labs(x = "Microbial.necromass", y = "No. Systerms") +
  theme(
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 为每个类别添加中位线，颜色随Category变化
  lapply(unique(filtered_data$Category), function(cat) {
    med_val <- medians$MedianValue[medians$Category == cat]
    color <- category_colors[cat] # 获取对应类别的颜色
    geom_segment(aes(x = med_val, xend = med_val, y = -Inf, yend = Inf), 
                 linetype = "dashed", colour = color,linewidth = 0.5) # 使用对应的颜色
  })


#Bacterial.necromass
filtered_data <- data_long %>% filter(Variable == "Bacterial.necromass")
filtered_data$Category <- factor(filtered_data$Category, levels = c("MF", "SM", "SG"))
summary(filtered_data)
sum(is.na(filtered_data$Value))

#00000000000000000000000000
# 计算每个类别的中位数
medians <- filtered_data %>%
  group_by(Category) %>%
  summarise(MedianValue = median(Value, na.rm = TRUE), .groups = 'drop')

# 绘制堆叠直方图，并去除网格线
ggplot(filtered_data, aes(x = Value, fill = Category)) +
  geom_histogram(binwidth = 0.04,boundary = 0, position = "stack", alpha = 1) +
  scale_fill_manual(values = category_colors) +
  coord_cartesian(xlim = c(0, 1.2))+
  scale_x_continuous(breaks = seq(0, 1.2, by = 0.1))+
  theme_minimal() +
  labs(x = "Bacterial.necromass", y = "No. Systerms") +
  theme(
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 为每个类别添加中位线，颜色随Category变化
  lapply(unique(filtered_data$Category), function(cat) {
    med_val <- medians$MedianValue[medians$Category == cat]
    color <- category_colors[cat] # 获取对应类别的颜色
    geom_segment(aes(x = med_val, xend = med_val, y = -Inf, yend = Inf), 
                 linetype = "dashed", colour = color,linewidth = 0.5) # 使用对应的颜色
  })


#Fungal.necromass
filtered_data <- data_long %>% filter(Variable == "Fungal.necromass")
filtered_data$Category <- factor(filtered_data$Category, levels = c("MF", "SM", "SG"))
summary(filtered_data)
sum(is.na(filtered_data$Value))

#00000000000000000000000000
# 计算每个类别的中位数
medians <- filtered_data %>%
  group_by(Category) %>%
  summarise(MedianValue = median(Value, na.rm = TRUE), .groups = 'drop')

# 绘制堆叠直方图，并去除网格线
ggplot(filtered_data, aes(x = Value, fill = Category)) +
  geom_histogram(binwidth = 0.15,boundary = 0, position = "stack", alpha = 1) +
  scale_fill_manual(values = category_colors) +
  coord_cartesian(xlim = c(0, 4.3))+
  scale_x_continuous(breaks = seq(0, 4.3, by = 0.8))+
  theme_minimal() +
  labs(x = "Fungal.necromass", y = "No. Systerms") +
  theme(
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 为每个类别添加中位线，颜色随Category变化
  lapply(unique(filtered_data$Category), function(cat) {
    med_val <- medians$MedianValue[medians$Category == cat]
    color <- category_colors[cat] # 获取对应类别的颜色
    geom_segment(aes(x = med_val, xend = med_val, y = -Inf, yend = Inf), 
                 linetype = "dashed", colour = color,linewidth = 0.5) # 使用对应的颜色
  })


#alpha.insim
filtered_data <- data_long %>% filter(Variable == "alpha.insim")
filtered_data$Category <- factor(filtered_data$Category, levels = c("MF", "SM", "SG"))
summary(filtered_data)
sum(is.na(filtered_data$Value))

#00000000000000000000000000
# 计算每个类别的中位数
medians <- filtered_data %>%
  group_by(Category) %>%
  summarise(MedianValue = median(Value, na.rm = TRUE), .groups = 'drop')

# 绘制堆叠直方图，并去除网格线
ggplot(filtered_data, aes(x = Value, fill = Category)) +
  geom_histogram(binwidth = 2,boundary = 0, position = "stack", alpha = 1) +
  scale_fill_manual(values = category_colors) +
  coord_cartesian(xlim = c(4, 63))+
  scale_x_continuous(breaks = seq(0, 63, by = 8))+
  theme_minimal() +
  labs(x = "alpha.insim", y = "No. Systerms") +
  theme(
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 为每个类别添加中位线，颜色随Category变化
  lapply(unique(filtered_data$Category), function(cat) {
    med_val <- medians$MedianValue[medians$Category == cat]
    color <- category_colors[cat] # 获取对应类别的颜色
    geom_segment(aes(x = med_val, xend = med_val, y = -Inf, yend = Inf), 
                 linetype = "dashed", colour = color,linewidth = 0.5) # 使用对应的颜色
  })


#alpha
filtered_data <- data_long %>% filter(Variable == "alpha")
filtered_data$Category <- factor(filtered_data$Category, levels = c("MF", "SM", "SG"))
summary(filtered_data)
sum(is.na(filtered_data$Value))

#00000000000000000000000000
# 计算每个类别的中位数
medians <- filtered_data %>%
  group_by(Category) %>%
  summarise(MedianValue = median(Value, na.rm = TRUE), .groups = 'drop')

# 绘制堆叠直方图，并去除网格线
ggplot(filtered_data, aes(x = Value, fill = Category)) +
  geom_histogram(binwidth = 0.1,boundary = 0, position = "stack", alpha = 1) +
  scale_fill_manual(values = category_colors) +
  coord_cartesian(xlim = c(2.4, 5))+
  scale_x_continuous(breaks = seq(2.4, 5, by = 0.5))+
  theme_minimal() +
  labs(x = "alpha", y = "No. Systerms") +
  theme(
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 为每个类别添加中位线，颜色随Category变化
  lapply(unique(filtered_data$Category), function(cat) {
    med_val <- medians$MedianValue[medians$Category == cat]
    color <- category_colors[cat] # 获取对应类别的颜色
    geom_segment(aes(x = med_val, xend = med_val, y = -Inf, yend = Inf), 
                 linetype = "dashed", colour = color,linewidth = 0.5) # 使用对应的颜色
  })


#Pcoa1.bin
filtered_data <- data_long %>% filter(Variable == "Pcoa1.bin")
filtered_data$Category <- factor(filtered_data$Category, levels = c("MF", "SM", "SG"))
summary(filtered_data)
sum(is.na(filtered_data$Value))

#00000000000000000000000000
# 计算每个类别的中位数
medians <- filtered_data %>%
  group_by(Category) %>%
  summarise(MedianValue = median(Value, na.rm = TRUE), .groups = 'drop')

# 绘制堆叠直方图，并去除网格线
ggplot(filtered_data, aes(x = Value, fill = Category)) +
  geom_histogram(binwidth = 0.03,boundary = 0, position = "stack", alpha = 1) +
  scale_fill_manual(values = category_colors) +
  coord_cartesian(xlim = c(-0.28, 0.48))+
  scale_x_continuous(breaks = seq(-0.28, 0.48, by = 0.2))+
  theme_minimal() +
  labs(x = "Pcoa1.bin", y = "No. Systerms") +
  theme(
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 为每个类别添加中位线，颜色随Category变化
  lapply(unique(filtered_data$Category), function(cat) {
    med_val <- medians$MedianValue[medians$Category == cat]
    color <- category_colors[cat] # 获取对应类别的颜色
    geom_segment(aes(x = med_val, xend = med_val, y = -Inf, yend = Inf), 
                 linetype = "dashed", colour = color,linewidth = 0.5) # 使用对应的颜色
  })


#Pcoa1.bin
filtered_data <- data_long %>% filter(Variable == "Pcoa1.bin")
filtered_data$Category <- factor(filtered_data$Category, levels = c("MF", "SM", "SG"))
summary(filtered_data)
sum(is.na(filtered_data$Value))

#00000000000000000000000000
# 计算每个类别的中位数
medians <- filtered_data %>%
  group_by(Category) %>%
  summarise(MedianValue = median(Value, na.rm = TRUE), .groups = 'drop')

# 绘制堆叠直方图，并去除网格线
ggplot(filtered_data, aes(x = Value, fill = Category)) +
  geom_histogram(binwidth = 0.03,boundary = 0, position = "stack", alpha = 1) +
  scale_fill_manual(values = category_colors) +
  coord_cartesian(xlim = c(-0.28, 0.48))+
  scale_x_continuous(breaks = seq(-0.28, 0.48, by = 0.2))+
  theme_minimal() +
  labs(x = "Pcoa1.bin", y = "No. Systerms") +
  theme(
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 为每个类别添加中位线，颜色随Category变化
  lapply(unique(filtered_data$Category), function(cat) {
    med_val <- medians$MedianValue[medians$Category == cat]
    color <- category_colors[cat] # 获取对应类别的颜色
    geom_segment(aes(x = med_val, xend = med_val, y = -Inf, yend = Inf), 
                 linetype = "dashed", colour = color,linewidth = 0.5) # 使用对应的颜色
  })


#Pcoa2.bin
filtered_data <- data_long %>% filter(Variable == "Pcoa2.bin")
filtered_data$Category <- factor(filtered_data$Category, levels = c("MF", "SM", "SG"))
summary(filtered_data)
sum(is.na(filtered_data$Value))

#00000000000000000000000000
# 计算每个类别的中位数
medians <- filtered_data %>%
  group_by(Category) %>%
  summarise(MedianValue = median(Value, na.rm = TRUE), .groups = 'drop')

# 绘制堆叠直方图，并去除网格线
ggplot(filtered_data, aes(x = Value, fill = Category)) +
  geom_histogram(binwidth = 0.025,boundary = 0, position = "stack", alpha = 1) +
  scale_fill_manual(values = category_colors) +
  coord_cartesian(xlim = c(-0.42, 0.31))+
  scale_x_continuous(breaks = seq(-0.42, 0.31, by = 0.1))+
  theme_minimal() +
  labs(x = "Pcoa2.bin", y = "No. Systerms") +
  theme(
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 为每个类别添加中位线，颜色随Category变化
  lapply(unique(filtered_data$Category), function(cat) {
    med_val <- medians$MedianValue[medians$Category == cat]
    color <- category_colors[cat] # 获取对应类别的颜色
    geom_segment(aes(x = med_val, xend = med_val, y = -Inf, yend = Inf), 
                 linetype = "dashed", colour = color,linewidth = 0.5) # 使用对应的颜色
  })


#alpha.gene.insim
filtered_data <- data_long %>% filter(Variable == "alpha.gene.insim")
filtered_data$Category <- factor(filtered_data$Category, levels = c("MF", "SM", "SG"))
summary(filtered_data)
sum(is.na(filtered_data$Value))

#00000000000000000000000000
# 计算每个类别的中位数
medians <- filtered_data %>%
  group_by(Category) %>%
  summarise(MedianValue = median(Value, na.rm = TRUE), .groups = 'drop')

# 绘制堆叠直方图，并去除网格线
ggplot(filtered_data, aes(x = Value, fill = Category)) +
  geom_histogram(binwidth = 20,boundary = 0, position = "stack", alpha = 1) +
  scale_fill_manual(values = category_colors) +
  coord_cartesian(xlim = c(1440, 1990))+
  scale_x_continuous(breaks = seq(1440, 1990, by = 100))+
  theme_minimal() +
  labs(x = "alpha.gene.insim", y = "No. Systerms") +
  theme(
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 为每个类别添加中位线，颜色随Category变化
  lapply(unique(filtered_data$Category), function(cat) {
    med_val <- medians$MedianValue[medians$Category == cat]
    color <- category_colors[cat] # 获取对应类别的颜色
    geom_segment(aes(x = med_val, xend = med_val, y = -Inf, yend = Inf), 
                 linetype = "dashed", colour = color,linewidth = 0.5) # 使用对应的颜色
  })


#alpha.gene
filtered_data <- data_long %>% filter(Variable == "alpha.gene")
filtered_data$Category <- factor(filtered_data$Category, levels = c("MF", "SM", "SG"))
summary(filtered_data)
sum(is.na(filtered_data$Value))

#00000000000000000000000000
# 计算每个类别的中位数
medians <- filtered_data %>%
  group_by(Category) %>%
  summarise(MedianValue = median(Value, na.rm = TRUE), .groups = 'drop')

# 绘制堆叠直方图，并去除网格线
ggplot(filtered_data, aes(x = Value, fill = Category)) +
  geom_histogram(binwidth = 0.012,boundary = 0, position = "stack", alpha = 1) +
  scale_fill_manual(values = category_colors) +
  coord_cartesian(xlim = c(7.73, 8.03))+
  scale_x_continuous(breaks = seq(7.73, 8.03, by = 0.1))+
  theme_minimal() +
  labs(x = "alpha.gene", y = "No. Systerms") +
  theme(
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 为每个类别添加中位线，颜色随Category变化
  lapply(unique(filtered_data$Category), function(cat) {
    med_val <- medians$MedianValue[medians$Category == cat]
    color <- category_colors[cat] # 获取对应类别的颜色
    geom_segment(aes(x = med_val, xend = med_val, y = -Inf, yend = Inf), 
                 linetype = "dashed", colour = color,linewidth = 0.5) # 使用对应的颜色
  })


#Pcoa1.gene
filtered_data <- data_long %>% filter(Variable == "Pcoa1.gene")
filtered_data$Category <- factor(filtered_data$Category, levels = c("MF", "SM", "SG"))
summary(filtered_data)
sum(is.na(filtered_data$Value))

#00000000000000000000000000
# 计算每个类别的中位数
medians <- filtered_data %>%
  group_by(Category) %>%
  summarise(MedianValue = median(Value, na.rm = TRUE), .groups = 'drop')

# 绘制堆叠直方图，并去除网格线
ggplot(filtered_data, aes(x = Value, fill = Category)) +
  geom_histogram(binwidth = 0.012,boundary = 0, position = "stack", alpha = 1) +
  scale_fill_manual(values = category_colors) +
  coord_cartesian(xlim = c(-0.16, 0.18))+
  scale_x_continuous(breaks = seq(-0.16, 0.18, by = 30))+
  theme_minimal() +
  labs(x = "Pcoa1.gene", y = "No. Systerms") +
  theme(
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 为每个类别添加中位线，颜色随Category变化
  lapply(unique(filtered_data$Category), function(cat) {
    med_val <- medians$MedianValue[medians$Category == cat]
    color <- category_colors[cat] # 获取对应类别的颜色
    geom_segment(aes(x = med_val, xend = med_val, y = -Inf, yend = Inf), 
                 linetype = "dashed", colour = color,linewidth = 0.5) # 使用对应的颜色
  })


#Pcoa2.gene
filtered_data <- data_long %>% filter(Variable == "Pcoa2.gene")
filtered_data$Category <- factor(filtered_data$Category, levels = c("MF", "SM", "SG"))
summary(filtered_data)
sum(is.na(filtered_data$Value))

#00000000000000000000000000
# 计算每个类别的中位数
medians <- filtered_data %>%
  group_by(Category) %>%
  summarise(MedianValue = median(Value, na.rm = TRUE), .groups = 'drop')

# 绘制堆叠直方图，并去除网格线
ggplot(filtered_data, aes(x = Value, fill = Category)) +
  geom_histogram(binwidth = 0.01,boundary = 0, position = "stack", alpha = 1) +
  scale_fill_manual(values = category_colors) +
  coord_cartesian(xlim = c(-0.22, 0.095))+
  scale_x_continuous(breaks = seq(-0.22, 0.095, by = 0.05))+
  theme_minimal() +
  labs(x = "Pcoa2.gene", y = "No. Systerms") +
  theme(
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 为每个类别添加中位线，颜色随Category变化
  lapply(unique(filtered_data$Category), function(cat) {
    med_val <- medians$MedianValue[medians$Category == cat]
    color <- category_colors[cat] # 获取对应类别的颜色
    geom_segment(aes(x = med_val, xend = med_val, y = -Inf, yend = Inf), 
                 linetype = "dashed", colour = color,linewidth = 0.5) # 使用对应的颜色
  })


#TOC
filtered_data <- data_long %>% filter(Variable == "TOC")
filtered_data$Category <- factor(filtered_data$Category, levels = c("MF", "SM", "SG"))
summary(filtered_data)
sum(is.na(filtered_data$Value))

#00000000000000000000000000
# 计算每个类别的中位数
medians <- filtered_data %>%
  group_by(Category) %>%
  summarise(MedianValue = median(Value, na.rm = TRUE), .groups = 'drop')

# 绘制堆叠直方图，并去除网格线
ggplot(filtered_data, aes(x = Value, fill = Category)) +
  geom_histogram(binwidth = 2,boundary = 0, position = "stack", alpha = 1) +
  scale_fill_manual(values = category_colors) +
  coord_cartesian(xlim = c(0, 55))+
  scale_x_continuous(breaks = seq(0, 55, by = 5))+
  theme_minimal() +
  labs(x = "TOC", y = "No. Systerms") +
  theme(
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 为每个类别添加中位线，颜色随Category变化
  lapply(unique(filtered_data$Category), function(cat) {
    med_val <- medians$MedianValue[medians$Category == cat]
    color <- category_colors[cat] # 获取对应类别的颜色
    geom_segment(aes(x = med_val, xend = med_val, y = -Inf, yend = Inf), 
                 linetype = "dashed", colour = color,linewidth = 0.5) # 使用对应的颜色
  })


#POC
filtered_data <- data_long %>% filter(Variable == "POC")
filtered_data$Category <- factor(filtered_data$Category, levels = c("MF", "SM", "SG"))
summary(filtered_data)
sum(is.na(filtered_data$Value))

#00000000000000000000000000
# 计算每个类别的中位数
medians <- filtered_data %>%
  group_by(Category) %>%
  summarise(MedianValue = median(Value, na.rm = TRUE), .groups = 'drop')

# 绘制堆叠直方图，并去除网格线
ggplot(filtered_data, aes(x = Value, fill = Category)) +
  geom_histogram(binwidth = 2,boundary = 0, position = "stack", alpha = 1) +
  scale_fill_manual(values = category_colors) +
  coord_cartesian(xlim = c(0, 31))+
  scale_x_continuous(breaks = seq(0, 31, by = 5))+
  theme_minimal() +
  labs(x = "POC", y = "No. Systerms") +
  theme(
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 为每个类别添加中位线，颜色随Category变化
  lapply(unique(filtered_data$Category), function(cat) {
    med_val <- medians$MedianValue[medians$Category == cat]
    color <- category_colors[cat] # 获取对应类别的颜色
    geom_segment(aes(x = med_val, xend = med_val, y = -Inf, yend = Inf), 
                 linetype = "dashed", colour = color,linewidth = 0.5) # 使用对应的颜色
  })


#MOC
filtered_data <- data_long %>% filter(Variable == "MOC")
filtered_data$Category <- factor(filtered_data$Category, levels = c("MF", "SM", "SG"))
summary(filtered_data)
sum(is.na(filtered_data$Value))

#00000000000000000000000000
# 计算每个类别的中位数
medians <- filtered_data %>%
  group_by(Category) %>%
  summarise(MedianValue = median(Value, na.rm = TRUE), .groups = 'drop')

# 绘制堆叠直方图，并去除网格线
ggplot(filtered_data, aes(x = Value, fill = Category)) +
  geom_histogram(binwidth = 1.2,boundary = 0, position = "stack", alpha = 1) +
  scale_fill_manual(values = category_colors) +
  coord_cartesian(xlim = c(0, 31.2))+
  scale_x_continuous(breaks = seq(0, 31.2, by = 5))+
  theme_minimal() +
  labs(x = "MOC", y = "No. Systerms") +
  theme(
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 为每个类别添加中位线，颜色随Category变化
  lapply(unique(filtered_data$Category), function(cat) {
    med_val <- medians$MedianValue[medians$Category == cat]
    color <- category_colors[cat] # 获取对应类别的颜色
    geom_segment(aes(x = med_val, xend = med_val, y = -Inf, yend = Inf), 
                 linetype = "dashed", colour = color,linewidth = 0.5) # 使用对应的颜色
  })


#MOC
filtered_data <- data_long %>% filter(Variable == "MOC")
filtered_data$Category <- factor(filtered_data$Category, levels = c("MF", "SM", "SG"))
summary(filtered_data)
sum(is.na(filtered_data$Value))

#00000000000000000000000000
# 计算每个类别的中位数
medians <- filtered_data %>%
  group_by(Category) %>%
  summarise(MedianValue = median(Value, na.rm = TRUE), .groups = 'drop')

# 绘制堆叠直方图，并去除网格线
ggplot(filtered_data, aes(x = Value, fill = Category)) +
  geom_histogram(binwidth = 1.2,boundary = 0, position = "stack", alpha = 1) +
  scale_fill_manual(values = category_colors) +
  coord_cartesian(xlim = c(0, 31.2))+
  scale_x_continuous(breaks = seq(0, 31.2, by = 5))+
  theme_minimal() +
  labs(x = "MOC", y = "No. Systerms") +
  theme(
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 为每个类别添加中位线，颜色随Category变化
  lapply(unique(filtered_data$Category), function(cat) {
    med_val <- medians$MedianValue[medians$Category == cat]
    color <- category_colors[cat] # 获取对应类别的颜色
    geom_segment(aes(x = med_val, xend = med_val, y = -Inf, yend = Inf), 
                 linetype = "dashed", colour = color,linewidth = 0.5) # 使用对应的颜色
  })


#TN
filtered_data <- data_long %>% filter(Variable == "TN")
filtered_data$Category <- factor(filtered_data$Category, levels = c("MF", "SM", "SG"))
summary(filtered_data)
sum(is.na(filtered_data$Value))

#00000000000000000000000000
# 计算每个类别的中位数
medians <- filtered_data %>%
  group_by(Category) %>%
  summarise(MedianValue = median(Value, na.rm = TRUE), .groups = 'drop')

# 绘制堆叠直方图，并去除网格线
ggplot(filtered_data, aes(x = Value, fill = Category)) +
  geom_histogram(binwidth = 0.12,boundary = 0, position = "stack", alpha = 1) +
  scale_fill_manual(values = category_colors) +
  coord_cartesian(xlim = c(0, 3.3))+
  scale_x_continuous(breaks = seq(0, 3.3, by = 0.5))+
  theme_minimal() +
  labs(x = "TN", y = "No. Systerms") +
  theme(
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 为每个类别添加中位线，颜色随Category变化
  lapply(unique(filtered_data$Category), function(cat) {
    med_val <- medians$MedianValue[medians$Category == cat]
    color <- category_colors[cat] # 获取对应类别的颜色
    geom_segment(aes(x = med_val, xend = med_val, y = -Inf, yend = Inf), 
                 linetype = "dashed", colour = color,linewidth = 0.5) # 使用对应的颜色
  })


#C.N
filtered_data <- data_long %>% filter(Variable == "C.N")
filtered_data$Category <- factor(filtered_data$Category, levels = c("MF", "SM", "SG"))
summary(filtered_data)
sum(is.na(filtered_data$Value))

#00000000000000000000000000
# 计算每个类别的中位数
medians <- filtered_data %>%
  group_by(Category) %>%
  summarise(MedianValue = median(Value, na.rm = TRUE), .groups = 'drop')

# 绘制堆叠直方图，并去除网格线
ggplot(filtered_data, aes(x = Value, fill = Category)) +
  geom_histogram(binwidth = 2.5,boundary = 0, position = "stack", alpha = 1) +
  scale_fill_manual(values = category_colors) +
  coord_cartesian(xlim = c(5, 71))+
  scale_x_continuous(breaks = seq(5, 71, by = 5))+
  theme_minimal() +
  labs(x = "C.N", y = "No. Systerms") +
  theme(
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 为每个类别添加中位线，颜色随Category变化
  lapply(unique(filtered_data$Category), function(cat) {
    med_val <- medians$MedianValue[medians$Category == cat]
    color <- category_colors[cat] # 获取对应类别的颜色
    geom_segment(aes(x = med_val, xend = med_val, y = -Inf, yend = Inf), 
                 linetype = "dashed", colour = color,linewidth = 0.5) # 使用对应的颜色
  })


#P
filtered_data <- data_long %>% filter(Variable == "P")
filtered_data$Category <- factor(filtered_data$Category, levels = c("MF", "SM", "SG"))
summary(filtered_data)
sum(is.na(filtered_data$Value))

#00000000000000000000000000
# 计算每个类别的中位数
medians <- filtered_data %>%
  group_by(Category) %>%
  summarise(MedianValue = median(Value, na.rm = TRUE), .groups = 'drop')

# 绘制堆叠直方图，并去除网格线
ggplot(filtered_data, aes(x = Value, fill = Category)) +
  geom_histogram(binwidth = 50,boundary = 0, position = "stack", alpha = 1) +
  scale_fill_manual(values = category_colors) +
  coord_cartesian(xlim = c(80, 1500))+
  scale_x_continuous(breaks = seq(80, 1500, by = 200))+
  theme_minimal() +
  labs(x = "P", y = "No. Systerms") +
  theme(
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 为每个类别添加中位线，颜色随Category变化
  lapply(unique(filtered_data$Category), function(cat) {
    med_val <- medians$MedianValue[medians$Category == cat]
    color <- category_colors[cat] # 获取对应类别的颜色
    geom_segment(aes(x = med_val, xend = med_val, y = -Inf, yend = Inf), 
                 linetype = "dashed", colour = color,linewidth = 0.5) # 使用对应的颜色
  })


#C.P
filtered_data <- data_long %>% filter(Variable == "C.P")
filtered_data$Category <- factor(filtered_data$Category, levels = c("MF", "SM", "SG"))
summary(filtered_data)
sum(is.na(filtered_data$Value))

#00000000000000000000000000
# 计算每个类别的中位数
medians <- filtered_data %>%
  group_by(Category) %>%
  summarise(MedianValue = median(Value, na.rm = TRUE), .groups = 'drop')

# 绘制堆叠直方图，并去除网格线
ggplot(filtered_data, aes(x = Value, fill = Category)) +
  geom_histogram(binwidth = 4.5,boundary = 0, position = "stack", alpha = 1) +
  scale_fill_manual(values = category_colors) +
  coord_cartesian(xlim = c(1, 124))+
  scale_x_continuous(breaks = seq(1, 124, by = 12))+
  theme_minimal() +
  labs(x = "C.P", y = "No. Systerms") +
  theme(
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 为每个类别添加中位线，颜色随Category变化
  lapply(unique(filtered_data$Category), function(cat) {
    med_val <- medians$MedianValue[medians$Category == cat]
    color <- category_colors[cat] # 获取对应类别的颜色
    geom_segment(aes(x = med_val, xend = med_val, y = -Inf, yend = Inf), 
                 linetype = "dashed", colour = color,linewidth = 0.5) # 使用对应的颜色
  })


#pH
filtered_data <- data_long %>% filter(Variable == "pH")
filtered_data$Category <- factor(filtered_data$Category, levels = c("MF", "SM", "SG"))
summary(filtered_data)
sum(is.na(filtered_data$Value))

#00000000000000000000000000
# 计算每个类别的中位数
medians <- filtered_data %>%
  group_by(Category) %>%
  summarise(MedianValue = median(Value, na.rm = TRUE), .groups = 'drop')

# 绘制堆叠直方图，并去除网格线
ggplot(filtered_data, aes(x = Value, fill = Category)) +
  geom_histogram(binwidth = 0.2,boundary = 0, position = "stack", alpha = 1) +
  scale_fill_manual(values = category_colors) +
  coord_cartesian(xlim = c(4, 9.2))+
  scale_x_continuous(breaks = seq(4, 9.2, by = 0.5))+
  theme_minimal() +
  labs(x = "pH", y = "No. Systerms") +
  theme(
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 为每个类别添加中位线，颜色随Category变化
  lapply(unique(filtered_data$Category), function(cat) {
    med_val <- medians$MedianValue[medians$Category == cat]
    color <- category_colors[cat] # 获取对应类别的颜色
    geom_segment(aes(x = med_val, xend = med_val, y = -Inf, yend = Inf), 
                 linetype = "dashed", colour = color,linewidth = 0.5) # 使用对应的颜色
  })


#EC
filtered_data <- data_long %>% filter(Variable == "EC")
filtered_data$Category <- factor(filtered_data$Category, levels = c("MF", "SM", "SG"))
summary(filtered_data)
sum(is.na(filtered_data$Value))

#00000000000000000000000000
# 计算每个类别的中位数
medians <- filtered_data %>%
  group_by(Category) %>%
  summarise(MedianValue = median(Value, na.rm = TRUE), .groups = 'drop')

# 绘制堆叠直方图，并去除网格线
ggplot(filtered_data, aes(x = Value, fill = Category)) +
  geom_histogram(binwidth = 0.5,boundary = 0, position = "stack", alpha = 1) +
  scale_fill_manual(values = category_colors) +
  coord_cartesian(xlim = c(0, 12.8))+
  scale_x_continuous(breaks = seq(0, 12.8, by = 2))+
  theme_minimal() +
  labs(x = "EC", y = "No. Systerms") +
  theme(
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 为每个类别添加中位线，颜色随Category变化
  lapply(unique(filtered_data$Category), function(cat) {
    med_val <- medians$MedianValue[medians$Category == cat]
    color <- category_colors[cat] # 获取对应类别的颜色
    geom_segment(aes(x = med_val, xend = med_val, y = -Inf, yend = Inf), 
                 linetype = "dashed", colour = color,linewidth = 0.5) # 使用对应的颜色
  })


#Fed
filtered_data <- data_long %>% filter(Variable == "Fed")
filtered_data$Category <- factor(filtered_data$Category, levels = c("MF", "SM", "SG"))
summary(filtered_data)
sum(is.na(filtered_data$Value))

#00000000000000000000000000
# 计算每个类别的中位数
medians <- filtered_data %>%
  group_by(Category) %>%
  summarise(MedianValue = median(Value, na.rm = TRUE), .groups = 'drop')

# 绘制堆叠直方图，并去除网格线
ggplot(filtered_data, aes(x = Value, fill = Category)) +
  geom_histogram(binwidth = 1.25, boundary = 0, position = "stack", alpha = 1) +
  scale_fill_manual(values = category_colors) +
  coord_cartesian(xlim = c(0, 34))+
  scale_x_continuous(breaks = seq(0, 34, by = 5))+
  theme_minimal() +
  labs(x = "Fed", y = "No. Systerms") +
  theme(
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 为每个类别添加中位线，颜色随Category变化
  lapply(unique(filtered_data$Category), function(cat) {
    med_val <- medians$MedianValue[medians$Category == cat]
    color <- category_colors[cat] # 获取对应类别的颜色
    geom_segment(aes(x = med_val, xend = med_val, y = -Inf, yend = Inf), 
                 linetype = "dashed", colour = color,linewidth = 0.5) # 使用对应的颜色
  })


#Feo
filtered_data <- data_long %>% filter(Variable == "Feo")
filtered_data$Category <- factor(filtered_data$Category, levels = c("MF", "SM", "SG"))
summary(filtered_data)
sum(is.na(filtered_data$Value))

#00000000000000000000000000
# 计算每个类别的中位数
medians <- filtered_data %>%
  group_by(Category) %>%
  summarise(MedianValue = median(Value, na.rm = TRUE), .groups = 'drop')

# 绘制堆叠直方图，并去除网格线
ggplot(filtered_data, aes(x = Value, fill = Category)) +
  geom_histogram(binwidth = 0.5,boundary = 0, position = "stack", alpha = 1) +
  scale_fill_manual(values = category_colors) +
  coord_cartesian(xlim = c(0, 15))+
  scale_x_continuous(breaks = seq(0, 15, by = 2))+
  theme_minimal() +
  labs(x = "Feo", y = "No. Systerms") +
  theme(
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 为每个类别添加中位线，颜色随Category变化
  lapply(unique(filtered_data$Category), function(cat) {
    med_val <- medians$MedianValue[medians$Category == cat]
    color <- category_colors[cat] # 获取对应类别的颜色
    geom_segment(aes(x = med_val, xend = med_val, y = -Inf, yend = Inf), 
                 linetype = "dashed", colour = color,linewidth = 0.5) # 使用对应的颜色
  })


#Fep
filtered_data <- data_long %>% filter(Variable == "Fep")
filtered_data$Category <- factor(filtered_data$Category, levels = c("MF", "SM", "SG"))
summary(filtered_data)
sum(is.na(filtered_data$Value))

#00000000000000000000000000
# 计算每个类别的中位数
medians <- filtered_data %>%
  group_by(Category) %>%
  summarise(MedianValue = median(Value, na.rm = TRUE), .groups = 'drop')

# 绘制堆叠直方图，并去除网格线
ggplot(filtered_data, aes(x = Value, fill = Category)) +
  geom_histogram(binwidth = 0.08,boundary = 0, position = "stack", alpha = 1) +
  scale_fill_manual(values = category_colors) +
  coord_cartesian(xlim = c(0, 2.22))+
  scale_x_continuous(breaks = seq(0, 2.22, by = 0.2))+
  theme_minimal() +
  labs(x = "Fep", y = "No. Systerms") +
  theme(
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 为每个类别添加中位线，颜色随Category变化
  lapply(unique(filtered_data$Category), function(cat) {
    med_val <- medians$MedianValue[medians$Category == cat]
    color <- category_colors[cat] # 获取对应类别的颜色
    geom_segment(aes(x = med_val, xend = med_val, y = -Inf, yend = Inf), 
                 linetype = "dashed", colour = color,linewidth = 0.5) # 使用对应的颜色
  })


#MeanSize
filtered_data <- data_long %>% filter(Variable == "MeanSize")
filtered_data$Category <- factor(filtered_data$Category, levels = c("MF", "SM", "SG"))
summary(filtered_data)
sum(is.na(filtered_data$Value))

#00000000000000000000000000
# 计算每个类别的中位数
medians <- filtered_data %>%
  group_by(Category) %>%
  summarise(MedianValue = median(Value, na.rm = TRUE), .groups = 'drop')

# 绘制堆叠直方图，并去除网格线
ggplot(filtered_data, aes(x = Value, fill = Category)) +
  geom_histogram(binwidth = 20,boundary = 0, position = "stack", alpha = 1) +
  scale_fill_manual(values = category_colors) +
  coord_cartesian(xlim = c(3, 604))+
  scale_x_continuous(breaks = seq(3, 604, by = 60))+
  theme_minimal() +
  labs(x = "MeanSize", y = "No. Systerms") +
  theme(
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 为每个类别添加中位线，颜色随Category变化
  lapply(unique(filtered_data$Category), function(cat) {
    med_val <- medians$MedianValue[medians$Category == cat]
    color <- category_colors[cat] # 获取对应类别的颜色
    geom_segment(aes(x = med_val, xend = med_val, y = -Inf, yend = Inf), 
                 linetype = "dashed", colour = color,linewidth = 0.5) # 使用对应的颜色
  })



#MeanSize
filtered_data <- data_long %>% filter(Variable == "MeanSize")
filtered_data$Category <- factor(filtered_data$Category, levels = c("MF", "SM", "SG"))
summary(filtered_data)
sum(is.na(filtered_data$Value))

#00000000000000000000000000
# 计算每个类别的中位数
medians <- filtered_data %>%
  group_by(Category) %>%
  summarise(MedianValue = median(Value, na.rm = TRUE), .groups = 'drop')

# 绘制堆叠直方图，并去除网格线
ggplot(filtered_data, aes(x = Value, fill = Category)) +
  geom_histogram(binwidth = 20,boundary = 0, position = "stack", alpha = 1) +
  scale_fill_manual(values = category_colors) +
  coord_cartesian(xlim = c(3, 604))+
  scale_x_continuous(breaks = seq(3, 604, by = 60))+
  theme_minimal() +
  labs(x = "MeanSize", y = "No. Systerms") +
  theme(
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 为每个类别添加中位线，颜色随Category变化
  lapply(unique(filtered_data$Category), function(cat) {
    med_val <- medians$MedianValue[medians$Category == cat]
    color <- category_colors[cat] # 获取对应类别的颜色
    geom_segment(aes(x = med_val, xend = med_val, y = -Inf, yend = Inf), 
                 linetype = "dashed", colour = color,linewidth = 0.5) # 使用对应的颜色
  })


#X.4
filtered_data <- data_long %>% filter(Variable == "X.4")
filtered_data$Category <- factor(filtered_data$Category, levels = c("MF", "SM", "SG"))
summary(filtered_data)
sum(is.na(filtered_data$Value))

#00000000000000000000000000
# 计算每个类别的中位数
medians <- filtered_data %>%
  group_by(Category) %>%
  summarise(MedianValue = median(Value, na.rm = TRUE), .groups = 'drop')

# 绘制堆叠直方图，并去除网格线
ggplot(filtered_data, aes(x = Value, fill = Category)) +
  geom_histogram(binwidth = 2,boundary = 0, position = "stack", alpha = 1) +
  scale_fill_manual(values = category_colors) +
  coord_cartesian(xlim = c(0, 50))+
  scale_x_continuous(breaks = seq(0, 50, by = 5))+
  theme_minimal() +
  labs(x = "X.4", y = "No. Systerms") +
  theme(
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 为每个类别添加中位线，颜色随Category变化
  lapply(unique(filtered_data$Category), function(cat) {
    med_val <- medians$MedianValue[medians$Category == cat]
    color <- category_colors[cat] # 获取对应类别的颜色
    geom_segment(aes(x = med_val, xend = med_val, y = -Inf, yend = Inf), 
                 linetype = "dashed", colour = color,linewidth = 0.5) # 使用对应的颜色
  })


#X.4
filtered_data <- data_long %>% filter(Variable == "X.4")
filtered_data$Category <- factor(filtered_data$Category, levels = c("MF", "SM", "SG"))
summary(filtered_data)
sum(is.na(filtered_data$Value))

#00000000000000000000000000
# 计算每个类别的中位数
medians <- filtered_data %>%
  group_by(Category) %>%
  summarise(MedianValue = median(Value, na.rm = TRUE), .groups = 'drop')

# 绘制堆叠直方图，并去除网格线
ggplot(filtered_data, aes(x = Value, fill = Category)) +
  geom_histogram(binwidth = 2,boundary = 0, position = "stack", alpha = 1) +
  scale_fill_manual(values = category_colors) +
  coord_cartesian(xlim = c(0, 50))+
  scale_x_continuous(breaks = seq(0, 50, by = 5))+
  theme_minimal() +
  labs(x = "X.4", y = "No. Systerms") +
  theme(
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 为每个类别添加中位线，颜色随Category变化
  lapply(unique(filtered_data$Category), function(cat) {
    med_val <- medians$MedianValue[medians$Category == cat]
    color <- category_colors[cat] # 获取对应类别的颜色
    geom_segment(aes(x = med_val, xend = med_val, y = -Inf, yend = Inf), 
                 linetype = "dashed", colour = color,linewidth = 0.5) # 使用对应的颜色
  })


#X4.63
filtered_data <- data_long %>% filter(Variable == "X4.63")
filtered_data$Category <- factor(filtered_data$Category, levels = c("MF", "SM", "SG"))
summary(filtered_data)
sum(is.na(filtered_data$Value))

#00000000000000000000000000
# 计算每个类别的中位数
medians <- filtered_data %>%
  group_by(Category) %>%
  summarise(MedianValue = median(Value, na.rm = TRUE), .groups = 'drop')

# 绘制堆叠直方图，并去除网格线
ggplot(filtered_data, aes(x = Value, fill = Category)) +
  geom_histogram(binwidth = 3.5,boundary = 0, position = "stack", alpha = 1) +
  scale_fill_manual(values = category_colors) +
  coord_cartesian(xlim = c(0, 82))+
  scale_x_continuous(breaks = seq(0, 82, by = 8))+
  theme_minimal() +
  labs(x = "X4.63", y = "No. Systerms") +
  theme(
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 为每个类别添加中位线，颜色随Category变化
  lapply(unique(filtered_data$Category), function(cat) {
    med_val <- medians$MedianValue[medians$Category == cat]
    color <- category_colors[cat] # 获取对应类别的颜色
    geom_segment(aes(x = med_val, xend = med_val, y = -Inf, yend = Inf), 
                 linetype = "dashed", colour = color,linewidth = 0.5) # 使用对应的颜色
  })


#X.63
filtered_data <- data_long %>% filter(Variable == "X.63")
filtered_data$Category <- factor(filtered_data$Category, levels = c("MF", "SM", "SG"))
summary(filtered_data)
sum(is.na(filtered_data$Value))

#00000000000000000000000000
# 计算每个类别的中位数
medians <- filtered_data %>%
  group_by(Category) %>%
  summarise(MedianValue = median(Value, na.rm = TRUE), .groups = 'drop')

# 绘制堆叠直方图，并去除网格线
ggplot(filtered_data, aes(x = Value, fill = Category)) +
  geom_histogram(binwidth = 4,boundary = 0, position = "stack", alpha = 1) +
  scale_fill_manual(values = category_colors) +
  coord_cartesian(xlim = c(0, 100))+
  scale_x_continuous(breaks = seq(0, 100, by = 10))+
  theme_minimal() +
  labs(x = "X.63", y = "No. Systerms") +
  theme(
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 为每个类别添加中位线，颜色随Category变化
  lapply(unique(filtered_data$Category), function(cat) {
    med_val <- medians$MedianValue[medians$Category == cat]
    color <- category_colors[cat] # 获取对应类别的颜色
    geom_segment(aes(x = med_val, xend = med_val, y = -Inf, yend = Inf), 
                 linetype = "dashed", colour = color,linewidth = 0.5) # 使用对应的颜色
  })


#Ideg
filtered_data <- data_long %>% filter(Variable == "Ideg")
filtered_data$Category <- factor(filtered_data$Category, levels = c("MF", "SM", "SG"))
summary(filtered_data)
sum(is.na(filtered_data$Value))

#00000000000000000000000000
# 计算每个类别的中位数
medians <- filtered_data %>%
  group_by(Category) %>%
  summarise(MedianValue = median(Value, na.rm = TRUE), .groups = 'drop')

# 绘制堆叠直方图，并去除网格线
ggplot(filtered_data, aes(x = Value, fill = Category)) +
  geom_histogram(binwidth = 0.03,boundary = 0, position = "stack", alpha = 1) +
  scale_fill_manual(values = category_colors) +
  coord_cartesian(xlim = c(0, 0.7))+
  scale_x_continuous(breaks = seq(0, 0.7, by = 0.1))+
  theme_minimal() +
  labs(x = "Ideg", y = "No. Systerms") +
  theme(
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 为每个类别添加中位线，颜色随Category变化
  lapply(unique(filtered_data$Category), function(cat) {
    med_val <- medians$MedianValue[medians$Category == cat]
    color <- category_colors[cat] # 获取对应类别的颜色
    geom_segment(aes(x = med_val, xend = med_val, y = -Inf, yend = Inf), 
                 linetype = "dashed", colour = color,linewidth = 0.5) # 使用对应的颜色
  })


#dG
filtered_data <- data_long %>% filter(Variable == "dG")
filtered_data$Category <- factor(filtered_data$Category, levels = c("MF", "SM", "SG"))
summary(filtered_data)
sum(is.na(filtered_data$Value))

#00000000000000000000000000
# 计算每个类别的中位数
medians <- filtered_data %>%
  group_by(Category) %>%
  summarise(MedianValue = median(Value, na.rm = TRUE), .groups = 'drop')

# 绘制堆叠直方图，并去除网格线
ggplot(filtered_data, aes(x = Value, fill = Category)) +
  geom_histogram(binwidth = 2.5,boundary = 0, position = "stack", alpha = 1) +
  scale_fill_manual(values = category_colors) +
  coord_cartesian(xlim = c(35, 92))+
  scale_x_continuous(breaks = seq(35, 92, by = 6))+
  theme_minimal() +
  labs(x = "dG", y = "No. Systerms") +
  theme(
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 为每个类别添加中位线，颜色随Category变化
  lapply(unique(filtered_data$Category), function(cat) {
    med_val <- medians$MedianValue[medians$Category == cat]
    color <- category_colors[cat] # 获取对应类别的颜色
    geom_segment(aes(x = med_val, xend = med_val, y = -Inf, yend = Inf), 
                 linetype = "dashed", colour = color,linewidth = 0.5) # 使用对应的颜色
  })


#Idex.terr
filtered_data <- data_long %>% filter(Variable == "Idex.terr")
filtered_data$Category <- factor(filtered_data$Category, levels = c("MF", "SM", "SG"))
summary(filtered_data)
sum(is.na(filtered_data$Value))

#00000000000000000000000000
# 计算每个类别的中位数
medians <- filtered_data %>%
  group_by(Category) %>%
  summarise(MedianValue = median(Value, na.rm = TRUE), .groups = 'drop')

# 绘制堆叠直方图，并去除网格线
ggplot(filtered_data, aes(x = Value, fill = Category)) +
  geom_histogram(binwidth = 0.035,boundary = 0, position = "stack", alpha = 1) +
  scale_fill_manual(values = category_colors) +
  coord_cartesian(xlim = c(0.07, 0.9))+
  scale_x_continuous(breaks = seq(0.07, 0.9, by = 0.1))+
  theme_minimal() +
  labs(x = "Idex.terr", y = "No. Systerms") +
  theme(
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 为每个类别添加中位线，颜色随Category变化
  lapply(unique(filtered_data$Category), function(cat) {
    med_val <- medians$MedianValue[medians$Category == cat]
    color <- category_colors[cat] # 获取对应类别的颜色
    geom_segment(aes(x = med_val, xend = med_val, y = -Inf, yend = Inf), 
                 linetype = "dashed", colour = color,linewidth = 0.5) # 使用对应的颜色
  })



#Df.AI.dom
filtered_data <- data_long %>% filter(Variable == "Df.AI.dom")
filtered_data$Category <- factor(filtered_data$Category, levels = c("MF", "SM", "SG"))
summary(filtered_data)
sum(is.na(filtered_data$Value))

#00000000000000000000000000
# 计算每个类别的中位数
medians <- filtered_data %>%
  group_by(Category) %>%
  summarise(MedianValue = median(Value, na.rm = TRUE), .groups = 'drop')

# 绘制堆叠直方图，并去除网格线
ggplot(filtered_data, aes(x = Value, fill = Category)) +
  geom_histogram(binwidth = 0.006,boundary = 0, position = "stack", alpha = 1) +
  scale_fill_manual(values = category_colors) +
  coord_cartesian(xlim = c(0.035, 0.168))+
  scale_x_continuous(breaks = seq(0.035, 0.168, by = 0.02))+
  theme_minimal() +
  labs(x = "Df.Al.dom", y = "No. Systerms") +
  theme(
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 为每个类别添加中位线，颜色随Category变化
  lapply(unique(filtered_data$Category), function(cat) {
    med_val <- medians$MedianValue[medians$Category == cat]
    color <- category_colors[cat] # 获取对应类别的颜色
    geom_segment(aes(x = med_val, xend = med_val, y = -Inf, yend = Inf), 
                 linetype = "dashed", colour = color,linewidth = 0.5) # 使用对应的颜色
  })


#Df.AI.dom
filtered_data <- data_long %>% filter(Variable == "Df.AI.dom")
filtered_data$Category <- factor(filtered_data$Category, levels = c("MF", "SM", "SG"))
summary(filtered_data)
sum(is.na(filtered_data$Value))

#00000000000000000000000000
# 计算每个类别的中位数
medians <- filtered_data %>%
  group_by(Category) %>%
  summarise(MedianValue = median(Value, na.rm = TRUE), .groups = 'drop')

# 绘制堆叠直方图，并去除网格线
ggplot(filtered_data, aes(x = Value, fill = Category)) +
  geom_histogram(binwidth = 0.0055,boundary = 0, position = "stack", alpha = 1) +
  scale_fill_manual(values = category_colors) +
  coord_cartesian(xlim = c(0.035, 0.168))+
  scale_x_continuous(breaks = seq(0.035, 0.168, by = 0.02))+
  theme_minimal() +
  labs(x = "Df.Al.dom", y = "No. Systerms") +
  theme(
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 为每个类别添加中位线，颜色随Category变化
  lapply(unique(filtered_data$Category), function(cat) {
    med_val <- medians$MedianValue[medians$Category == cat]
    color <- category_colors[cat] # 获取对应类别的颜色
    geom_segment(aes(x = med_val, xend = med_val, y = -Inf, yend = Inf), 
                 linetype = "dashed", colour = color,linewidth = 0.5) # 使用对应的颜色
  })



#Df.NOSC.dom
filtered_data <- data_long %>% filter(Variable == "Df.NOSC.dom")
filtered_data$Category <- factor(filtered_data$Category, levels = c("MF", "SM", "SG"))
summary(filtered_data)
sum(is.na(filtered_data$Value))

#00000000000000000000000000
# 计算每个类别的中位数
medians <- filtered_data %>%
  group_by(Category) %>%
  summarise(MedianValue = median(Value, na.rm = TRUE), .groups = 'drop')

# 绘制堆叠直方图，并去除网格线
ggplot(filtered_data, aes(x = Value, fill = Category)) +
  geom_histogram(binwidth = 0.013,boundary = 0, position = "stack", alpha = 1) +
  scale_fill_manual(values = category_colors) +
  coord_cartesian(xlim = c(0.07, 0.40))+
  scale_x_continuous(breaks = seq(0.07, 0.40, by = 0.05))+
  theme_minimal() +
  labs(x = "Df.NOSC.dom", y = "No. Systerms") +
  theme(
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 为每个类别添加中位线，颜色随Category变化
  lapply(unique(filtered_data$Category), function(cat) {
    med_val <- medians$MedianValue[medians$Category == cat]
    color <- category_colors[cat] # 获取对应类别的颜色
    geom_segment(aes(x = med_val, xend = med_val, y = -Inf, yend = Inf), 
                 linetype = "dashed", colour = color,linewidth = 0.5) # 使用对应的颜色
  })



#Df.HC.dom
filtered_data <- data_long %>% filter(Variable == "Df.HC.dom")
filtered_data$Category <- factor(filtered_data$Category, levels = c("MF", "SM", "SG"))
summary(filtered_data)
sum(is.na(filtered_data$Value))

#00000000000000000000000000
# 计算每个类别的中位数
medians <- filtered_data %>%
  group_by(Category) %>%
  summarise(MedianValue = median(Value, na.rm = TRUE), .groups = 'drop')

# 绘制堆叠直方图，并去除网格线
ggplot(filtered_data, aes(x = Value, fill = Category)) +
  geom_histogram(binwidth = 0.008,boundary = 0, position = "stack", alpha = 1) +
  scale_fill_manual(values = category_colors) +
  coord_cartesian(xlim = c(0.062, 0.26))+
  scale_x_continuous(breaks = seq(0.062, 0.26, by = 0.05))+
  theme_minimal() +
  labs(x = "Df.HC.dom", y = "No. Systerms") +
  theme(
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 为每个类别添加中位线，颜色随Category变化
  lapply(unique(filtered_data$Category), function(cat) {
    med_val <- medians$MedianValue[medians$Category == cat]
    color <- category_colors[cat] # 获取对应类别的颜色
    geom_segment(aes(x = med_val, xend = med_val, y = -Inf, yend = Inf), 
                 linetype = "dashed", colour = color,linewidth = 0.5) # 使用对应的颜色
  })


#补充 2025.01.17
#DOC
filtered_data <- data_long %>% filter(Variable == "DOC")
filtered_data$Category <- factor(filtered_data$Category, levels = c("MF", "SM", "SG"))
summary(filtered_data)
sum(is.na(filtered_data$Value))

#00000000000000000000000000
# 计算每个类别的中位数
medians <- filtered_data %>%
  group_by(Category) %>%
  summarise(MedianValue = median(Value, na.rm = TRUE), .groups = 'drop')

# 绘制堆叠直方图，并去除网格线
ggplot(filtered_data, aes(x = Value, fill = Category)) +
  geom_histogram(binwidth = 5,boundary = 0, position = "stack", alpha = 1) +
  scale_fill_manual(values = category_colors) +
  coord_cartesian(xlim = c(5, 131.5))+
  scale_x_continuous(breaks = seq(5, 131.5, by = 20))+
  theme_minimal() +
  labs(x = "DOC", y = "No. Systerms") +
  theme(
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 为每个类别添加中位线，颜色随Category变化
  lapply(unique(filtered_data$Category), function(cat) {
    med_val <- medians$MedianValue[medians$Category == cat]
    color <- category_colors[cat] # 获取对应类别的颜色
    geom_segment(aes(x = med_val, xend = med_val, y = -Inf, yend = Inf), 
                 linetype = "dashed", colour = color,linewidth = 0.5) # 使用对应的颜色
  })



#Fe.OC
filtered_data <- data_long %>% filter(Variable == "Fe.OC")
filtered_data$Category <- factor(filtered_data$Category, levels = c("MF", "SM", "SG"))
summary(filtered_data)
sum(is.na(filtered_data$Value))

#00000000000000000000000000
# 计算每个类别的中位数
medians <- filtered_data %>%
  group_by(Category) %>%
  summarise(MedianValue = median(Value, na.rm = TRUE), .groups = 'drop')

# 绘制堆叠直方图，并去除网格线
ggplot(filtered_data, aes(x = Value, fill = Category)) +
  geom_histogram(binwidth = 0.45,boundary = 0, position = "stack", alpha = 1) +
  scale_fill_manual(values = category_colors) +
  coord_cartesian(xlim = c(0, 11.3))+
  scale_x_continuous(breaks = seq(0, 11.3, by = 1))+
  theme_minimal() +
  labs(x = "Fe.OC", y = "No. Systerms") +
  theme(
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 为每个类别添加中位线，颜色随Category变化
  lapply(unique(filtered_data$Category), function(cat) {
    med_val <- medians$MedianValue[medians$Category == cat]
    color <- category_colors[cat] # 获取对应类别的颜色
    geom_segment(aes(x = med_val, xend = med_val, y = -Inf, yend = Inf), 
                 linetype = "dashed", colour = color,linewidth = 0.5) # 使用对应的颜色
  })



#Ca.OC
filtered_data <- data_long %>% filter(Variable == "Ca.OC")
filtered_data$Category <- factor(filtered_data$Category, levels = c("MF", "SM", "SG"))
summary(filtered_data)
sum(is.na(filtered_data$Value))

#00000000000000000000000000
# 计算每个类别的中位数
medians <- filtered_data %>%
  group_by(Category) %>%
  summarise(MedianValue = median(Value, na.rm = TRUE), .groups = 'drop')

# 绘制堆叠直方图，并去除网格线
ggplot(filtered_data, aes(x = Value, fill = Category)) +
  geom_histogram(binwidth = 0.05,boundary = 0, position = "stack", alpha = 1) +
  scale_fill_manual(values = category_colors) +
  coord_cartesian(xlim = c(0.01, 1.2))+
  scale_x_continuous(breaks = seq(0.01, 1.2, by = 0.2))+
  theme_minimal() +
  labs(x = "Ca.OC", y = "No. Systerms") +
  theme(
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 为每个类别添加中位线，颜色随Category变化
  lapply(unique(filtered_data$Category), function(cat) {
    med_val <- medians$MedianValue[medians$Category == cat]
    color <- category_colors[cat] # 获取对应类别的颜色
    geom_segment(aes(x = med_val, xend = med_val, y = -Inf, yend = Inf), 
                 linetype = "dashed", colour = color,linewidth = 0.5) # 使用对应的颜色
  })


#Ald
filtered_data <- data_long %>% filter(Variable == "Ald")
filtered_data$Category <- factor(filtered_data$Category, levels = c("MF", "SM", "SG"))
summary(filtered_data)
sum(is.na(filtered_data$Value))

#00000000000000000000000000
# 计算每个类别的中位数
medians <- filtered_data %>%
  group_by(Category) %>%
  summarise(MedianValue = median(Value, na.rm = TRUE), .groups = 'drop')

# 绘制堆叠直方图，并去除网格线
ggplot(filtered_data, aes(x = Value, fill = Category)) +
  geom_histogram(binwidth = 0.13,boundary = 0, position = "stack", alpha = 1) +
  scale_fill_manual(values = category_colors) +
  coord_cartesian(xlim = c(0, 3.2))+
  scale_x_continuous(breaks = seq(0, 3.2, by = 0.4))+
  theme_minimal() +
  labs(x = "Ald", y = "No. Systerms") +
  theme(
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 为每个类别添加中位线，颜色随Category变化
  lapply(unique(filtered_data$Category), function(cat) {
    med_val <- medians$MedianValue[medians$Category == cat]
    color <- category_colors[cat] # 获取对应类别的颜色
    geom_segment(aes(x = med_val, xend = med_val, y = -Inf, yend = Inf), 
                 linetype = "dashed", colour = color,linewidth = 0.5) # 使用对应的颜色
  })


#Alo
filtered_data <- data_long %>% filter(Variable == "Alo")
filtered_data$Category <- factor(filtered_data$Category, levels = c("MF", "SM", "SG"))
summary(filtered_data)
sum(is.na(filtered_data$Value))

#00000000000000000000000000
# 计算每个类别的中位数
medians <- filtered_data %>%
  group_by(Category) %>%
  summarise(MedianValue = median(Value, na.rm = TRUE), .groups = 'drop')

# 绘制堆叠直方图，并去除网格线
ggplot(filtered_data, aes(x = Value, fill = Category)) +
  geom_histogram(binwidth = 0.23,boundary = 0, position = "stack", alpha = 1) +
  scale_fill_manual(values = category_colors) +
  coord_cartesian(xlim = c(0, 5.6))+
  scale_x_continuous(breaks = seq(0, 5.6, by = 0.6))+
  theme_minimal() +
  labs(x = "Alo", y = "No. Systerms") +
  theme(
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 为每个类别添加中位线，颜色随Category变化
  lapply(unique(filtered_data$Category), function(cat) {
    med_val <- medians$MedianValue[medians$Category == cat]
    color <- category_colors[cat] # 获取对应类别的颜色
    geom_segment(aes(x = med_val, xend = med_val, y = -Inf, yend = Inf), 
                 linetype = "dashed", colour = color,linewidth = 0.5) # 使用对应的颜色
  })


#Alp
filtered_data <- data_long %>% filter(Variable == "Alp")
filtered_data$Category <- factor(filtered_data$Category, levels = c("MF", "SM", "SG"))
summary(filtered_data)
sum(is.na(filtered_data$Value))

#00000000000000000000000000
# 计算每个类别的中位数
medians <- filtered_data %>%
  group_by(Category) %>%
  summarise(MedianValue = median(Value, na.rm = TRUE), .groups = 'drop')

# 绘制堆叠直方图，并去除网格线
ggplot(filtered_data, aes(x = Value, fill = Category)) +
  geom_histogram(binwidth = 0.06,boundary = 0, position = "stack", alpha = 1) +
  scale_fill_manual(values = category_colors) +
  coord_cartesian(xlim = c(0, 1.4))+
  scale_x_continuous(breaks = seq(0, 1.4, by = 0.2))+
  theme_minimal() +
  labs(x = "Alp", y = "No. Systerms") +
  theme(
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 为每个类别添加中位线，颜色随Category变化
  lapply(unique(filtered_data$Category), function(cat) {
    med_val <- medians$MedianValue[medians$Category == cat]
    color <- category_colors[cat] # 获取对应类别的颜色
    geom_segment(aes(x = med_val, xend = med_val, y = -Inf, yend = Inf), 
                 linetype = "dashed", colour = color,linewidth = 0.5) # 使用对应的颜色
  })


#Longitude
filtered_data <- data_long %>% filter(Variable == "Lon")
filtered_data$Category <- factor(filtered_data$Category, levels = c("MF", "SM", "SG"))
summary(filtered_data)
sum(is.na(filtered_data$Value))

#00000000000000000000000000
# 计算每个类别的中位数
medians <- filtered_data %>%
  group_by(Category) %>%
  summarise(MedianValue = median(Value, na.rm = TRUE), .groups = 'drop')

# 绘制堆叠直方图，并去除网格线
ggplot(filtered_data, aes(x = Value, fill = Category)) +
  geom_histogram(binwidth = 0.7,boundary = 0, position = "stack", alpha = 1) +
  scale_fill_manual(values = category_colors) +
  coord_cartesian(xlim = c(108, 125))+
  scale_x_continuous(breaks = seq(108, 125, by = 2))+
  theme_minimal() +
  labs(x = "Lon", y = "No. Systerms") +
  theme(
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 为每个类别添加中位线，颜色随Category变化
  lapply(unique(filtered_data$Category), function(cat) {
    med_val <- medians$MedianValue[medians$Category == cat]
    color <- category_colors[cat] # 获取对应类别的颜色
    geom_segment(aes(x = med_val, xend = med_val, y = -Inf, yend = Inf), 
                 linetype = "dashed", colour = color,linewidth = 0.5) # 使用对应的颜色
  })


#Latitude
filtered_data <- data_long %>% filter(Variable == "Lat")
filtered_data$Category <- factor(filtered_data$Category, levels = c("MF", "SM", "SG"))
summary(filtered_data)
sum(is.na(filtered_data$Value))

#00000000000000000000000000
# 计算每个类别的中位数
medians <- filtered_data %>%
  group_by(Category) %>%
  summarise(MedianValue = median(Value, na.rm = TRUE), .groups = 'drop')

# 绘制堆叠直方图，并去除网格线
ggplot(filtered_data, aes(x = Value, fill = Category)) +
  geom_histogram(binwidth = 1,boundary = 0, position = "stack", alpha = 1) +
  scale_fill_manual(values = category_colors) +
  coord_cartesian(xlim = c(18, 41.2))+
  scale_x_continuous(breaks = seq(18, 41.2, by = 2))+
  theme_minimal() +
  labs(x = "Lat", y = "No. Systerms") +
  theme(
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 为每个类别添加中位线，颜色随Category变化
  lapply(unique(filtered_data$Category), function(cat) {
    med_val <- medians$MedianValue[medians$Category == cat]
    color <- category_colors[cat] # 获取对应类别的颜色
    geom_segment(aes(x = med_val, xend = med_val, y = -Inf, yend = Inf), 
                 linetype = "dashed", colour = color,linewidth = 0.5) # 使用对应的颜色
  })


#MAT
filtered_data <- data_long %>% filter(Variable == "MAT")
filtered_data$Category <- factor(filtered_data$Category, levels = c("MF", "SM", "SG"))
summary(filtered_data)
sum(is.na(filtered_data$Value))

#00000000000000000000000000
# 计算每个类别的中位数
medians <- filtered_data %>%
  group_by(Category) %>%
  summarise(MedianValue = median(Value, na.rm = TRUE), .groups = 'drop')

# 绘制堆叠直方图，并去除网格线
ggplot(filtered_data, aes(x = Value, fill = Category)) +
  geom_histogram(binwidth = 0.7,boundary = 0, position = "stack", alpha = 1) +
  scale_fill_manual(values = category_colors) +
  coord_cartesian(xlim = c(8.5, 25))+
  scale_x_continuous(breaks = seq(8.5, 25, by = 2))+
  theme_minimal() +
  labs(x = "MAT", y = "No. Systerms") +
  theme(
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 为每个类别添加中位线，颜色随Category变化
  lapply(unique(filtered_data$Category), function(cat) {
    med_val <- medians$MedianValue[medians$Category == cat]
    color <- category_colors[cat] # 获取对应类别的颜色
    geom_segment(aes(x = med_val, xend = med_val, y = -Inf, yend = Inf), 
                 linetype = "dashed", colour = color,linewidth = 0.5) # 使用对应的颜色
  })



#MAP
filtered_data <- data_long %>% filter(Variable == "MAP")
filtered_data$Category <- factor(filtered_data$Category, levels = c("MF", "SM", "SG"))
summary(filtered_data)
sum(is.na(filtered_data$Value))

#00000000000000000000000000
# 计算每个类别的中位数
medians <- filtered_data %>%
  group_by(Category) %>%
  summarise(MedianValue = median(Value, na.rm = TRUE), .groups = 'drop')

# 绘制堆叠直方图，并去除网格线
ggplot(filtered_data, aes(x = Value, fill = Category)) +
  geom_histogram(binwidth = 70,boundary = 0, position = "stack", alpha = 1) +
  scale_fill_manual(values = category_colors) +
  coord_cartesian(xlim = c(560, 2190))+
  scale_x_continuous(breaks = seq(560, 2190, by = 300))+
  theme_minimal() +
  labs(x = "MAP", y = "No. Systerms") +
  theme(
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 为每个类别添加中位线，颜色随Category变化
  lapply(unique(filtered_data$Category), function(cat) {
    med_val <- medians$MedianValue[medians$Category == cat]
    color <- category_colors[cat] # 获取对应类别的颜色
    geom_segment(aes(x = med_val, xend = med_val, y = -Inf, yend = Inf), 
                 linetype = "dashed", colour = color,linewidth = 0.5) # 使用对应的颜色
  })


#人口密度
filtered_data <- data_long %>% filter(Variable == "人口密度")
filtered_data$Category <- factor(filtered_data$Category, levels = c("MF", "SM", "SG"))
summary(filtered_data)
sum(is.na(filtered_data$Value))

#00000000000000000000000000
# 计算每个类别的中位数
medians <- filtered_data %>%
  group_by(Category) %>%
  summarise(MedianValue = median(Value, na.rm = TRUE), .groups = 'drop')

# 绘制堆叠直方图，并去除网格线
ggplot(filtered_data, aes(x = Value, fill = Category)) +
  geom_histogram(binwidth = 600,boundary = 0, position = "stack", alpha = 1) +
  scale_fill_manual(values = category_colors) +
  coord_cartesian(xlim = c(0, 17000))+
  scale_x_continuous(breaks = seq(0, 17000, by = 2000))+
  theme_minimal() +
  labs(x = "Human", y = "No. Systerms") +
  theme(
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 为每个类别添加中位线，颜色随Category变化
  lapply(unique(filtered_data$Category), function(cat) {
    med_val <- medians$MedianValue[medians$Category == cat]
    color <- category_colors[cat] # 获取对应类别的颜色
    geom_segment(aes(x = med_val, xend = med_val, y = -Inf, yend = Inf), 
                 linetype = "dashed", colour = color,linewidth = 0.5) # 使用对应的颜色
  })


#地区生产总值
filtered_data <- data_long %>% filter(Variable == "地区生产总值")
filtered_data$Category <- factor(filtered_data$Category, levels = c("MF", "SM", "SG"))
summary(filtered_data)
sum(is.na(filtered_data$Value))

#00000000000000000000000000
# 计算每个类别的中位数
medians <- filtered_data %>%
  group_by(Category) %>%
  summarise(MedianValue = median(Value, na.rm = TRUE), .groups = 'drop')

# 绘制堆叠直方图，并去除网格线
ggplot(filtered_data, aes(x = Value, fill = Category)) +
  geom_histogram(binwidth = 650,boundary = 0, position = "stack", alpha = 1) +
  scale_fill_manual(values = category_colors) +
  coord_cartesian(xlim = c(79, 16200))+
  scale_x_continuous(breaks = seq(79, 16200, by = 3000))+
  theme_minimal() +
  labs(x = "GDP", y = "No. Systerms") +
  theme(
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 为每个类别添加中位线，颜色随Category变化
  lapply(unique(filtered_data$Category), function(cat) {
    med_val <- medians$MedianValue[medians$Category == cat]
    color <- category_colors[cat] # 获取对应类别的颜色
    geom_segment(aes(x = med_val, xend = med_val, y = -Inf, yend = Inf), 
                 linetype = "dashed", colour = color,linewidth = 0.5) # 使用对应的颜色
  })



#设施农业种植占地面积
filtered_data <- data_long %>% filter(Variable == "设施农业种植占地面积")
filtered_data$Category <- factor(filtered_data$Category, levels = c("MF", "SM", "SG"))
summary(filtered_data)
sum(is.na(filtered_data$Value))

#00000000000000000000000000
# 计算每个类别的中位数
medians <- filtered_data %>%
  group_by(Category) %>%
  summarise(MedianValue = median(Value, na.rm = TRUE), .groups = 'drop')

# 绘制堆叠直方图，并去除网格线
ggplot(filtered_data, aes(x = Value, fill = Category)) +
  geom_histogram(binwidth = 10,boundary = 0, position = "stack", alpha = 1) +
  scale_fill_manual(values = category_colors) +
  coord_cartesian(xlim = c(0, 235))+
  scale_x_continuous(breaks = seq(0, 235, by = 20))+
  theme_minimal() +
  labs(x = "Area", y = "No. Systerms") +
  theme(
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 为每个类别添加中位线，颜色随Category变化
  lapply(unique(filtered_data$Category), function(cat) {
    med_val <- medians$MedianValue[medians$Category == cat]
    color <- category_colors[cat] # 获取对应类别的颜色
    geom_segment(aes(x = med_val, xend = med_val, y = -Inf, yend = Inf), 
                 linetype = "dashed", colour = color,linewidth = 0.5) # 使用对应的颜色
  })


#总种植面积
filtered_data <- data_long %>% filter(Variable == "总种植面积")
filtered_data$Category <- factor(filtered_data$Category, levels = c("MF", "SM", "SG"))
summary(filtered_data)
sum(is.na(filtered_data$Value))

#00000000000000000000000000
# 计算每个类别的中位数
medians <- filtered_data %>%
  group_by(Category) %>%
  summarise(MedianValue = median(Value, na.rm = TRUE), .groups = 'drop')

# 绘制堆叠直方图，并去除网格线
ggplot(filtered_data, aes(x = Value, fill = Category)) +
  geom_histogram(binwidth = 8000,boundary = 0, position = "stack", alpha = 1) +
  scale_fill_manual(values = category_colors) +
  coord_cartesian(xlim = c(0, 222000))+
  scale_x_continuous(breaks = seq(0, 222000, by = 50000))+
  theme_minimal() +
  labs(x = "Area", y = "No. Systerms") +
  theme(
    strip.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 为每个类别添加中位线，颜色随Category变化
  lapply(unique(filtered_data$Category), function(cat) {
    med_val <- medians$MedianValue[medians$Category == cat]
    color <- category_colors[cat] # 获取对应类别的颜色
    geom_segment(aes(x = med_val, xend = med_val, y = -Inf, yend = Inf), 
                 linetype = "dashed", colour = color,linewidth = 0.5) # 使用对应的颜色
  })

