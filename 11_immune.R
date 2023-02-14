### 处理表达量数据数据
rm(list = ls())
library(dplyr)
library(tidyr)
### 1.加载marker
load(file = "cellMarker_ssGSEA.Rdata")
### 2.加载表达量
load(file = "exprSet_diff.Rda")
exprSet <- exprSet_diff
test <- exprSet[1:10,1:10]

expr <- exprSet

expr <- as.matrix(expr)

library(GSVA)
gsva_data <- gsva(expr,cellMarker, method = "ssgsea")
##转置
tcga_gsva <- as.data.frame(t(gsva_data))

## 添加分组信息
clin <- as.data.frame(rownames(tcga_gsva)) 
clin[,2] <- c(rep("normal",41),rep("tumor",471)) 
colnames(clin) <- c("sample","group")
tcga_gsva <- cbind(rownames(tcga_gsva),tcga_gsva)
colnames(tcga_gsva)[1] <- "sample"

tcga_gsva <- inner_join(clin,tcga_gsva, by= "sample")
save(tcga_gsva,file = "tcga_gsva.Rdata")



### 特定癌种中，癌和癌旁的差异
rm(list = ls())
load(file = "tcga_gsva.Rdata")
test2 <- tcga_gsva[1:10,1:10]

dd <- tcga_gsva
### 调整数据
library(dplyr)
library(tidyr)
dd1 <- dd %>% 
  pivot_longer(cols=3:30,
               names_to= "celltype",
               values_to = "NES")
dd1 <- dd1[,-1]
colnames(dd1)[1] <- "sample"
library(ggplot2)
library(ggpubr)


### 箱线图
ggplot(data =dd1, aes(x = celltype, y = NES))+
  geom_boxplot(aes(fill = sample),outlier.shape = NA)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1, colour = "black"))+
  stat_compare_means(aes(group=sample), label = "p.signif")

### 小提琴
ggplot(data =dd1, aes(x = celltype, y = NES))+
  geom_violin(aes(fill = sample),position = position_dodge(1),scale = "width")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1, colour = "black"))+
  stat_compare_means(aes(group=sample), label = "p.signif")

### 混合叠加
ggplot(data =dd1, aes(x = celltype, y = NES))+
  geom_boxplot(aes(fill = sample),position = position_dodge(1),width=.3,outlier.shape = NA)+
  geom_violin(aes(colour = sample),position = position_dodge(1),scale = "width",fill=NA)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1, colour = "black"))+
  stat_compare_means(aes(group=sample), label = "p.signif")
