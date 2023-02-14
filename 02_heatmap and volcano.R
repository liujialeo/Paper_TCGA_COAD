#用行名提取数据
rm(list = ls())
load(file = "exprSet_diff.Rda")
load(file = "diffLab.Rda")
if(! require("pheatmap")) install.packages("pheatmap")
制作一个分组信息
heatdata <- exprSet_diff [rownames(diffLab),]
class <- c(rep("normal",41),rep("tumor",471)) 
annotation_col <- data.frame(class)
rownames(annotation_col) <- colnames(heatdata)

#如果注释出界，可以通过调整格子比例和字体修正
pheatmap(heatdata, #热图的数据
         cluster_rows = TRUE,#行聚类
         cluster_cols = F,#列聚类，可以看出样本之间的区分度
         annotation_col =annotation_col, #标注样本分类
         annotation_legend=F, # 显示注释
         show_rownames = F,# 显示行名
         show_colnames = F, 
         scale = "row", #以行来标准化，这个功能很不错
         color =colorRampPalette(c("green", "black","red"))(100),#调色
         #filename = "heatmap_F.pdf",#是否保存
         cellwidth = 0.4, cellheight = 0.2,# 格子比例
         fontsize = 10)

##volcano火山图
###############
###############
##用ggplot2
library(ggplot2)
library(ggrepel)
library(dplyr)
load(file = "allDiff.Rda")
data <- allDiff
data$significant <- as.factor(data$P.Value<0.05 & abs(data$logFC) > 1)
data$gene <- rownames(data)

ggplot(data=data, aes(x=logFC, y =-log10(P.Value),color=significant)) +
  geom_point(alpha=0.8, size=1.2)+
  scale_color_manual(values =c("black","red"))+
  labs(title="Volcanoplot", x="log2 (fold change)",y="-log10 (q-value)")+
  theme(plot.title = element_text(hjust = 0.4))+
  geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.6,alpha=0.8)+
  geom_vline(xintercept = c(1,-1),lty=4,lwd=0.6,alpha=0.8)+
  #theme(legend.position='none')
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),   
        axis.line = element_line(colour = "black"))+
  #geom_text(data=subset(data, abs(logFC) > 3), aes(label=gene),col="red",alpha = 1)+
  geom_text_repel(data=subset(data, abs(logFC) > 3), aes(label=gene),col="black",alpha = 0.8)
#ggsave("vocanol.pdf",,width = 7.09, height =5.6,dpi = 300)
