rm(list = ls())

load(file = "exprSet_samples.Rda")
test <- exprSet_merge[1:10,1:10]

library(ggplot2)

exprSet <- exprSet_merge
#看BRCA1基因的是癌和癌旁的表达
#成功后，尝试一下 TP53，ERBB2，ESR1，修改y的值就可以
ggplot(exprSet,aes(x=sample,y=SQLE))+
  geom_boxplot()


#试试你知道的gene，要注意这里没有非编码的gene

## 美图是个刚需，任何时候都需要，但是不用刻意去学
## 试试下面这个更加简洁的代码，你会有新的认识。
# 需要安装这个包
if(! require("ggstatsplot")) install.packages("ggstatsplot")

## 开始！！
library(ggstatsplot)
ggbetweenstats(data = exprSet, 
               x = sample, 
               y = WNT11,
               title = "WNT11",
               mean.label.size = 20,
               outlier.coef = 1.5,
               outlier.shape = 19
               )
#library(export)
#graph2ppt(file="SQLE_expression.pptx", width=7, height=5)



###多个基因的表达量
rownames(exprSet) <- exprSet[,1]
exprSet <- exprSet[,-1]
test <- exprSet[1:10,1:10]
index <- which(colnames(exprSet) %in% c("NFKBIE","IL1A","IL6","TNF","TNFAIP3",
                                        "NFKB1","REL","CXCL10","TLR2","ICAM1"))
#index <- which(colnames(exprSet) %in% c("NFKBIA","NFKBIB","NFKBIE"))
expr8 <- exprSet[c(1,index)]

expr_gather <- tidyr::gather(expr8,key = Genenames, value = Geneexpr,-c(sample))

library(ggplot2)
# 箱线图,分基因查看癌和癌旁
ggplot(
  expr_gather, 
  aes(x = sample, y = Geneexpr,color=sample)
) +
  geom_boxplot()+
  facet_wrap(~Genenames,nrow = 2)
library(ggridges)
# 用山脊图，只考虑亚型
ggplot(
  expr_gather, 
  aes(x = Geneexpr, y = sample)
) +
  geom_density_ridges_gradient(
    aes(fill = ..x..), scale = 1, size = 0.3
  ) +
  scale_fill_gradientn(
    colours = c("#0D0887FF", "#CC4678FF", "#F0F921FF"),
    name = "Temp. [F]"
  )+
  facet_wrap(~Genenames,nrow = 2)

####################################调整配色R包！
if(!require(paletteer))install.packages("paletteer")
if(!require(scico))install.packages('scico')
if(!require(nord))install.packages('nord')
library(paletteer)
paletteer_c("scico::berlin", n = 10)
library(ggpubr)

library(ggsignif)
my_comparisons=list(c("Normal","Tumor"))
# 箱线图,分基因查看癌和癌旁

ggplot(
  expr_gather, 
  aes(x = sample, y = Geneexpr,color=sample)
) +
  geom_boxplot()+
  facet_wrap(~Genenames,nrow = 2) +
  #scale_color_paletteer_d("basetheme::minimal")
  scale_color_manual(values = c("deepskyblue", "orangered"))+
  #scale_color_manual(values = c("#63B8FF", "#FF6EB4"))+
  geom_jitter(position = position_jitter(0.1), alpha = 0.1, size=2)+
  #theme_classic()+
  #geom_signif(color="black",comparisons=my_comparisons,map_signif_level=T,test=t.test) 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif",method = "t.test")
#保存到pdf文件
ggsave("TNF_signaling_gene_expression.pdf", width = 10, height = 6)
ggsave("TNF_signaling_gene_expression.tiff", width = 10, height = 6)

#ggsave("IKBs_gene_expression.pdf", width = 5, height = 3)