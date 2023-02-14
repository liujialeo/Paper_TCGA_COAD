rm(list = ls())
load(file = "exprSet_diff.Rda")
test <- exprSet_diff[1:10,1:10]
exprSet <- exprSet_diff

batch_cor <- function(gene){
  y <- as.numeric(exprSet[gene,])
  rownames <- rownames(exprSet)
  do.call(rbind,future_lapply(rownames, function(x){
    dd  <- cor.test(as.numeric(exprSet[x,]),y,type="spearman")
    data.frame(gene=gene,mRNAs=x,cor=dd$estimate,p.value=dd$p.value )
  }))
}
library(future.apply)
plan(multiprocess)
system.time(dd <- batch_cor("SQLE"))


###将与SQLE的相关性基因进行GSEA分析
gene <- dd$mRNAs
## 转换
library(clusterProfiler)
gene = bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
## 去重
gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)

gene_df <- data.frame(logFC=dd$cor,
                      SYMBOL = dd$mRNAs)
gene_df <- merge(gene_df,gene,by="SYMBOL")

## geneList 三部曲
## 1.获取基因logFC
geneList <- gene_df$logFC
## 2.命名
names(geneList) = gene_df$ENTREZID
## 3.排序很重要
geneList = sort(geneList, decreasing = TRUE)

library(clusterProfiler)
## 读入hallmarks gene set，从哪来？
hallmarks <- read.gmt("h.all.v7.1.entrez.gmt")
# 需要网络
y <- GSEA(geneList,TERM2GENE =hallmarks)
### 看整体分布
library(ggplot2)
dotplot(y,showCategory=12,split=".sign")+facet_grid(~.sign)
##将图片输出为PDF
library(export)
graph2ppt(file="SQLE-GSEA.pptx", width=10.5, height=7)

yd <- data.frame(y)


library(enrichplot)
gseaplot2(y,"HALLMARK_PI3K_AKT_MTOR_SIGNALING",color = "red",pvalue_table = T, base_size=21)

gseaplot2(y,"HALLMARK_MTORC1_SIGNALING",color = "red",pvalue_table = T, base_size=21)

