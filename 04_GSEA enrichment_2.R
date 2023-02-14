rm(list = ls())

##############################################################
load(file = "allDiff.Rda")
allDiff$gene <- rownames(allDiff)
library(dplyr)
gene_df <- allDiff %>%
  dplyr::select(logFC,gene) %>%
  ## 去掉NA
  filter(gene!="") %>%
  ## 去掉重复
  distinct(gene,.keep_all = T)

### 1.获取基因logFC
geneList <- gene_df$logFC
### 2.命名
names(geneList) = gene_df$gene
## 3.排序很重要
geneList = sort(geneList, decreasing = TRUE)

head(geneList)
library(clusterProfiler)


### 准备基因集
library(msigdbr)
msigdbr_species()
msigdbr_collections()
library(dplyr)
## Hallmarks gene set
h_df <- msigdbr(species = "Homo sapiens") %>% 
  filter(gs_cat == "H") %>% 
  dplyr::select(gs_name,gene_symbol)

y <- GSEA(geneList,TERM2GENE =h_df)
yd <- as.data.frame(y)
library(ggplot2)
dotplot(y,showCategory=15,split=".sign")+facet_grid(~.sign)
ggsave("output/GSEA_dotploot.pdf", width = 9, height = 10)

### 自定义画图
ggplot(y, showCategory = 30, aes(NES, forcats::fct_reorder(Description, NES))) + 
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_colour_gradientn(colours=c("#7fc97f","#984ea3")) +
  #scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
  #scale_color_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE))+
  scale_size_continuous(range=c(2, 10)) +
  theme_bw() + 
  xlab("Normalized Enrichment Score") +
  ylab(NULL)
ggsave("GSEA_dotploot_2.pdf", width = 9, height = 10)


library(enrichplot)
pathway.id = "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
gseaplot2(y,color = "green",geneSetID = pathway.id,pvalue_table = T)
ggsave("GSEA_TNFA_SIGNALING_VIA_NFKB.pdf", width = 12, height = 8)

