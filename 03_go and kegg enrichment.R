#########
# GO分析，KEGG分析
###########################
###########################
rm(list = ls())
library(clusterProfiler)

#获得基因列表
load(file = "diffLab.Rda")
gene <- rownames(diffLab)
#基因名称转换，返回的是数据框
gene = bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(gene)


#**GO分析**

ego_CC <- enrichGO(gene = gene$ENTREZID,
                   OrgDb= org.Hs.eg.db,
                   ont = "all",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.01,
                   readable = TRUE)

save(ego_CC, file = "go_enrichment.Rda")
load(file = "go_enrichment.Rda")
abc <- as.data.frame(ego_CC)

#**作图
p <- dotplot(ego_CC, split="ONTOLOGY") +facet_grid(ONTOLOGY~., scale="free")
p


#**KEGG分析**
kk <- enrichKEGG(gene = gene$ENTREZID,
                 organism ="human",
                 pvalueCutoff = 0.01,
                 qvalueCutoff = 0.01,
                 minGSSize = 1,
                 #readable = TRUE ,
                 use_internal_data =FALSE)


save(kk,file = "kegg_enrichment.Rda")
load(file = "kegg_enrichment.Rda")
abc <- as.data.frame(kk)
dotplot(kk)
###将图导出为PPT格式
library(export)
graph2ppt(file="KEGG.pptx", width=7, height=5)

###富集到某一个通路上的图
browseKEGG(kk, 'hsa04151')
