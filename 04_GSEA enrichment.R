####GSEA分析来啦

#获得基因列表
rm(list = ls())
library(clusterProfiler)
load(file = "allDiff.Rda")
eg <- rownames(allDiff)
#基因名称转换，返回的是数据框
eg = bitr(eg, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

eg <- dplyr::distinct(eg,SYMBOL,.keep_all=TRUE)
# geneList 三部曲
geneList = allDiff[,1]
names(geneList) = eg$ENTREZID
geneList = sort(geneList, decreasing = TRUE)
geneList[1:5]

library(clusterProfiler)
library(export)
install.packages("R.utils")
R.utils::setOption("clusterProfiler.download.method",'auto') 
############# KEGG
gseaKEGG <- gseKEGG(geneList     = geneList,
                    organism     = "hsa",
                    nPerm        = 1000,
                    minGSSize    = 20,
                    pvalueCutoff = 1,
                    verbose      = FALSE)
abc <- as.data.frame(gseaKEGG)
library(ggplot2)
dotplot(gseaKEGG,showCategory=12,split=".sign")+facet_grid(~.sign)
test <- as.data.frame(gseaKEGG)
graph2ppt(file="GSEAKEGG.pptx", width=7, height=5)
################单个通路作图

###1.中心碳代谢
library(enrichplot)
pathway.id = "hsa04668"
gseaplot2(gseaKEGG, 
          color = "green",
          geneSetID = pathway.id,
          pvalue_table = T)
###pathview
library(pathview)
pathway.id = "hsa04668"
pv.out <- pathview(gene.data  = geneList,
                   pathway.id = pathway.id,
                   species    = "hsa")
library(export)
graph2ppt(file="TNF signaling pathway.pptx", width=7, height=5)


###另一种方式的pathview
library(pathview)
pathway.id = "hsa05230"
pv.out <- pathview(gene.data  = geneList,
                   pathway.id = pathway.id,
                   species    = "hsa",
                   kegg.native = F)


############# GO
require(DOSE)
library(clusterProfiler)
gseaGO <- gseGO(geneList     = geneList,
                    OrgDb        = org.Hs.eg.db,
                    ont          = "all",
                    nPerm        = 1000,
                    minGSSize    = 20,
                    pvalueCutoff = 1,
                    verbose      = FALSE)
abc <- as.data.frame(gseaGO)
library(ggplot2)
dotplot(gseaGO,showCategory=12,split=".sign")+facet_grid(~.sign)
test <- as.data.frame(gseaGO)

