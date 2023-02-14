### 更详细的GSEA，leading edge analysis
### xPierGSEA
rm(list = ls())

### 基于pi的GSEA分析
options("repos"=c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
if(!require("RCircos")) install.packages("RCircos",update = F,ask = F)
if(!require("ggnetwork")) install.packages("ggnetwork",update = F,ask = F)
if(!require("tibble")) install.packages("tibble",update = F,ask = F)
if(!require("gridExtra")) install.packages("gridExtra",update = F,ask = F)
install.packages("https://cran.r-project.org/src/contrib/Archive/XGR/XGR_1.1.7.tar.gz",repos=NULL,type="source")
if(!require("Pi")) BiocManager::install("Pi",update = F,ask = F)


### 数据是什么？
load(file = "allDiff.Rda")

library(dplyr)
library(tibble)
data <- allDiff %>%
  rownames_to_column("gene") %>% 
  select(gene,logFC,adj.P.Val) %>%
  mutate(
    offset = logFC,
    offset_abs = abs(logFC),
  ) %>% 
  mutate(rank = rank(-offset), priority = offset) %>%
  select(priority, rank, gene) %>% 
  column_to_rownames("gene")


# GSEA using MsigdbH
####################
### GSEA 分析
library(Pi)
eGSEA <-
  Pi::xPierGSEA(
    data,
    ontology = "MsigdbH",
    size.range = c(20, 5000),
    nperm = 20000,
    fast = F,
    RData.location = "http://galahad.well.ox.ac.uk/bigdata")


###发现连接不上上边的网址，需要本地下载放到指定文件夹那里
library(Pi)
index <- c("GWAS2EF", "GWAS_LD", "IlluminaHumanHT",
           "IlluminaOmniExpress", "ig.DO",
           "ig.EF", "ig.GOBP", "ig.GOCC", "ig.GOMF", "ig.HPCM", "ig.HPMA",
           "ig.HPMI", "ig.HPPA",
           "ig.MP", "org.Hs.eg", "org.Hs.egDGIdb", "org.Hs.egDO", "org.Hs.egGOBP",
           "org.Hs.egGOCC", "org.Hs.egGOMF", "org.Hs.egHPCM", "org.Hs.egHPMA",
           "org.Hs.egHPMI",
           "org.Hs.egHPPA", "org.Hs.egMP", "org.Hs.egMsigdbC1",
           "org.Hs.egMsigdbC2BIOCARTA",
           "org.Hs.egMsigdbC2CGP", "org.Hs.egMsigdbC2CPall",
           "org.Hs.egMsigdbC2CP",
           "org.Hs.egMsigdbC2KEGG", "org.Hs.egMsigdbC2REACTOME",
           "org.Hs.egMsigdbC3MIR",
           "org.Hs.egMsigdbC3TFT", "org.Hs.egMsigdbC4CGN", "org.Hs.egMsigdbC4CM",
           "org.Hs.egMsigdbC5BP", "org.Hs.egMsigdbC5CC", "org.Hs.egMsigdbC5MF",
           "org.Hs.egMsigdbC6", "org.Hs.egMsigdbC7", "org.Hs.egMsigdbH",
           "org.Hs.egPS",
           "org.Hs.egSF", "org.Hs.egPfam", "org.Hs.string", "org.Hs.PCommons_DN",
           "org.Hs.PCommons_UN")

###下载是由xRDataLoader这个函数执行的，创建一个文件夹叫作，“ontology_Rdata”用来存放数据，for循环执行操作
dir.create("ontology_Rdata/")
for (i in index) {
  print(i)
  geneset = xRDataLoader(RData=i)
  assign(i,geneset)
  save(list = i,file = paste0("ontology_Rdata/",i,".Rdata"))
}

###重新分析，把这个文件夹放在当前工作目录，修改RData.location参数就可以了
eGSEA <-
  Pi::xPierGSEA(
    data,
    ontology = "MsigdbH",
    size.range = c(20, 5000),
    nperm = 20000,
    fast = F,
    RData.location = paste0(getwd(),"/ontology_Rdata"))
save(eGSEA,file = "eGSEA.Rda")
####################
### 使用xGSEAdotplot来画图
### leading edge 分析
Pi::xGSEAdotplot(
  eGSEA,
  top = 1,
  leading = T,
  leading.edge.only = F,
  leading.size = 5,
  colormap = "spectral",
  zlim = c(-3, 3),
  peak = T,
  peak.color = 'black',
  clab = 'LFC\n(ipiNivo - pembro)',
  signature = FALSE) + theme(
    plot.title = element_text(hjust = 0.5, size = 8),
    plot.subtitle = element_text(hjust = 0.5, size = 8)
  )
##########################################################################################
##将脂肪酸代谢通路相关基因提取出来
fatty_acid_genes <- as.data.frame(eGSEA[["leading"]][["HALLMARK_FATTY_ACID_METABOLISM"]])
colnames(fatty_acid_genes)[1] <- "ENTREZID"
fatty_acid_genes_COAD <- fatty_acid_genes
save(fatty_acid_genes_COAD,file = "eGSEA_fatty_acid_genes_COAD.Rda")

load(file = "exprSet_diff.Rda")

if(! require("pheatmap")) install.packages("pheatmap")
制作一个分组信息
heatdata <- exprSet_diff [rownames(fatty_acid_genes),]
class <- c(rep("normal",41),rep("tumor",471)) 
annotation_col <- data.frame(class)
rownames(annotation_col) <- colnames(heatdata)

#如果注释出界，可以通过调整格子比例和字体修正
pheatmap(heatdata, #热图的数据
         cluster_rows = T,#行聚类
         cluster_cols = F,#列聚类，可以看出样本之间的区分度
         annotation_col =annotation_col, #标注样本分类
         annotation_legend=T, # 显示注释
         show_rownames = T,# 显示行名
         show_colnames = F, 
         border_color = "grey60",
         scale = "row", #以行来标准化，这个功能很不错
         color =colorRampPalette(c("blue", "white","red"))(100),#调色
         #filename = "heatmap_F.pdf",#是否保存
         cellwidth = 1, cellheight = 9,# 格子比例
         fontsize = 10)

#####################################################################################




### 批量作图画
ls_gp <- xGSEAdotplot(eGSEA, top=1:4, signature=F,colormap = "gbr",
                      subtitle ="both",leading = F)
library(gridExtra)
grid.arrange(grobs=ls_gp, ncol=1)

####################
### 画个火山图
df_summary <- eGSEA$df_summary %>%
  mutate(label = gsub('HALLMARK_', '', setID)) %>%
  mutate(flag = ifelse(adjp < 0.05, 'Y', 'N'))

### 火山图
ggplot(df_summary, aes(x = nes, y = -log10(adjp))) +
  geom_point(alpha = 0.6, shape = 16) +
  xlab("NES") +
  ylab(expression(-log[10]("FDR"))) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggrepel::geom_text_repel(
    data = subset(df_summary, flag == 'Y'),
    aes(label = label),
    size = 4,
    show.legend = F,
    segment.alpha = 1,
    segment.color = "RED",
    segment.size = 1,
    arrow = arrow(length = unit(0.01, 'npc'))
  )

### 使用火山图呈现GSEA富集分析的结果
### https://mp.weixin.qq.com/s/jolWmKoLic5m_M5F1E8K2g
### 使用举例
### https://www.nature.com/articles/s41591-019-0734-6