rm(list = ls())
library(dplyr)
library(tidyr)
library(tibble)
###从UCSC-XENA网站下载的乳腺癌转录组数据，数据读入
exprSet_FPKM = data.table::fread(file = "TCGA-COAD.htseq_fpkm.tsv")
test <- exprSet_FPKM[1:100,1:100]
colnames(exprSet_FPKM)[1] <- "gene_id"

#将第一列的gene_id中的点去掉，用的是separate函数
exprSet_FPKM <- exprSet_FPKM%>%
  separate(col = gene_id, into = c("id","dot"))%>%
  select(-dot)

test2 <- exprSet_FPKM[1:100,1:100]
colnames(exprSet_FPKM)[1] <- "gene_id"  
  
  
##ID转换，提取编码RNA，注释基因
load(file = "gtf_df.Rda")

ABC <- gtf_df %>%
  filter(type=="gene",gene_biotype=="protein_coding")
ABC <- ABC[,c(10,12)]
##将表达数据跟上边的ID转换文件进行合并
exprSet <- inner_join(ABC,exprSet_FPKM,by = "gene_id")
test <- exprSet[1:10,1:10]
exprSet <- exprSet %>%
  distinct(gene_name,.keep_all = T)
exprSet <- exprSet[,-1]
rownames(exprSet) <- exprSet[,1]
exprSet <- exprSet[,-1]
range(exprSet)
##保存处理过的数据
save(exprSet,file = "exprSet.Rda")


rm(list = ls())
library(tidyr)
library(dplyr)
library(limma)
load(file = "exprSet.Rda")
test <- exprSet[1:10,1:10]
##将癌组织和癌旁组织进行标注
metadata <- data.frame(names(exprSet)) #转换成数据框
for (i in 1:length(metadata[,1])) {
  num <- as.numeric(substring(metadata[i,1],14,15)) #substring的用法，这是元素获取
  if (num %in% seq(1,9)) {metadata[i,2] <- "Tumor"} #如果是肿瘤，就给第2列加上Tumor
  if (num %in% seq(10,29)) {metadata[i,2] <- "Normal"} #如果是正常组织，就给第2列加上Normal
}
metadata <-  arrange(metadata, V2)
names(metadata) <- c("TCGA_id","sample")
##将metadata文件和乳腺癌表达数据进行合并
#首先，将表达数据进行转置，列名变行名
exprSet_t <- as.data.frame(t(exprSet)) 
test <- exprSet_t[1:10,1:10]
exprSet_t <- cbind(rownames(exprSet_t),exprSet_t)
colnames(exprSet_t)[1] <- "TCGA_id"

##合并样本信息和表达数据
exprSet_merge <- inner_join(metadata,exprSet_t, by = "TCGA_id")
test <- exprSet_merge[,1:10]
##计算癌组织和癌旁样本的个数
table(exprSet_merge$sample) ##癌旁41,癌组织471
save(exprSet_merge, file = "exprSet_samples.Rda")
##再将sample信息删除，然后再转置
exprSet_diff <- exprSet_merge[,-2]
test <- exprSet_diff[1:10,1:10]
rownames(exprSet_diff) <- exprSet_diff[,1]
exprSet_diff <- exprSet_diff[,-1]
exprSet_diff <- as.data.frame(t(exprSet_diff))

save(exprSet_diff,file = "exprSet_diff.Rda")



