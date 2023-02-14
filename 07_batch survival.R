rm(list = ls())
##加载06中已经处理好的包含生存分析的表达数据
library(survival)
load(file = "exprSet_survival.Rda")
rt <- exprSet_survival
test <- rt[,1:10]
rownames(rt) <- rt[,1]
rt <- rt[,-1]
save(rt,file = "exprSet_batchsurvival.Rda")
###批量生存分析共有三种方法，是从果子老师的帖子中学到的：
###https://mp.weixin.qq.com/s?src=11&timestamp=1592568413&ver=2410&signature=PC4*LOG-Xv*L6qgTU0v5t6xkbF8PaDXI1kaho-EsVOUO0qzsVNjmBs2Hj45478yEnktMfDMcPkR8tbE6mlifPmQfUY4p161Jq0o8viuAgz6xv9mUGacrq5qstmkFlHE*&new=1

##我们采用的第一种方法，for循环
## 创建空的数据框
res2 <- data.frame()
## 获取基因列表
genes <- colnames(rt)[-c(1:2)]
## for循环开始
for (i in 1:length(genes)) {
  print(i)
  # 中位数分组
  group = ifelse(rt[,genes[i]]>median(rt[,genes[i]]),"high","low")
  if(length(table(group))==1) next
  surv =as.formula(paste('Surv(futime, fustat)~', "group"))
  data = cbind(rt[,1:2],group)
  #生存分析求差异
  x = survdiff(surv, data = data)
  # 获取p值
  pValue=1-pchisq(x$chisq,df=1) 
  #第一列基因名称
  res2[i,1] = genes[i]
  #第二列是p值
  res2[i,2] = pValue
}
names(res2) <- c("ID","pValue_log")
#作图
library(dplyr)
res2_significant <- res2 %>%
  filter(pValue_log < 0.05)%>%
  arrange(pValue_log)


rt <- as.data.frame(rt) ###下一步group容易报错，原因可能是rt不是dataframe，将其转化为数据框
index <- res2$ID[93]
group = ifelse(rt[,index]>median(rt[,index]),"high","low")
surv =as.formula(paste('Surv(futime, fustat)~', "group"))
data = cbind(rt[,1:2],group)
my.surv <- Surv(rt$futime, rt$fustat)
fit <- survfit(my.surv ~ group)
library(survminer)
ggsurvplot(fit, data = data, pval=TRUE,
           risk.table = TRUE,tables.height = 0.2,tables.theme = theme_cleantable(),
           ggtheme = theme_bw(), title = index)
##保存批量生存分析结果
save(res2, file = "batch_survival_results.Rda")
save(res2_significant,file = "batch_survival_results_sgnificant.Rda")



######批量生存分析作图
for (i in 1:nrow(res2)) {
  print(i)
  geneindex <- res2$ID[i]
  group = ifelse(rt[,index]>median(rt[,index]),"high","low")
  surv =as.formula(paste('Surv(futime, fustat)~', "group"))
  data = cbind(rt[,1:2],group)
  my.surv <- Surv(rt$futime, rt$fustat)
  fit <- survfit(my.surv ~ group)
  library(survminer)
  plot <- ggsurvplot(fit, data = data, pval=TRUE,
             risk.table = TRUE,tables.height = 0.2,tables.theme = theme_cleantable(),
             ggtheme = theme_bw(), title = index)
  plot
  ggsave(plot, filename = paste0("survival/survival_", geneindex, ".pdf"),width = 12, height = 8)
  
}





##将生存分析有统计学意义的基因和乳腺癌中高表达的基因取交集
rm(list = ls())
library(tidyr)
library(dplyr)
load(file = "batch_survival_results_sgnificant.Rda")
load(file = "diffLab.Rda")

diff_gene <- cbind(rownames(diffLab),diffLab)
colnames(diff_gene)[1] <- "ID"
diff_gene <- diff_gene %>%
  filter(logFC>1)%>%
  arrange(logFC)
diff_gene <- diff_gene[,c(1:2)]

upgene_survival <- inner_join(res2_significant,diff_gene,by = "ID")
save(upgene_survival,file = "upgene_survival.Rda")

##单独作图
load(file = "exprSet_batchsurvival.Rda")
rt <- as.data.frame(rt)
index <- upgene_survival$ID[1]
group = ifelse(rt[,index]>median(rt[,index]),"high","low")
surv =as.formula(paste('Surv(futime, fustat)~', "group"))
data = cbind(rt[,1:2],group)
my.surv <- Surv(rt$futime, rt$fustat)
fit <- survfit(my.surv ~ group)
library(survminer)
ggsurvplot(fit, data = data, pval=TRUE,
           risk.table = TRUE,tables.height = 0.2,tables.theme = theme_cleantable(),
           ggtheme = theme_bw() , title = index)


#######################################################################################
###将COAD和colitis的脂肪酸代谢相关基因交集进行批量生存分析作图
load(file = "innergene_eGSEA_fatty_acid.Rda")
##挑一个基因试一下

index <- innergene_eGSEA_fatty_acid$gene[2]
group = ifelse(rt[,index]>median(rt[,index]),"high","low")
surv =as.formula(paste('Surv(futime, fustat)~', "group"))
data = cbind(rt[,1:2],group)
my.surv <- Surv(rt$futime, rt$fustat)
fit <- survfit(my.surv ~ group)
library(survminer)
surplot <- ggsurvplot(fit, data = data, pval=TRUE,conf.int=T,  #置信区间
           risk.table = TRUE,tables.height = 0.2,tables.theme = theme_cleantable(),
           ggtheme = theme_bw() , title = index)
surplot

pdf(file = "survival_CA2.pdf" ,width=6,height=6,)
print(surplot)
dev.off()
##批量


###将高表达低生存的基因挑选出来,进行富集分析，找关键通路
rm(list = ls())
load(file = "upgene_survival.Rda")

upgene_survival_up <- as.data.frame(upgene_survival[,1]) 


upgene_survival_up <- data.table::fread(file = "survival_upgene_up.txt")
colnames(upgene_survival_up)[1] <- "ID"
upgene_survival_up <- inner_join(upgene_survival_up,upgene_survival,by = "ID")

rownames(upgene_survival_up) <- upgene_survival_up[,1]

gene <- rownames(upgene_survival_up)
#基因名称转换，返回的是数据框
library(clusterProfiler)
gene = bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(gene)


#**GO分析**

ego_CC <- enrichGO(gene = gene$ENTREZID,
                   OrgDb= org.Hs.eg.db,
                   ont = "all",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)


abc <- as.data.frame(ego_CC)

#**作图
p <- dotplot(ego_CC, split="ONTOLOGY") +facet_grid(ONTOLOGY~., scale="free")
p


#**KEGG分析**
kk <- enrichKEGG(gene = gene$ENTREZID,
                 organism ="human",
                 pvalueCutoff = 1,
                 qvalueCutoff = 1,
                 minGSSize = 1,
                 #readable = TRUE ,
                 use_internal_data =FALSE)



abc <- as.data.frame(kk)
dotplot(kk)

###富集到某一个通路上的图
browseKEGG(kk, 'hsa00100')
