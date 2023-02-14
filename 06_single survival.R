rm(list = ls())

##读取生存数据
clin = data.table::fread(file = "TCGA-COAD.survival.tsv")
colnames(clin) <- c("TCGA_id","fustat","sample","futime")
clin <- clin[,c(1,4,2)]
##载入表达数据
load(file = "exprSet_samples.Rda")
test <- exprSet_merge[1:10,1:10]

##合并表达数据和生存数据
library(tidyr)
library(dplyr)
exprSet_survival <- inner_join(clin,exprSet_merge,by = "TCGA_id")
test <- exprSet_survival[1:10,1:10]
exprSet_survival <- arrange(exprSet_survival, sample)

##将非肿瘤组织去掉
table(exprSet_survival$sample) ##共有39个正常样本，去掉
test <- exprSet_survival[,1:10]
exprSet_survival <- exprSet_survival[-(1:39),]
exprSet_survival <- exprSet_survival[,-4] ##去掉sample一列
save(exprSet_survival,file = "exprSet_survival.Rda")
##数据格式已经处理成为做生存分析所用的格式
##下面进行生存分析
library(survival)
library(survminer)
library(ggplot2)
##TNF signaling pathway 相关基因IL1A","TNF","IL6","TNFAIP3","NFKB1","REL","NFKBIE","CXCL10","TLR2","ICAM1"
group <- ifelse(exprSet_survival$ICAM1>median(exprSet_survival$ICAM1),'high','low')  ##修改基因名称
sfit <- survfit(Surv(futime, fustat)~group, data=exprSet_survival)
sfit
summary(sfit)
p <- ggsurvplot(sfit, conf.int=F, pval=TRUE,
                surv.median.line = "hv",  # 添加中位生存时间线
                xlab = "Follow up time(d)", # 指定x轴标签
                legend = c(0.8,0.9), # 指定图例位置
                legend.labs = c("High", "Low"), # 指定图例分组标签
                risk.table = TRUE,tables.height = 0.2,tables.theme = theme_cleantable(),
                ggtheme = theme_bw() )
p
ggsave(p,filename = "ICAM1_survival.pdf",width = 8,height = 6)

library(export)
graph2pdf(file = "ICAM1_survival.pdf",width = 8,height = 6)
graph2tif(file = "ICAM1_survival.tiff",width = 8,height = 6)
graph2ppt(file="ICAM1_survival.pptx", width=7, height=5)



###################################################
### 有没有更好的策略？化腐朽为神奇！
### Best separation
rm(list = ls())
### 导入上一步保存的rt
load(file = "exprSet_survival.Rda")
rt <- exprSet_survival
library(survival)
library(survminer)
res.cut <- surv_cutpoint(rt, 
                         time ="futime", 
                         event = "fustat", 
                         variables = "ICAM1", 
                         minprop = 0.3)  
##按照bestSeparation分高低表达
risk <- surv_categorize(res.cut)
rt$risk <- risk$ICAM1
## 构建生存对象Surv
surv_object = Surv(rt$futime, rt$fustat)
## 生存数据拟合survfit
fit1 <- survfit(surv_object ~ risk, data = rt)
summary(fit1)
ggsurvplot(fit1, , conf.int=F, pval=TRUE,
           surv.median.line = "hv",  # 添加中位生存时间线
           xlab = "Follow up time(d)", # 指定x轴标签
           legend = c(0.8,0.9), # 指定图例位置
           legend.labs = c("High", "Low"), # 指定图例分组标签
           risk.table = TRUE,tables.height = 0.2,tables.theme = theme_cleantable(),
           ggtheme = theme_bw() )
### 这个里面的是pval如何获取呢？survdiff
if(T){
  x = survdiff(surv_object ~ risk, data = rt)
  pValue=1-pchisq(x$chisq,df=1)
  pValue= round(pValue,4)
  pValue
}
library(export)
graph2pdf(file = "ICAM1_survival_best.pdf",width = 8,height = 6)

