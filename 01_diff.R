############差异表达分析
##构建分组矩阵,这一段有两种方法，但是初学的时候特别容易误解，这一步分完全就是独立的
rm(list = ls())
library(limma)
load(file = "exprSet_diff.Rda")

class <- c(rep("normal",41),rep("tumor",471)) 
design <- model.matrix(~factor(class))
colnames(design) <- c("normal","tumor")
design

#线性模型拟合
fit <- lmFit(exprSet_diff,design)
#贝叶斯检验
fit2 <- eBayes(fit)
#输出基因
allDiff=topTable(fit2,adjust='fdr',coef=2,number=Inf) 
save(allDiff,file = "allDiff.Rda")

#写入表格
write.table(allDiff,file="limmaTab.xls",sep="\t",quote=F)

#找出差异两倍以上，pvalue小于0.05，1078个
diffLab <- subset(allDiff,(logFC> 1 | logFC< (-1)) & adj.P.Val < 0.05)
save(diffLab,file = "diffLab.Rda")
write.table(diffLab,file="diffExp.xls",sep="\t",quote=F)
