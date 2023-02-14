###单因素COX分析

rm(list = ls())
outTab=data.frame()

library(survival)
load(file = "exprSet_survival.Rda")
test <- exprSet_survival[1:10,1:10]
rownames(exprSet_survival) <- exprSet_survival[,1]
exprSet_survival <- exprSet_survival[,-1]
rt <- exprSet_survival

for(i in colnames(rt[,3:ncol(rt)])){
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  outTab=rbind(outTab,cbind(gene=i,HR=coxSummary$coefficients[,"exp(coef)"],
                            z=coxSummary$coefficients[,"z"],
                            pvalue=coxSummary$coefficients[,"Pr(>|z|)"]))
}

save(outTab,file = "single_cox.Rda")

##多因素COX分析
