rm(list = ls())
#install.packages("survival")
#install.packages("survminer")

library(survival)
library("survminer")
rt=read.table("output/risk.txt",header=T,sep="\t")
diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
pdf(file="output/survival.pdf",onefile = FALSE,
    width = 8,
    height =8)
ggsurvplot(fit, 
           data=rt,
           conf.int=TRUE,
           pval=paste0("p=",pValue),
           pval.size=4,
           risk.table=T,
           legend.labs=c("High risk", "Low risk"),
           legend.title="Risk",
           xlab="Time(days)",
           break.time.by = 365,
           ggtheme = theme_light(),
           risk.table.y.text.col = T,
           risk.table.height = 0.18,
           risk.table.y.text = F,
           ncensor.plot = T,
           ncensor.plot.height = 0.18,
           conf.int.style = "ribbon")
dev.off()
#summary(fit)
