rm(list = ls())
library(dplyr)
load(file = "output/rt_plot.Rdata")
library(dplyr)
rt <- rt_plot %>% 
  # 去掉小于30天的
  filter(futime >= 30) %>% 
  mutate(futime = futime/365)

### 我们做单个基因的生存分析，这里是单个甲基化位点UCN
geneindex <- "UCN"
rt <- rt[,c("futime","fustat",geneindex)]
save(rt,file = "output/rt.Rdata")
#ifelse 联合median快速二分类
rt$risk <- ifelse(rt[,geneindex] > median(rt[,geneindex]),"High","Low")

### logrank的方法
### 首先Surv函数用于创建生存数据对象
library(survival)
surv_object = Surv(rt$futime, rt$fustat)
## 生存数据拟合survfit
fit1 <- survfit(surv_object ~ risk, data = rt)
summary(fit1)
library(survminer)
p <- ggsurvplot(fit1,
           data=rt,
           surv.median.line = "hv",  # 添加中位生存时间线
           legend = c(0.8,0.85), # 指定图例位置
           conf.int=TRUE,
           pval=TRUE,
           pval.size=6,
           risk.table=F,
           legend.labs=c("High", "Low"),
           legend.title=paste0(geneindex," expression"),
           xlab="Time(years)",
           break.time.by = 1,
           ggtheme = theme_light(),
           risk.table.y.text.col = T,
           risk.table.height = 0.2,
           risk.table.y.text = F,
           ncensor.plot = F,
           ncensor.plot.height = 0.2,
           conf.int.style = "ribbon")
p
library(export)
graph2pdf(file = paste0("./output/","survival_",geneindex,".pdf"),width = 4,height = 3)

### 这个里面的是pval如何获取呢？survdiff函数
x = survdiff(surv_object ~ risk, data = rt)
pValue=1-pchisq(x$chisq,df=1)
pValue= round(pValue,3)
pValue

###################################################
### 有没有更好的策略？化腐朽为神奇！
### Best separation
rm(list = ls())
### 导入上一步保存的rt
load(file = "output/rt.Rdata")
library(survival)
library(survminer)
res.cut <- surv_cutpoint(rt, 
                         time ="futime", 
                         event = "fustat", 
                         variables = "UCN", 
                         minprop = 0.3)  
##按照bestSeparation分高低表达
risk <- surv_categorize(res.cut)
rt$risk <- risk$UCN
## 构建生存对象Surv
surv_object = Surv(rt$futime, rt$fustat)
## 生存数据拟合survfit
fit1 <- survfit(surv_object ~ risk, data = rt)
summary(fit1)
geneindex <- "UCN"
ggsurvplot(fit1, data = rt, pval = TRUE)

p <- ggsurvplot(fit1,
           data=rt,
           surv.median.line = "hv",  # 添加中位生存时间线
           legend = c(0.8,0.85), # 指定图例位置
           conf.int=TRUE,
           pval=TRUE,
           pval.size=6,
           risk.table=T,
           legend.labs=c("High", "Low"),
           legend.title=paste0(geneindex," expression"),
           xlab="Time(days)",
           break.time.by = 1,
           ggtheme = theme_light(),
           risk.table.y.text.col = T,
           risk.table.height = 0.2,
           risk.table.y.text = F,
           ncensor.plot = T,
           ncensor.plot.height = 0.2,
           conf.int.style = "ribbon")
p

library(export)
graph2pdf(file = paste0("./output/","survival_",geneindex,"_best.pdf"),width = 8,height = 6)

### 这个里面的是pval如何获取呢？survdiff
if(T){
  x = survdiff(surv_object ~ risk, data = rt)
  pValue=1-pchisq(x$chisq,df=1)
  pValue= round(pValue,4)
  pValue
}

