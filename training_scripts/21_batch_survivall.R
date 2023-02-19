### 本节任务: 批量生存分析
#############################################################################=
rm(list = ls())
load("output/rt_plot.Rdata")
test <- rt_plot[1:10,1:10]
library(dplyr)
coxdata <- rt_plot %>% 
  # 去掉小于30天的
  filter(futime >= 30) %>% 
  mutate(futime = futime/365)
test <- coxdata[,1:10]
## 单因素cox分析，正常情况下这里所有的因素都会小于0.05
## gsub
colnames(coxdata) <- gsub("-","_",colnames(coxdata))
colnames(coxdata) <- gsub(":","_",colnames(coxdata))
genes <- colnames(coxdata)[-c(1:2)][1:1000]
library(survival)
res <- data.frame()
for (i in 1:length(genes)) {
  print(i)
  surv =as.formula(paste('Surv(futime, fustat)~', genes[i]))
  x = coxph(surv, data = coxdata)
  x = summary(x)
  p.value=signif(x$wald["pvalue"], digits=2)
  HR =signif(x$coef[2], digits=2);#exp(beta)
  HR.confint.lower = signif(x$conf.int[,"lower .95"], 2)
  HR.confint.upper = signif(x$conf.int[,"upper .95"],2)
  CI <- paste0("(", 
               HR.confint.lower, "-", HR.confint.upper, ")")
  res[i,1] = genes[i]
  res[i,2] = HR
  res[i,3] = CI
  res[i,4] = p.value
}
names(res) <- c("ID","HR","95% CI","p.value")
save(res,file = "output/res_HR.Rdata")

### 批量生存分析log-rank
rm(list = ls())
library(survival)
library(survminer)
load("output/rt_plot.Rdata")
library(dplyr)
rt <- rt_plot %>% 
  filter(futime >= 30) %>% # 去掉小于30天的
  mutate(futime = futime/365)

res.cut <- surv_cutpoint(rt, 
                         time = "futime", 
                         event = "fustat", 
                         variables = names(rt)[3:ncol(rt)], 
                         minprop = F) 
res.cat <- surv_categorize(res.cut)

test <- res.cat[1:10,1:10]
### 也可以使用median来批量做
colnames(rt) <- gsub("-","_",colnames(rt))
colnames(rt) <- gsub(":","_",colnames(rt))
genes <- colnames(rt)[-c(1:2)][1:1000]

res2 <- data.frame()
for (i in 1:length(genes)) {
  print(i)
  surv =as.formula(paste('Surv(futime, fustat)~', genes[i]))
  x = survdiff(surv, data = res.cat)
  pValue=1-pchisq(x$chisq,df=1)
  res2[i,1] = genes[i]
  res2[i,2] = pValue
}
names(res2) <- c("ID","pValue_log")
load(file = "output/res_HR.Rdata")
res1 <- res
res <- merge(res1,res2,by="ID")

### 联合筛选
library(dplyr)
res_filter <- res %>% 
  filter(p.value < 0.01) %>% 
  filter(pValue_log < 0.05)
## 最终得到少量位点，进行下一步筛选。
save(res_filter,file = "output/res_filter.Rdata")

##################################################################################
##突发状况如何处理???
### 批量生存分析log-rank
rm(list = ls())
if(T){
  load("resources/dd1.Rda")
  library(survival)
  library(survminer)
  library(dplyr)
  rt <- dd1 %>% 
    filter(futime >= 30) %>% # 去掉小于30天的
    mutate(futime = futime/365)
  test <- rt[1:10,1:10]
  #colnames(rt) <- gsub("-","_",colnames(rt))
  #colnames(rt) <- gsub(":","_",colnames(rt))
}

genes <- colnames(rt)[-c(1:2)][1:2000]
res2 <- data.frame()
for (i in 1:length(genes)) {
  print(i)
  surv =as.formula(paste('Surv(futime, fustat)~', "group"))
  group = ifelse(rt[,genes[i]]>median(rt[,genes[i]]),"high","low")
  #if(length(table(group))==1) next
  data = cbind(rt[,1:2],group)
  x = survdiff(surv, data = data)
  pValue=1-pchisq(x$chisq,df=1) 
  res2[i,1] = genes[i]
  res2[i,2] = pValue
}
names(res2) <- c("ID","pValue_log")

#######################
### 去掉表达量均一的基因
dd <- rt[,-c(1,2)]
dd1 <- dd[,apply(dd,2,var)!=0]

rt <- cbind(rt[,c(1,2)],dd1)

## 8秒完成2万个基因的生存分析，人人都可以！
## https://mp.weixin.qq.com/s/o4e1HzG4zPIQoGT6-7D0ug