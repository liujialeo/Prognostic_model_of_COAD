rm(list = ls())
#install.packages("rms") 
#install.packages("foreign") 
#install.packages("survival") 
library(rms)
library(foreign)
library(survival)
library(tidyr)
library(dplyr)

load(file = "./output/COAD_clinical.Rdata")
rt=read.table("output/risk.txt",header=T,sep="\t")
rt <- rt[,c(1,17)]
colnames(rt)[1] <- "TCGA_id"

rt <- inner_join(rt,clinical, by = "TCGA_id")
rt <- rt[-c(502:503),]

rt$age <- factor(rt$age,labels=c())
rt$gender <- factor(rt$gender,labels=c("male", "female"))
rt$T <- factor(rt$T,labels=c("T1", "T2", "T3", "T4", "T4a","T4b","Tis"))
rt$N <- factor(rt$N,labels=c("N0", "N1", "N1a","N1b","N1c", "N2","N2a","N2b"))
rt$M <- factor(rt$M,labels=c("M0", "M1","M1a","M1b", "MX",""))
rt$risk <- factor(rt$risk,labels=c("high", "low"))
rt$stage <- factor(rt$stage,labels = c("not reported","stage i",
                                       "stage ia","stage ii","stage iia",
                                       "stage iib","stage iic","stage iii",
                                       "stage iiia","stage iiib","stage iiic",
                                       "stage iv","stage iva","stage ivb"))
ddist <- datadist(rt)
options(datadist='ddist')

#step1，筛选P<0.1者
f <- as.formula(Surv(futime,fustat) ~ age+gender+M+N+T+risk+stage)
cox <- coxph(f,data=rt)
summary(cox)

#step2，筛选P<0.05者
f1 <- as.formula(Surv(futime,fustat) ~ age+gender+M+N+T+risk+stage)
cox1 <- coxph(f1,data=rt)
summary(cox1)


#传统Nomogram
cox <- cph(Surv(futime,fustat) ~ T+N+risk,surv=T,x=T, y=T,data=rt) 
surv <- Survival(cox)
sur_1_year<-function(x)surv(1*365,lp=x)
sur_3_year<-function(x)surv(1*365*3,lp=x)
sur_5_year<-function(x)surv(1*365*5,lp=x)
nom_sur <- nomogram(cox,fun=list(sur_1_year,sur_3_year,sur_5_year),lp= F,funlabel=c('1-Year survival','3-Year survival','5-Year survival'),maxscale=100,fun.at= c('1.0',"0.99","0.95",'0.9','0.7','0.5','0.3','0.1',"0.01",'0'))
pdf("output/nomogram.pdf",width = 10, height = 5.5)
plot(nom_sur,xfrac=0.4)
dev.off()

