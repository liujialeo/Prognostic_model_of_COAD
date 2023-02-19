rm(list = ls())
library(ggplot2)
#install.packages("survival")
#install.packages("survminer")
library(tidyr)
library(dplyr)
library(survival)
library("survminer")

rt=read.table("output/risk.txt",header=T,sep="\t")
colnames(rt)[1] <- "TCGA_id"
load(file = "./output/COAD_clinical.Rdata")
clinical <- clinical[,c(1,4,5)]
rt <- inner_join(clinical,rt, by = "TCGA_id")


s = as.formula(paste("Surv(futime, fustat) ~ ",paste(colnames(rt)[c(3,6:17)],collapse = "+")))
model <- coxph(s, data = rt )


options(scipen=1)
p <- ggforest(model, data =rt, 
         main = "Hazard ratio", 
         cpositions = c(0.10, 0.22, 0.4), 
         fontsize = 1.0, 
         refLabel = "1", noDigits = 4)
p
ggsave(p,file = "output/multivariateCOX_forest.pdf")
