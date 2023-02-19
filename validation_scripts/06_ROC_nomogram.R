rm(list = ls())
library(survivalROC)
library(survival)
##合并手动修改后的临床信息和多因素COX
dataclinical <- data.table::fread(file = "./output/dataclinical.csv",data.table = F) #保存csv格式时第一列有信息丢失
colnames(dataclinical)[1] <- "TCGA_id"

load(file = "./test_output/sample_id.Rdata")
dataclinical <- inner_join(sample_id,dataclinical, by = "TCGA_id")
rownames(dataclinical) <- dataclinical[,1]

rt=read.table("test_output/risk.txt",header=T,sep="\t") ##读入多因素COX分析结果，将riskScore作为一种临床信息
colnames(rt)[1] <- "TCGA_id"

rt <- rt[,c(1,17)]
dataclinical <- inner_join(dataclinical,rt,by = "TCGA_id")
rownames(dataclinical) <- dataclinical[,1]
dataclinical <- dataclinical[,-1]
## risk中low改为0,high改为1， for 循环解决问题
for (i in 1:nrow(dataclinical)) {
  if(dataclinical$risk[i]=="low"){
    dataclinical$risk[i]= 0
  }else{
    dataclinical$risk[i]= 1
  }
}
rt <- dataclinical
colnames(rt)[9] <- "risk_score"

#计算Nomogram模型的risk score，并分为高低危组
cox_m <- coxph(Surv(futime,fustat) ~age+T+M+N+risk_score, data = rt)
cox_m1<-step(cox_m,direction = "both")
risk_score<-predict(cox_m1,type="risk",newdata=rt)
risk_level<-as.vector(ifelse(risk_score>median(risk_score),"High","Low"))
write.table(cbind(id=rownames(cbind(rt[,1:2],risk_score,risk_level)),cbind(rt[,1:2],risk_score,risk_level)),"test_output/risk_score.txt",sep="\t",quote=F,row.names=F)

#绘制ROC曲线
rt=read.table("test_output/risk_score.txt",header=T,sep="\t",check.names=F,row.names=1)
#1年ROC
pdf(file="test_output/nomogram_ROC-1.pdf",width=6,height=6)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$risk_score, 
                predict.time =365, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
     xlab="1-Specificity", ylab="Sensitivity",
     main=paste("ROC curve (", "AUC = ",round(roc$AUC,3),")"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)
dev.off()

#3年ROC
pdf(file="test_output/nomogram_ROC-3.pdf",width=6,height=6)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$risk_score, 
                predict.time =365*3, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
     xlab="1-Specificity", ylab="Sensitivity",
     main=paste("ROC curve (", "AUC = ",round(roc$AUC,3),")"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)
dev.off()

#5年ROC
pdf(file="test_output/nomogram_ROC-5.pdf",width=6,height=6)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$risk_score, 
                predict.time =365*5, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
     xlab="1-Specificity", ylab="Sensitivity",
     main=paste("ROC curve (", "AUC = ",round(roc$AUC,3),")"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)
dev.off()

#整合1，3，5年ROC
pdf(file="test_output/nomogram_ROC.pdf",width=6,height=6)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)

roc1=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$risk_score, 
                 predict.time =365, method="KM")
plot(roc1$FP, roc1$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
     xlab="1-Specificity", ylab="Sensitivity",
     main=paste("ROC curve"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)

roc2=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$risk_score, 
                 predict.time =365*3, method="KM")   #在此更改时间，单位为年
lines(roc2$FP,roc2$TP,type="l",xlim=c(0,1),ylim=c(0,1),col="blue",lwd=2)
#text(locator(1), paste("1 year",round(roc2$AUC,3),sep=":"),col="blue")

roc3=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$risk_score, 
                 predict.time =365*5, method="KM")   #在此更改时间，单位为年
lines(roc3$FP,roc3$TP,type="l",xlim=c(0,1),ylim=c(0,1),col="green",lwd=2)
#text(locator(1), paste("2 year",round(roc2$AUC,3),sep=":"),col="green")

legend("bottomright", 
       c("1-year AUC:0.781","3-year AUC:0.765","5-year AUC:0.861"),
       lwd=2,
       col=c("red","blue","green"))
dev.off()
