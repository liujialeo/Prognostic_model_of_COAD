rm(list = ls())
#install.packages("survivalROC")

library(survivalROC)
rt=read.table("output/risk.txt",header=T,sep="\t",check.names=F,row.names=1)

#1年ROC
pdf(file="output/ROC-1.pdf",width=6,height=6)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, 
                predict.time =365, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
     xlab="1-Specificity", ylab="Sensitivity",
     main=paste("ROC curve (", "AUC = ",round(roc$AUC,3),")"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)
dev.off()

#3年ROC
pdf(file="output/ROC-3.pdf",width=6,height=6)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, 
                predict.time =365*3, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
     xlab="1-Specificity", ylab="Sensitivity",
     main=paste("ROC curve (", "AUC = ",round(roc$AUC,3),")"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)
dev.off()

#5年ROC
pdf(file="output/ROC-5.pdf",width=6,height=6)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, 
                predict.time =365*5, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
     xlab="1-Specificity", ylab="Sensitivity",
     main=paste("ROC curve (", "AUC = ",round(roc$AUC,3),")"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)
dev.off()

#整合1，3，5年ROC
pdf(file="output/ROC.pdf",width=6,height=6)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)

roc1=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, 
                 predict.time =365, method="KM")
plot(roc1$FP, roc1$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
     xlab="1-Specificity", ylab="Sensitivity",
     main=paste("ROC curve"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)

roc2=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, 
                 predict.time =365*3, method="KM")   #在此更改时间，单位为年
lines(roc2$FP,roc2$TP,type="l",xlim=c(0,1),ylim=c(0,1),col="blue",lwd=2)
#text(locator(1), paste("1 year",round(roc2$AUC,3),sep=":"),col="blue")

roc3=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, 
                 predict.time =365*5, method="KM")   #在此更改时间，单位为年
lines(roc3$FP,roc3$TP,type="l",xlim=c(0,1),ylim=c(0,1),col="green",lwd=2)
#text(locator(1), paste("2 year",round(roc2$AUC,3),sep=":"),col="green")

legend("bottomright", 
       c("1-year AUC:0.757","3-year AUC:0.726","5-year AUC:0.742"),
       lwd=2,
       col=c("red","blue","green"))
dev.off()

