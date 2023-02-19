rm(list = ls())
#install.packages("rms") 
#install.packages("foreign") 
#install.packages("survival") 
library(rms)
library(foreign)
library(survival)

load(file = "./output/COAD_clinical.Rdata")
rt=read.table("test_output/risk.txt",header=T,sep="\t")
rt <- rt[,c(1,17)]
colnames(rt)[1] <- "TCGA_id"

rt <- inner_join(rt,clinical, by = "TCGA_id")

colnames(rt)[2] <- "risk_score"


rt$age <- factor(rt$age,labels=c())
rt$gender <- factor(rt$gender,labels=c("male", "female"))
rt$T <- factor(rt$T,labels=c("T1", "T2", "T3", "T4","T4a","T4b"))
rt$N <- factor(rt$N,labels=c("N0", "N1", "N1a","N1b", "N2","N2a","N2b"))
rt$M <- factor(rt$M,labels=c("M0", "M1","M1a", "MX",""))
rt$risk <- factor(rt$risk,labels=c("high", "low"))
rt$stage <- factor(rt$stage,labels = c("not reported","stage i","stage ii","stage iia",
                                       "stage iib","stage iii",
                                       "stage iiia","stage iiib","stage iiic",
                                       "stage iv","stage iva"))
ddist <- datadist(rt)
options(datadist='ddist')



#1-year
cox1 <- cph(Surv(futime,fustat) ~ T+M+N+risk_score,surv=T,x=T, y=T,time.inc = 1*365,data=rt) 
cal <- calibrate(cox1, cmethod="KM", method="boot", u=1*365, m= 60, B=1000)  ##样本数变了，m值也要变，应该小于200/3，改为60
pdf("test_output/calibration1.pdf",10,8)
par(mar = c(10,5,3,2),cex = 1.0)
plot(cal,lwd=3,lty=2,errbar.col="black",xlim = c(0,1.0),ylim = c(0,1.0),xlab ="Nomogram-Predicted Probability of 1-Year Survival",ylab="Actual 1-Year Survival",col="blue")
lines(cal,c('mean.predicted','KM'),type = 'a',lwd = 3,col ="black" ,pch = 16)
box(lwd = 1)
abline(0,1,lty = 3,lwd = 3,col = "black")
dev.off()

#3-year
cox1 <- cph(Surv(futime,fustat) ~ T+M+N+risk_score,surv=T,x=T, y=T,time.inc = 3*365,data=rt) 
cal <- calibrate(cox1, cmethod="KM", method="boot", u=3*365, m= 60, B=1000)
pdf("test_output/calibration3.pdf",10,8)
par(mar = c(10,5,3,2),cex = 1.0)
plot(cal,lwd=3,lty=2,errbar.col="black",xlim = c(0,1.0),ylim = c(0,1.0),xlab ="Nomogram-Predicted Probability of 3-Year Survival",ylab="Actual 3-Year Survival",col="blue")
#plot(cal,add=T, lwd=3,lty=2,errbar.col="red",xlim = c(0,1.0),ylim = c(0,1.0),xlab ="Nomogram-Predicted Probability of 3-Year Survival",ylab="Actual 3-Year Survival",col="blue")
lines(cal,c('mean.predicted','KM'),type = 'a',lwd = 3,col ="black" ,pch = 16)
box(lwd = 1)
abline(0,1,lty = 3,lwd = 3,col = "black")
dev.off()

#5-year
cox1 <- cph(Surv(futime,fustat) ~ T+M+N+risk_score,surv=T,x=T, y=T,time.inc = 5*365,data=rt) 
cal <- calibrate(cox1, cmethod="KM", method="boot", u=5*365, m= 60, B=1000)
pdf("test_output/calibration5.pdf",10,8)
par(mar = c(10,5,3,2),cex = 1.0)
plot(cal,lwd=3,lty=2,errbar.col="black",xlim = c(0,1.0),ylim = c(0,1.0),xlab ="Nomogram-Predicted Probability of 5-Year Survival",ylab="Actual 5-Year Survival",col="blue")
lines(cal,c('mean.predicted','KM'),type = 'a',lwd = 3,col ="black" ,pch = 16)
box(lwd = 1)
abline(0,1,lty = 3,lwd = 3,col = "black")
dev.off()
