rm(list = ls())
library(survival)

rt=read.table("output/lassoSigExp.txt",header=T,sep="\t",check.names=F,row.names=1)
rt <- rt[1:200,]
sample_id <- as.data.frame(row.names(rt))
colnames(sample_id)[1] <- "TCGA_id"

genes <- c("futime","fustat","FAM132B","DLX4","UCN","GABRD","MS4A2","DYNC1I1",
           "GRIK3","VWC2","LEP","SNCB","CDH10","GABRG1")
rt <- rt[,genes]

cox_m <- coxph(Surv(futime,fustat) ~FAM132B+DLX4+UCN+GABRD+MS4A2+DYNC1I1+
                 GRIK3+VWC2+LEP+SNCB+CDH10+GABRG1, data = rt)
cox_m1<-step(cox_m,direction = "both")
riskScore<-predict(cox_m1,type="risk",newdata=rt)
risk<-as.vector(ifelse(riskScore>median(riskScore),"high","low"))
dir.create("test_output")
write.table(cbind(id=rownames(cbind(rt[,1:14],riskScore,risk)),cbind(rt[,1:14],riskScore,risk)),"test_output/risk.txt",sep="\t",quote=F,row.names=F)
save(sample_id,file = "./test_output/sample_id.Rdata")
