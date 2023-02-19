rm(list = ls())
library("glmnet")
library("survival")

rt=read.table("output/diffgene_uniSigExp.txt",header=T,sep="\t",row.names=1,check.names=F)
rt$futime[rt$futime<=0]=1
rt$futime=rt$futime/1
genes=colnames(rt)
# gene=grep(paste0("\\|",""),genes,value=T)
# geneLength=ifelse(length(gene)>20,20,length(gene))
# rt=rt[,c("futime","fustat",gene[1:geneLength])]

x=as.matrix(rt[,c(3:ncol(rt))])
y=data.matrix(Surv(rt$futime,rt$fustat))

fit <- glmnet(x, y, family = "cox", maxit = 1000)
pdf("output/lambda.pdf")
plot(fit, xvar = "lambda", label = TRUE)
dev.off()

cvfit <- cv.glmnet(x, y, family="cox", maxit = 1000)
pdf("output/cvfit.pdf")
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
dev.off()

coef <- coef(fit, s = cvfit$lambda.min)
index <- which(coef != 0)
actCoef <- coef[index]
lassoGene=row.names(coef)[index]
lassoGene=c("futime","fustat",lassoGene)
lassoSigExp=rt[,lassoGene]
lassoSigExp=cbind(id=row.names(lassoSigExp),lassoSigExp)
write.table(lassoSigExp,file="output/diffgene_lassoSigExp.txt",sep="\t",row.names=F,quote=F)
