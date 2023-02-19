rm(list = ls())
#install.packages('survival')
library(survival)
library(tidyr)
library(dplyr)

load(file = "./output/COAD_allDiff.Rdata")
load(file = "./output/COAD_clinical.Rdata")
diffgene <- allDiff %>% 
  filter(gene !="") %>% 
  filter(adj.P.Val < 0.05) %>% 
  filter(abs(logFC) >2)%>%
  filter(gene_type =="protein_coding" )
load("output/COAD_normalized_counts.Rdata")
normalized_counts <- log(normalized_counts+0.1)
range(normalized_counts)

##将行名中点去掉
term=as.character(row.names(normalized_counts)[1]) 
term =  unlist(strsplit(term, split=".", fixed=T))[1]
rownames(normalized_counts) <- as.character(rownames(normalized_counts))
for (i in 1:nrow(normalized_counts)) {
  print(i)
  term =rownames(normalized_counts)[i]
  term= unlist(strsplit(term, split=".", fixed=T))[1]
  rownames(normalized_counts)[i] =term
}
expro <- cbind(gene_id=rownames(normalized_counts),normalized_counts)
diffgene <- diffgene[,1:2]
expro <- inner_join(diffgene, expro, by = "gene_id")
expro <- expro%>%
  mutate(newcolumn = rowMeans(.[,-c(1:2)])) %>% 
  arrange(desc(newcolumn)) %>% 
  distinct(gene,.keep_all = T) 

rownames(expro) <- expro$gene
expro <- expro[,-(1:2)]

dim(expro)
datExpr=as.data.frame(t(expro));
datExpr <- cbind(sample=rownames(datExpr),datExpr)

for (i in 1:nrow(datExpr)) {
  print(i)
  datExpr$sample[i] <- as.character(substring(datExpr$sample[i],1,16)) 
}
test <- datExpr[1:10,1:10]
datExpr <- datExpr%>%
  mutate(newcolumn = rowMeans(.[,-1])) %>% 
  arrange(desc(newcolumn)) %>% 
  distinct(sample,.keep_all = T) 
colnames(datExpr)[1] <- "TCGA_id"

###合并表达数据和临床数据
clinical <- clinical[,c(1,3,2)]
rt <- inner_join(clinical,datExpr, by = "TCGA_id")
rownames(rt) <- rt[,1]
rt <- rt[,-1]
outTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  outTab=rbind(outTab,
               cbind(id=i,
                     z=coxSummary$coefficients[,"z"],
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
}
outTab = outTab[is.na(outTab$pvalue)==FALSE,]
outTab=outTab[order(as.numeric(as.vector(outTab$pvalue))),]
write.table(outTab,file="output/diffgene_uniCoxResult.txt",sep="\t",row.names=F,quote=F)

sigTab=outTab[as.numeric(as.vector(outTab$pvalue))<0.05,] #P值的筛选阈值可设为0.10
write.table(sigTab,file="output/diffgene_uniCoxResult.Sig.txt",sep="\t",row.names=F,quote=F)

sigGenes=c("futime","fustat")
sigGenes=c(sigGenes,as.vector(sigTab[,1]))
uniSigExp=rt[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
write.table(uniSigExp,file="output/diffgene_uniSigExp.txt",sep="\t",row.names=F,quote=F)
