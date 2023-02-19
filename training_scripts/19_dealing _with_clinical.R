rm(list = ls())
library(tidyr)
library(dplyr)

load(file = "output/COAD_clinical.Rdata")
load(file = "./output/COAD_allDiff.Rdata")
diffgene <- allDiff %>% 
  filter(gene !="") %>% 
  filter(adj.P.Val < 0.01) %>% 
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
rt_plot <- inner_join(clinical, datExpr, by = "TCGA_id")
rt_plot <- rt_plot[,-c(4:9)]

## 修改列名
colnames(rt_plot)[1:3] <- c("TCGA_id", "fustat", "futime")
rownames(rt_plot) <- rt_plot[,1]
rt_plot <- rt_plot[,-1]
rt_plot <- rt_plot[,c(2,1,3:2060)]

## 保存数据
save(rt_plot,file = "output/rt_plot.Rdata")

