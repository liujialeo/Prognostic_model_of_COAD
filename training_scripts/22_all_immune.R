rm(list = ls())
library(tidyr)
library(dplyr)
library(GSVA)
#########################################制作marker gene文件
options("repos"=c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
if(!require("data.table")) install.packages("data.table",update = F,ask = F)
if(!require("GSVA")) BiocManager::install("GSVA",update = F,ask = F)

### 直接讲解：ssGSEA
### 需要表达量数据和marker数据
### marker数据的获取
### Pan-cancer Immunogenomic Analyses Reveal Genotype-Immunophenotype Relationships and Predictors of Response to Checkpoint Blockade
### https://www.cell.com/cms/10.1016/j.celrep.2016.12.019/attachment/f353dac9-4bf5-4a52-bb9a-775e74d5e968/mmc3.xlsx

### 1.准备细胞marker
cellMarker <- data.table::fread("resource/cellMarker_allsamples.csv",data.table = F)
colnames(cellMarker)[2] <- "celltype"

cellMarker <- lapply(split(cellMarker,cellMarker$celltype), function(x){
  dd = x$Metagene
  unique(dd)
})

save(cellMarker,file = "resource/cellMarker_allsamples_ssGSEA.Rdata")
#########################################
load(file = "./output/COAD_allDiff.Rdata")
diffgene <- allDiff %>% 
  filter(gene !="") %>% 
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
  datExpr$sample[i] <- as.character(substring(datExpr$sample[i],1,14)) 
}
test <- datExpr[1:10,1:10]
datExpr <- datExpr%>%
  mutate(newcolumn = rowMeans(.[,-1])) %>% 
  arrange(desc(newcolumn)) %>% 
  distinct(sample,.keep_all = T) 
rownames(datExpr) <- datExpr$sample
datExpr <- datExpr[,-1]
expr <- as.matrix(t(datExpr))

### 挺耗时间的，调用了12个线程，17:12开始, 17:37结束
gsva_data <- gsva(expr,cellMarker, method = "ssgsea")
## save(gsva_data,file = "gsva_data_TCGA.Rdata")
## load(file = "gsva_data_TCGA.Rdata")
test <- gsva_data[1:10,1:10]
tcga_gsva <- as.data.frame(t(gsva_data))
test <- tcga_gsva[1:10,1:10]
################## 添加分组信息
##根据基因表达量进行分组
risk_group <- datExpr[,c("UCN","GABRD")]
risk_group$expression_UCN <- ifelse(datExpr$UCN > median(datExpr$UCN),"High","Low")
risk_group$expression_GABRD <- ifelse(datExpr$GABRD > median(datExpr$GABRD),"High","Low")
risk_group <- cbind(TCGA_id=rownames(risk_group),risk_group[,3:4])
tcga_gsva <- cbind(TCGA_id = rownames(tcga_gsva),tcga_gsva)
tcga_gsva <- inner_join(risk_group,tcga_gsva, by = "TCGA_id")

### 调整数据
library(dplyr)
library(tidyr)
dd1 <- tcga_gsva %>% 
  pivot_longer(cols=4:17,
               names_to= "celltype",
               values_to = "NES")

#########根据UCN基因表达量分组
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
display.brewer.all()
###获取自定义颜色
colset <- brewer.pal(10,"Paired") ##brewer.pal选择颜色的最小种类是3
### 箱线图 
##UCN
ggplot(data =dd1, aes(x = celltype, y = NES))+
  geom_boxplot(aes(fill = expression_UCN),outlier.shape = NA)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1, colour = "black"))+
  stat_compare_means(aes(group=expression_UCN), label = "p.signif")+
  scale_fill_manual(values= colset[c(1,2)] )
ggsave(filename = "output/immnue_allsamples_scRNA_UCN.pdf")
##GABRD
ggplot(data =dd1, aes(x = celltype, y = NES))+
  geom_boxplot(aes(fill = expression_GABRD),outlier.shape = NA)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1, colour = "black"))+
  stat_compare_means(aes(group=expression_GABRD), label = "p.signif")+
  scale_fill_manual(values= colset[c(1,2)] )
ggsave(filename = "output/immnue_allsamples_scRNA_GABRD.pdf")


# ### 小提琴
# ggplot(data =dd1, aes(x = celltype, y = NES))+
#   geom_violin(aes(fill = expression_UCN),position = position_dodge(1),scale = "width")+
#   theme_bw()+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1, colour = "black"))+
#   stat_compare_means(aes(group=expression_UCN), label = "p.signif")+
#   scale_fill_manual(values= colset[c(1,2)] )

# ### 混合叠加
# ggplot(data =dd1, aes(x = celltype, y = NES))+
#   geom_boxplot(aes(fill = expression_UCN),position = position_dodge(1),width=.3,outlier.shape = NA)+
#   geom_violin(aes(colour = expression_UCN),position = position_dodge(1),scale = "width",fill=NA)+
#   theme_bw()+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1, colour = "black"))+
#   stat_compare_means(aes(group=expression_UCN), label = "p.signif")+
#   scale_fill_manual(values= colset[c(1,2)] )
