### 本节任务：KEGG分析
rm(list = ls())
library(clusterProfiler)
load(file = "output/COAD_allDiff.Rdata")
### 筛选差异基因
library(dplyr)
diffgene <- allDiff %>% 
  filter(gene !="") %>% 
  filter(adj.P.Val < 0.05) %>% 
  filter(abs(logFC) >1)

##############################################################
#**KEGG分析**
##############################################################
library(KEGG.db)
EGG <- enrichKEGG(gene = diffgene$entrez,
                  organism = 'hsa',
                  pvalueCutoff = 0.05,
                  use_internal_data =T)
KEGG_df <- as.data.frame(EGG)
## 画图
barplot(EGG)
dotplot(EGG)

######################################################################
### KEGG的富集分析比较特殊，他的背后是个网站

browseKEGG(EGG, 'hsa04151')


library(clusterProfiler)
diffgene <- allDiff %>% 
  filter(gene !="") %>% 
  filter(adj.P.Val < 0.05) %>% 
  filter(abs(logFC) >1)

## geneList 三部曲
## 1.获取基因logFC
geneList <- diffgene$logFC
## 2.命名
names(geneList) = diffgene$entrez
## 3.排序很重要
geneList = sort(geneList, decreasing = TRUE)

head(geneList)
############################################################################
### pathview可视化
library(pathview)
pathway.id = "hsa04151"
pv.out <- pathview(gene.data  = geneList,
                   pathway.id = pathway.id,
                   species    = "hsa")

#########################################################################
### pathview批量画图
### 新建文件夹，
dir.create("output/pathview_out")
### 设置工作目录到想要的地方

getwd()
setwd("D:/00_doctor/03_data/06_TCGA_WGCNA_molecular_cox/02_COAD/output/pathview_out")
### 然后循环绘图
for (pathway.id in KEGG_df$ID ){
  pathview(gene.data  = geneList,
           pathway.id = pathway.id,
           species    = "hsa"
  )
}

### 结束后切换回原来的工作目录
### 确认
setwd("D:/00_doctor/03_data/06_TCGA_WGCNA_molecular_cox/02_COAD")

getwd()


