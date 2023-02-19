#加载需要的R包
rm(list = ls())
library(Seurat)
library(msigdbr)
library(GSVA)
library(tidyverse)
library(clusterProfiler)
library(patchwork)
library(fgsea)
library(dplyr)
library(ggplot2)
#创建一个文件夹用于写分析结果
sam.name <- "GSEA_GSVA"
if(!dir.exists(sam.name)){
  dir.create(sam.name)
}

#加载单个样本单细胞数据分析结果
load(file = "./multi/multi_experiment.rename.Rdata")
#挑选B细胞相对于CD8T细胞特意性高表达的marker基因
B_vs_CD8T=FindMarkers(experiment.rename, ident.1="NKT cell",ident.2="B cell",only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
#制作geneList
#值为log2FC
geneList= B_vs_CD8T$avg_log2FC 
#name为基因名字
names(geneList)= rownames(B_vs_CD8T)
#按FC从大到小排序
geneList=sort(geneList,decreasing = T)
head(geneList)


#从GSEA官网下载GSEA分析需要的基因集
#http://www.gsea-msigdb.org/gsea/index.jsp

#下载免疫相关的基因集，c7: immunologic signature gene sets
gmtfile ='GSEA_data/c7.all.v7.4.symbols.gmt'

#读取gmt文件中的pathway信息
pathway<-read.gmt(gmtfile)
#进行GSEA分析
y <- GSEA(geneList,TERM2GENE =pathway)

#保存GSEA分析结果
write.csv(file="B_vs_CD8T_GSEA_result.csv",data.frame(y))

#气泡图展示显著富集的前20条通路
pdf(file="GSEA_dotplot.pdf",width=15)
dotplot(y,showCategory=20)
dev.off()

#绘制具体通路的GSEA图
library(enrichplot)
pdf(file="CD8_TCELL_VS_BCELL_NAIVE_DN.pdf",width=12)
gseaplot2(y,geneSetID = 'GSE22886_CD8_TCELL_VS_BCELL_NAIVE_DN',pvalue_table=T)
dev.off()

pdf(file="CD8_TCELL_VS_BCELL_NAIVE_UP.pdf",width=12)
gseaplot2(y,geneSetID = 'GSE22886_CD8_TCELL_VS_BCELL_NAIVE_UP',pvalue_table=T) 
dev.off()


##########################################
#GSVA
##########################################
#获取细胞类型
Idents(experiment.rename)
#计算每个基因在每个细胞亚群中的平均表达值
expr <- AverageExpression(experiment.rename, assays = "RNA", slot = "data")[[1]]
View(expr)
#选取非零基因
expr <- expr[rowSums(expr)>0,] 
#转换成矩阵
expr <- as.matrix(expr)


#从GSEA官网下载GSEA分析需要的基因集
#http://www.gsea-msigdb.org/gsea/index.jsp

#下载免疫相关的基因集，h: hallmark gene sets
gmtfile ='GSEA_data/h.all.v7.4.symbols.gmt'

#读取gmt文件中的pathway信息
pathway<-read.gmt(gmtfile)[,c(2,1)]
#去堆叠,转换成list
genesets=unstack(pathway)

#进行GSVA分析
gsva.res <- gsva(expr, genesets, method="ssgsea") 
#保存GSVA分析结果
gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
write.csv(gsva.df, "gsva_res.csv", row.names = F)

#绘制热图
library(pheatmap)
pdf(paste0("./",sam.name,"/GSVA_heatmap.pdf"),width = 8,height = 10)
pheatmap(gsva.res, show_colnames = T, scale = "row")
dev.off()

