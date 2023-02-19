
rm(list = ls())
library(Seurat)
library(tidyverse)
#创建一个文件夹用于写分析结果
sam.name <- "cell_cycle"
if(!dir.exists(result.name)){
  dir.create(result.name)
}
load(file = "./multi/multi_experiment.rename.Rdata")
#创建seurat矩阵
test.seu=experiment.rename

###############################################
#为了节省时间这里抽取5000个基因进行演示
#真实情况下请删除这一部分，用完整数据进行分析
###############################################
set.seed(123)
#cc.genes为细胞周期相关的基因，为list，有两个元素s.genes，g2m.genes
cc_gene=unlist(cc.genes)
#从数据中抽取5000个基因，加上97个细胞周期相关的基因，去除重复基因
subgene=unique(c(sample(rownames(test.seu),5000),as.character(cc_gene)))
#抽取相应的seurat对象
test.seu=test.seu[subgene,]
################################################

#下面内容跟分析单个单细胞样本一样
test.seu <- NormalizeData(test.seu)
test.seu <- FindVariableFeatures(test.seu, selection.method = "vst")
test.seu <- ScaleData(test.seu, features = rownames(test.seu))

test.seu <- RunPCA(test.seu, npcs = 50, verbose = FALSE)
test.seu <- FindNeighbors(test.seu, dims = 1:30)
test.seu <- FindClusters(test.seu, resolution = 0.4)

test.seu <- RunUMAP(test.seu, dims = 1:30)
test.seu <- RunTSNE(test.seu, dims = 1:30)

#绘制tsne图
pdf(paste0("./",sam.name,"/tsne.pdf"),width = 7,height = 8)
p0=DimPlot(test.seu,reduction = "tsne",group.by = "seurat_clusters")
p0
dev.off()

#featureplot查看s期基因PCNA和G2/M期基因MKI67表达情况
#s期基因PCNA，G2/M期基因MKI67均高表达，说明细胞均处于细胞周期
pdf(paste0("./",sam.name,"/PCNA_MKI67_expression.pdf"),width = 12,height = 8)
FeaturePlot(test.seu,features = c("PCNA","MKI67"),reduction = "tsne")
dev.off()

#计算细胞周期分值，判断这些细胞分别处于什么期
s.genes=cc.genes$s.genes
g2m.genes=cc.genes$g2m.genes
test.seu <- CellCycleScoring(test.seu, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

#这一步之后test.seu@meta.data数据框会多4列
#S.Score、G2M.Score、Phase、old.ident
head(test.seu@meta.data,2)

#评分和细胞周期phase的关系
#S.Score较高的为S期，G2M.Score较高的为G2M期，都比较低的为G1期
pdf(paste0("./",sam.name,"/cellcycle_score.pdf"),width = 7,height = 8)
test.seu@meta.data  %>% ggplot(aes(S.Score,G2M.Score))+geom_point(aes(color=Phase))+theme_minimal()
dev.off()

#校正之前，细胞分群与细胞周期的关系
p1=DimPlot(test.seu,reduction = "tsne")
pdf(paste0("./",sam.name,"/before_correction.pdf"),width = 7,height = 8)
p1
dev.off()

#去除细胞周期影响
test.seu2 <- ScaleData(test.seu, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(test.seu))

#重新跑PCA，分群聚类
test.seu2 <- RunPCA(test.seu2, npcs = 50, verbose = FALSE)
test.seu2 <- FindNeighbors(test.seu2, dims = 1:30)
test.seu2 <- FindClusters(test.seu2, resolution = 0.4)

test.seu2 <- RunUMAP(test.seu2, dims = 1:30)
test.seu2 <- RunTSNE(test.seu2, dims = 1:30)
#校正之后，细胞分群与细胞周期的关系

p2=DimPlot(test.seu2,reduction = "tsne")
p3=DimPlot(test.seu2,reduction = "tsne",group.by = "Phase")
pdf(paste0("./",sam.name,"/before_after_cellcycle_correction.pdf"),width = 7,height = 8)
p0+p1+p2+p3
dev.off()

