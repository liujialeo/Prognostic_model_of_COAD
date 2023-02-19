####  Code Description              ####
#---  1. Written by WoLin @ 2019.03.15，last update 19.07.14 ---#
#---  2. Analysis for single sample    ---#
#---  3. Support 10X data & expression matrix ---#
#---  4. Need to change:sample data, sam.name, dims ---#
rm(list = ls())
#### 1. 加载分析使用的工具包 ####
library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(ggsci)

#### 2. 读入原始表达数据 ####
#以下两种方式二选一
#10X 数据
getwd()
setwd("/home/data/t170309/public_data/01_colon_cancer_scRNA/02_GSE161277_normal_adenoma_carcinoma")
normal_1 <- Read10X("./data/normal_1/")
normal_2 <- Read10X("./data/normal_2/")
normal_3 <- Read10X("./data/normal_3/")

adenoma_1 <- Read10X("./data/adenoma_1/")
adenoma_2 <- Read10X("./data/adenoma_2/")
adenoma_3 <- Read10X("./data/adenoma_3/")
adenoma_4 <- Read10X("./data/adenoma_4/")

carcinoma_1 <- Read10X("./data/carcinoma_1/")
carcinoma_2 <- Read10X("./data/carcinoma_2/")
carcinoma_3 <- Read10X("./data/carcinoma_3/")
carcinoma_4 <- Read10X("./data/carcinoma_4/")

# 这里的列名就是barcode
colnames(normal_1) <- paste(colnames(normal_1),"normal_1",sep = "_")
colnames(normal_2) <- paste(colnames(normal_2),"normal_2",sep = "_")
colnames(normal_3) <- paste(colnames(normal_3),"normal_3",sep = "_")

colnames(adenoma_1) <- paste(colnames(adenoma_1),"adenoma_1",sep = "_")
colnames(adenoma_2) <- paste(colnames(adenoma_2),"adenoma_2",sep = "_")
colnames(adenoma_3) <- paste(colnames(adenoma_3),"adenoma_3",sep = "_")
colnames(adenoma_4) <- paste(colnames(adenoma_4),"adenoma_4",sep = "_")

colnames(carcinoma_1) <- paste(colnames(carcinoma_1),"carcinoma_1",sep = "_")
colnames(carcinoma_2) <- paste(colnames(carcinoma_2),"carcinoma_2",sep = "_")
colnames(carcinoma_3) <- paste(colnames(carcinoma_3),"carcinoma_3",sep = "_")
colnames(carcinoma_4) <- paste(colnames(carcinoma_4),"carcinoma_4",sep = "_")

#创建一个文件夹用于写分析结果
sam.name <- "multi"
if(!dir.exists(sam.name)){
  dir.create(sam.name)
}

#### 3. 创建Seurat分析对象 ####
normal_1 <- CreateSeuratObject(
  normal_1,
  project = "multi", 
  min.cells = 10,
  min.features = 200,
  names.field = 2,
  names.delim = "_")

normal_2 <- CreateSeuratObject(
  normal_2,
  project = "multi", 
  min.cells = 10,
  min.features = 200,
  names.field = 2,
  names.delim = "_")

normal_3 <- CreateSeuratObject(
  normal_3,
  project = "multi", 
  min.cells = 10,
  min.features = 200,
  names.field = 2,
  names.delim = "_")
##合并
normal <- merge(normal_1, y = c(normal_2,normal_3), 
                add.cell.ids = c("normal_1","normal_2", "normal_3"), 
                project = "normal")

adenoma_1 <- CreateSeuratObject(
  adenoma_1,
  project = "multi", 
  min.cells = 10,
  min.features = 200,
  names.field = 2,
  names.delim = "_")

adenoma_2 <- CreateSeuratObject(
  adenoma_2,
  project = "multi", 
  min.cells = 10,
  min.features = 200,
  names.field = 2,
  names.delim = "_")

adenoma_3 <- CreateSeuratObject(
  adenoma_3,
  project = "multi", 
  min.cells = 10,
  min.features = 200,
  names.field = 2,
  names.delim = "_")
adenoma_4 <- CreateSeuratObject(
  adenoma_4,
  project = "multi", 
  min.cells = 10,
  min.features = 200,
  names.field = 2,
  names.delim = "_")
##合并
adenoma <- merge(adenoma_1, y = c(adenoma_2,adenoma_3,adenoma_4), 
                add.cell.ids = c("adenoma_1","adenoma_2", "adenoma_3","adenoma_4"), 
                project = "adenoma")

carcinoma_1 <- CreateSeuratObject(
  carcinoma_1,
  project = "multi", 
  min.cells = 10,
  min.features = 200,
  names.field = 2,
  names.delim = "_")

carcinoma_2 <- CreateSeuratObject(
  carcinoma_2,
  project = "multi", 
  min.cells = 10,
  min.features = 200,
  names.field = 2,
  names.delim = "_")

carcinoma_3 <- CreateSeuratObject(
  carcinoma_3,
  project = "multi", 
  min.cells = 10,
  min.features = 200,
  names.field = 2,
  names.delim = "_")

carcinoma_4 <- CreateSeuratObject(
  carcinoma_4,
  project = "multi", 
  min.cells = 10,
  min.features = 200,
  names.field = 2,
  names.delim = "_")
##合并
carcinoma <- merge(carcinoma_1, y = c(carcinoma_2,carcinoma_3,carcinoma_4), 
                 add.cell.ids = c("carcinoma_1","carcinoma_2", "carcinoma_3","carcinoma_4"), 
                 project = "carcinoma")

View(normal[1:20,1:20])
View(adenoma[1:20,1:20])
View(carcinoma[1:20,1:20])

experiment.aggregate <- merge(normal, y = c(adenoma,carcinoma), 
                              add.cell.ids = c("normal","adenoma", "carcinoma"), 
                              project = "experiment.aggregate")

#将数据写到文件中一边后续分析使用
save(experiment.aggregate,file=paste0("./",sam.name,"/",sam.name,"_raw_SeuratObject.RData"))

#### 4. 数据概览 & QC ####
#查看SeuratObject中的对象
slotNames(experiment.aggregate)
#assay
experiment.aggregate@assays
#细胞及细胞中基因与RNA数量
dim(experiment.aggregate@meta.data)
View(experiment.aggregate@meta.data)

#转换成表达矩阵
experiment.aggregate.matrix <- as.matrix(experiment.aggregate@assays$RNA@counts)

##QC：统计线粒体基因在每个细胞中的占比
experiment.aggregate[["percent.mt"]] <- PercentageFeatureSet(experiment.aggregate, 
                                                             pattern = "^MT-")
pdf(paste0("./",sam.name,"/QC-VlnPlot.pdf"),width = 8,height = 4.5)
VlnPlot(experiment.aggregate, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3)
dev.off()

##QC：统计基因数，RNA，线粒体基因分布
gene.freq <- do.call("cbind", tapply(experiment.aggregate@meta.data$nFeature_RNA,experiment.aggregate@meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
rna.freq <- do.call("cbind", tapply(experiment.aggregate@meta.data$nCount_RNA,experiment.aggregate@meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
mt.freq <- do.call("cbind", tapply(experiment.aggregate@meta.data$percent.mt,experiment.aggregate@meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
freq.combine <- as.data.frame(cbind(gene.freq,rna.freq,mt.freq))
colnames(freq.combine) <- c(paste(colnames(gene.freq),"Gene",sep = "_"),
                            paste(colnames(rna.freq),"RNA",sep = "_"),
                            paste(colnames(mt.freq),"MT",sep = "_"))
write.table(freq.combine,file = paste0(sam.name,"/QC-gene_frequency.txt"),quote = F,sep = "\t")
rm(gene.freq,rna.freq,mt.freq)
View(freq.combine)

##QC：基因数与线粒体基因以及RNA数量的分布相关性
plot1 <- FeatureScatter(experiment.aggregate, pt.size = 0.1,feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(experiment.aggregate, pt.size = 0.1,feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pdf(paste0("./",sam.name,"/QC-FeatureScatter.pdf"),width = 8,height = 4.5)
CombinePlots(plots = list(plot1, plot2),legend = "none")
dev.off()
rm(plot1,plot2)

#### 5. 筛选细胞 ####
cat("Before filter :",nrow(experiment.aggregate@meta.data),"cells\n")
experiment.aggregate <- subset(experiment.aggregate, 
                               subset = 
                                 nFeature_RNA > 500 & 
                                 nFeature_RNA < 5000 & 
                                 percent.mt < 50)
cat("After filter :",nrow(experiment.aggregate@meta.data),"cells\n")

#### 6. 表达量标准化 ####
experiment.aggregate <- NormalizeData(experiment.aggregate, 
                                      normalization.method = "LogNormalize",
                                      scale.factor = 10000)

#计算表达量变化显著的基因FindVariableFeatures
experiment.aggregate <- FindVariableFeatures(experiment.aggregate, 
                                             selection.method = "vst",
                                             nfeatures = 1000)

#展示标准化之后的整体表达水平
top10 <- head(x = VariableFeatures(experiment.aggregate), 10)
plot1 <- VariableFeaturePlot(experiment.aggregate)
plot2 <- LabelPoints(plot = plot1, points = top10)
pdf(file = paste0(sam.name,"/Norm-feature_variable_plot.pdf"),width = 8,height = 5)
CombinePlots(plots = list(plot1, plot2),legend = "none")
dev.off()

#### 7. 均一化与PCA ####
#均一化（需要一点时间）
experiment.aggregate <- ScaleData(

    object = experiment.aggregate,
  do.scale = FALSE,
  do.center = FALSE,
  vars.to.regress = c("orig.ident","percent.mt"))

#PCA计算
experiment.aggregate <- RunPCA(object = experiment.aggregate, 
                               features = VariableFeatures(experiment.aggregate),
                               verbose = F,npcs = 50)

#PCA结果展示-1
pdf(paste0("./",sam.name,"/PCA-VizDimLoadings.pdf"),width = 7,height = 5)
VizDimLoadings(experiment.aggregate, dims = 1:2, reduction = "pca")
dev.off()

#PCA结果展示-2
pdf(paste0("./",sam.name,"/PCA-DimPlot.pdf"),width = 5,height = 4)
DimPlot(experiment.aggregate, reduction = "pca")
dev.off()

#PCA结果展示-3
pdf(paste0("./",sam.name,"/PCA-DimHeatmap.pdf"),width = 5,height = 4)
DimHeatmap(experiment.aggregate, dims = 1:6, cells = 500, balanced = TRUE)
dev.off()

#### 8. 确定细胞类群分析PC ####
#耗时较久
experiment.aggregate <- JackStraw(experiment.aggregate, num.replicate = 100,dims = 40)
experiment.aggregate <- ScoreJackStraw(experiment.aggregate, dims = 1:40)
pdf(paste0("./",sam.name,"/PCA-JackStrawPlot_40.pdf"),width = 6,height = 5)
JackStrawPlot(object = experiment.aggregate, dims = 1:40)
dev.off()

#碎石图
pdf(paste0("./",sam.name,"/PCA-ElbowPlot.pdf"),width = 6,height = 5)
ElbowPlot(experiment.aggregate,ndims = 40)
dev.off()

#确定用于细胞分群的PC
dim.use <- 1:20

#### 9. 细胞分群TSNE算法 ####
#TSNE算法
experiment.aggregate <- FindNeighbors(experiment.aggregate, dims = dim.use)
experiment.aggregate <- FindClusters(experiment.aggregate, resolution = 0.5)

experiment.aggregate <- RunTSNE(experiment.aggregate, dims = dim.use,
                                do.fast = TRUE)
# experiment.aggregate <- RunUMAP(experiment.aggregate, dims = dim.use, 
#                                 do.fast = TRUE)
pdf(paste0("./",sam.name,"/CellCluster-TSNEPlot_res0.5_",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = experiment.aggregate, pt.size=0.5,label = T,reduction = "tsne")
dev.off()

#按照数据来源分组展示细胞异同--画在一张图中
experiment.aggregate$orig.ident <- factor(experiment.aggregate$orig.ident,levels = c("normal","adenoma","carcinoma")) #按照样本名称排顺序

pdf(paste0("./",sam.name,"/CellCluster-TSNEPlot_SamGroup_",max(dim.use),"PC.pdf"),width = 6,height = 4.6)
DimPlot(object = experiment.aggregate, 
        group.by="orig.ident", 
        pt.size=0.5,reduction = "tsne")+
  ggsci::scale_color_jama()
dev.off()

#按照数据来源分组展示细胞异同--画在多张图中
experiment.aggregate$orig.ident <- factor(experiment.aggregate$orig.ident,levels = c("normal","adenoma","carcinoma"))
pdf(paste0("./",sam.name,"/CellCluster-TSNEPlot_SamGroup_slipt_",max(dim.use),"PC.pdf"),width = 8,height = 4)
DimPlot(object = experiment.aggregate, 
        split.by ="orig.ident", 
        pt.size=0.5,reduction = "tsne")
dev.off()

table(experiment.aggregate@meta.data$orig.ident)

#### 10. 计算marker基因 ####
#这一步计算的时候可以把min.pct以及logfc.threshold调的比较低，然后再基于结果手动筛选
all.markers <- FindAllMarkers(experiment.aggregate, only.pos = TRUE, 
                              min.pct = 0.3, logfc.threshold = 0.25)
write.table(all.markers,
            file=paste0("./",sam.name,"/",sam.name,"_total_marker_genes_tsne_",max(dim.use),"PC.txt"),
            sep="\t",quote = F,row.names = F)

# 遍历每一个cluster然后展示其中前4个基因
marker.sig <- all.markers %>% 
  mutate(Ratio = round(pct.1/pct.2,3)) %>%
  filter(p_val_adj <= 0.05)  # 本条件为过滤统计学不显著的基因

for(cluster_id in unique(marker.sig$cluster)){
  # cluster.markers <- FindMarkers(experiment.aggregate, ident.1 = cluster, min.pct = 0.3)
  # cluster.markers <- as.data.frame(cluster.markers) %>% 
  #   mutate(Gene = rownames(cluster.markers))
  cl4.genes <- marker.sig %>% 
    filter(cluster == cluster_id) %>%
    arrange(desc(avg_log2FC))
  cl4.genes <- cl4.genes[1:min(nrow(cl4.genes),4),"gene"]
  
  #VlnPlot
  pvn <- VlnPlot(experiment.aggregate, features = cl4.genes,ncol = 2)
  pdf(paste0("./",sam.name,"/MarkerGene-VlnPlot_cluster",cluster_id,"_tsne_",max(dim.use),"PC.pdf"),width = 7,height = 6)
  print(pvn)
  dev.off()
  
  #Feather plot 
  pvn <- FeaturePlot(experiment.aggregate,features=cl4.genes,ncol = 2)
  pdf(paste0("./",sam.name,"/MarkerGene-FeaturePlot_cluster",cluster_id,"_tsne_",max(dim.use),"PC.pdf"),width = 7,height = 6)
  print(pvn)
  dev.off()
  
  #RidgePlot
  pvn<-RidgePlot(experiment.aggregate, features = cl4.genes, ncol = 2)
  pdf(paste0("./",sam.name,"/MarkerGene-RidgePlot_cluster",cluster_id,"_tsne_",max(dim.use),"PC.pdf"),width = 7,height = 6)
  print(pvn)
  dev.off()
}


#热图展示Top marker基因
#筛选top5的marker基因，可以通过参数改为其他数值
top5 <- marker.sig %>% group_by(cluster) %>% 
  top_n(n = 5, wt = avg_log2FC)

#top-marker基因热图
pdf(paste0("./",sam.name,"/MarkerGene-Heatmap_all_cluster_tsne_",max(dim.use),"PC.pdf"))
DoHeatmap(experiment.aggregate, features = top5$gene,size = 2) +
  theme(legend.position = "none", 
        axis.text.y = element_text(size = 6))
dev.off()

#top-marker基因dotplot
pdf(paste0("./",sam.name,"/MarkerGene-DotPlot_all_cluster_tsne_",max(dim.use),"PC.pdf"),width = 50,height = 5)
DotPlot(experiment.aggregate, features = unique(top5$gene))+
  RotatedAxis()
dev.off()



#筛选top5的marker基因，可以通过参数改为其他数值
top10 <- marker.sig %>% group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)
#整理成表格，只显示基因名字
top10_table=unstack(top10, gene ~ cluster)
names(top10_table)=gsub("X","cluster",names(top10_table))
write.csv(top10,paste0("./",sam.name,"/top10_marker_genes.csv"),row.names=F)
write.csv(top10_table,paste0("./",sam.name,"/top10_table_marker_genes.csv"),row.names=F)
##保存数据
save(experiment.aggregate,file = "./multi/multi_experiment.aggregate.Rdata")

