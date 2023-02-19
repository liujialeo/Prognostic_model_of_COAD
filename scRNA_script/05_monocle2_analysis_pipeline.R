# 该代码用于进行拟时序分析

#加载分析使用的包
rm(list = ls())
library(Seurat)
library(monocle)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
load(file = "./multi/multi_experiment.rename.Rdata")
#创建一个文件夹用于写分析结果
result.name <- "monocle2"
if(!dir.exists(result.name)){
  dir.create(result.name)
}
experiment.rename$celltype <- Idents(experiment.rename)
experiment.aggregate <- experiment.rename
##使用monocle2进行拟时序分析
#构造表达及注释数据，提取CCA之后的数据
exp.matrix<-as(as.matrix(experiment.aggregate@assays$RNA@data), 'sparseMatrix')
feature_ann<-data.frame(gene_id=rownames(exp.matrix),gene_short_name=rownames(exp.matrix))
rownames(feature_ann)<-rownames(exp.matrix)
exp_fd<-new("AnnotatedDataFrame", data = feature_ann)
sample_ann<-experiment.aggregate@meta.data
rownames(sample_ann)<-colnames(exp.matrix)
exp_pd<-new("AnnotatedDataFrame", data =sample_ann)

#生成monocle对象
exp.monocle<-newCellDataSet(exp.matrix,phenoData =exp_pd,featureData =exp_fd,expressionFamily=negbinomial.size())
head(pData(exp.monocle))
head(fData(exp.monocle))

#计算sizefactor
exp.monocle <- estimateSizeFactors(exp.monocle)
exp.monocle <- estimateDispersions(exp.monocle)

#根据seurat cluster计算差异表达基因并挑选用于构建拟时序轨迹的基因
diff_test_res<-differentialGeneTest(exp.monocle,fullModelFormulaStr = "~seurat_clusters") 
ordering_genes<-row.names (subset(diff_test_res, qval < 0.01))
exp.monocle<-setOrderingFilter(exp.monocle, ordering_genes)
plot_ordering_genes(exp.monocle)

#DDRTree的方法降维并构建拟时序
exp.monocle<-reduceDimension(exp.monocle, max_components = 2, reduction_method = "DDRTree")
exp.monocle<-orderCells(exp.monocle)
colnames(pData(exp.monocle))

#修改monocle对象中的列名示例
names(pData(exp.monocle))[names(pData(exp.monocle))=="seurat_clusters"]="Cluster"

#########################绘图
#以state来展示
pdf(paste0("./",result.name,"/trajectory_state.pdf"),width = 16,height = 10)
plot_cell_trajectory(exp.monocle)
dev.off()

#以细胞类型来展示
pdf(paste0("./",result.name,"/trajectory_celltype.pdf"),width = 16,height = 10)
plot_cell_trajectory(exp.monocle, color_by = 'celltype') 
dev.off()

#分别展示每一个细胞类型的发育轨迹
pdf(paste0("./",result.name,"/trajectory_celltype_separate.pdf"),width = 16,height = 10)
plot_cell_trajectory(exp.monocle, color_by = 'celltype')+ facet_wrap(~celltype, nrow = 3) + NoLegend()
dev.off()

#查看特定的gene在发育轨迹上的表达情况
pdf(paste0("./",result.name,"/5gene_trajectory.pdf"),width = 12,height = 6)
plot_cell_trajectory(exp.monocle,markers=c("NFKBIE","NFKB1","LGR5"),use_color_gradient=T)
dev.off()

#heatmap
top5=all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
sig_gene_names <- unique(top5$gene)
pseudotemporalplot<- plot_pseudotime_heatmap(exp.monocle[sig_gene_names,],
                                             num_clusters = 19,  #亚群数需要对应修改
                                             cores = 4,
                                             hmcols = NULL,
                                             show_rownames = T,
                                             return_heatmap = T)
pdf(paste0("./",result.name,"/monocle_heatmap.pdf"),width = 16,height = 10)
pseudotemporalplot
dev.off()



#将不同分组情况的拟时序轨迹图画到一起
plot1<-plot_cell_trajectory(exp.monocle, color_by = "celltype",cell_size=1)
# plot2<-plot_cell_trajectory(exp.monocle, color_by = "sample",cell_size=1)
# plot3<-plot_cell_trajectory(exp.monocle, color_by = "batch",cell_size=1)
plot4<-plot_cell_trajectory(exp.monocle, color_by = "State",cell_size=1)
plot5<-plot_cell_trajectory(exp.monocle, color_by = "Pseudotime",cell_size=1)
plot6 <- plot_cell_trajectory(exp.monocle,markers=c("NFKBIE"),use_color_gradient=T)
pdf(paste0("./",result.name,"/trajectory_plot.pdf"),width = 16,height = 10)
CombinePlots(plots = list(plot1, plot4,plot5,plot6),legend = NULL) 
dev.off()
rm(plot1,plot2,plot3,plot4,plot5,plot6)




###保存拟时序计算结果
save(exp.monocle, file=paste0("./",result.name,"/",result.name,"_exp.monocle.RData"))


