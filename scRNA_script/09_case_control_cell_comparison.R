rm(list = ls())
library(Seurat)
library(dplyr)
library(gtable)
library(cowplot)
library(gridExtra)
#创建一个文件夹用于写分析结果
sam.name <- "case_control"
# if(!dir.exists(sam.name)){
#   dir.create(sam.name)
# }

##加载数据
load(file = "./multi/multi_experiment.rename.Rdata")
obj.combined <- experiment.rename

#计算整体的细胞数目和比例
cell.num <- table(Idents(obj.combined))
cell.freq <- round(prop.table(table(Idents(obj.combined)))*100,2)
cell.combined <- rbind(cell.num, cell.freq)
write.csv(cell.combined,paste0("./",sam.name,"/combined_cell_counts_freq.csv"))


#分组计算细胞数目和比例
cell.num.orig.ident <- table(Idents(obj.combined), obj.combined$orig.ident) 
colnames(cell.num.orig.ident) <- paste0(colnames(cell.num.orig.ident),'_cell_counts')
cell.freq.orig.ident <- round(prop.table(table(Idents(obj.combined), obj.combined$orig.ident), margin = 2) *100,2)
colnames(cell.freq.orig.ident) <- paste0(colnames(cell.freq.orig.ident),'_cell_Freq')
cell.orig.ident <- cbind(cell.num.orig.ident, cell.freq.orig.ident)
write.csv(cell.orig.ident,paste0("./",sam.name,"/orig.ident_cell_counts_freq.csv"))

#表格和UMAP展示case和control组的细胞数目变化
pdf(paste0("./",sam.name,"/tsne_multi_samples_split_anno.pdf"),width = 25, height = 8)
p<-DimPlot(obj.combined, reduction = "tsne", split.by = "orig.ident",label=T,repel=T) + NoLegend()
tb <- tableGrob(cell.orig.ident)
plot_grid(p, tb,ncol=2,rel_widths=c(0.6,0.4))
dev.off()


