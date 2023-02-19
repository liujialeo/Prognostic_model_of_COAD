####  Code Description              ####
#---  1. Written by WoLin @ 2019.03.15，last update 19.08.26 ---#
#---  2. private use for scRNA    ---#

# 修改具体类别的名???(需要基于marker基因来定???)

library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(ggsci)
sam.name <- "multi"
dim.use <- 1:20

load(file = "./multi/multi_experiment.aggregate.Rdata")

experiment.rename <- RenameIdents(
  object = experiment.aggregate,
  "0" = "NKT cell",
  "1" = "B cell",
  "2" = "CD4+ cytotoxic T cell",
  "3" = "Paneth cell",
  "4" = "Goblet cell",
  "5" = "CD8+ T cell ",
  "6" = "Th17 cell",
  "7" = "LGR5+ stem cell",
  "8" = "Enterocyte progenitor cell",
  "9" = "Paneth cell",
  "10" = "Goblet cell",
  "11" = "LGR5+ stem cell",
  "12" = "LGR5+ stem cell",
  "13" = "Enterocyte",
  "14" = "NKT cell",
  "15" = "Enterocyte",
  "16" = "Plasma cell",
  "17" = "LGR5+ stem cell",
  "18" = "Enterocyte",
  "19" = "Enterocyte",
  "20" = "Goblet cell",
  "21" = "B cell",
  "22" = "NK cell",
  "23" = "DCLK1+ progenitor cell",
  "24" = "MKI67+ progenitor cell",
  "25" = "B cell"
)
experiment.rename$celltype <- Idents(experiment.rename)
experiment.rename$orig.ident <- factor(experiment.rename$orig.ident,levels = c("normal","adenoma","carcinoma"))

save(experiment.rename, file = "./multi/multi_experiment.rename.Rdata")

#命名后的heatmap
pdf(paste0("./",sam.name,"/markers_heatmap_celltype_rename.pdf"), width = 14, height = 16)
DoHeatmap(experiment.rename, features = top5$gene, size = 3, angle = -50,hjust=0.8) + NoLegend()
dev.off()

pdf(paste0("./",sam.name,"/CellCluster-TSNEPlot_Rename_",max(dim.use),"PC.pdf"),width = 8,height = 5)
DimPlot(object = experiment.rename, pt.size=0.5,
        reduction = "tsne",label = F) +
  ggsci::scale_color_igv()
dev.off()

pdf(paste0("./",sam.name,"/CellCluster-TSNEPlot_Rename_split",max(dim.use),"PC.pdf"),width = 15,height = 5)
DimPlot(object = experiment.rename, pt.size=0.5,split.by ="orig.ident", 
        reduction = "tsne",label =F) +
  ggsci::scale_color_igv()
dev.off()

# 计算特定两组细胞之间的差异基???
sub.markers <- FindMarkers(experiment.rename,
                           ident.1 = "NKT cell",
                           ident.2 = "Plasma cell")

# 只展示特定类群细胞中top-marker基因热图
pdf(paste0("./",sam.name,"/MarkerGene-Heatmap_all_cluster_tsne_",max(dim.use),"PC.pdf"),width = 8,height = 5)
DoHeatmap(experiment.rename, features = top5$gene,size = 2,
          cells = rownames(subset(experiment.rename@meta.data, seurat_clusters %in% c("Plasma B cell","Memory B cell")))) +
  theme(legend.position = "none")
dev.off()



# 计算每一类样本中不同细胞的比例并画图
library(RColorBrewer)
meta_data <- experiment.rename@meta.data 
plot_data <- data.frame(table(meta_data$orig.ident,meta_data$celltype))
plot_data$Total <- apply(plot_data,1,function(x)sum(plot_data[plot_data$Var1 == x[1],3]))
plot_data <- plot_data %>% mutate(Percentage = round(Freq/Total,3) * 100)
ggplot(plot_data,aes(x = Var1,y = Percentage,fill = Var2)) +
  geom_bar(stat = "identity",position = "stack") +
  theme_classic() + 
  scale_fill_manual(values = colorRampPalette(brewer.pal(11, "Spectral"))(14)) + # 设置填充颜色
  theme(axis.title.x = element_blank()) + labs(fill = "Cluster")
ggsave(paste0("./",sam.name,"/CellCluster-TSNEPlot_Rename_stackplot_",max(dim.use),"PC.pdf"),width = 8,height = 5)

# 提取子矩阵重新聚类分???
cells_sub <- subset(experiment.aggregate@meta.data, 
                    seurat_clusters %in% c(4,6))
scRNA_sub <- subset(experiment.aggregate, 
                   cells=row.names(cells_sub))