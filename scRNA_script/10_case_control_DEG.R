
library(Seurat)
library(dplyr)
library(cowplot)
library(patchwork)
library(ggplot2)
library(tidyr)

#将默认检测设置为RNA
DefaultAssay(obj.combined) <- "RNA"

#保存每个细胞的细胞类型
obj.combined$celltype <- Idents(obj.combined)

#细胞类型前面加上分组信息
obj.combined$celltype.orig.ident <- paste(Idents(obj.combined), obj.combined$orig.ident, sep = "_")
#将每个细胞的identity转换成celltype.orig.ident
Idents(obj.combined) <- "celltype.orig.ident"

#鉴定Paneth cell type 1细胞在case和control组中差异表达的基因
Paneth_cell_type_1_adenoma_normal.diff <- FindMarkers(obj.combined, ident.1 = "NKT cell_adenoma", 
                                                    ident.2 = "NKT cell_normal", 
                                                    verbose = FALSE)

Paneth_cell_type_1_carcinoma_normal.diff <- FindMarkers(obj.combined, ident.1 = "NKT cell_carcinoma", 
                                                      ident.2 = "NKT cell_normal", 
                                                      verbose = FALSE)
##差异基因取交集
Paneth_cell_type_1_adenoma_normal.diff <- cbind(rownames(Paneth_cell_type_1_adenoma_normal.diff),
                                                Paneth_cell_type_1_adenoma_normal.diff)
Paneth_cell_type_1_carcinoma_normal.diff <- cbind(rownames(Paneth_cell_type_1_carcinoma_normal.diff),
                                                  Paneth_cell_type_1_carcinoma_normal.diff)
colnames(Paneth_cell_type_1_adenoma_normal.diff)[1] <- "gene"
colnames(Paneth_cell_type_1_carcinoma_normal.diff)[1] <- "gene"
Paneth_cell.diff <- inner_join(Paneth_cell_type_1_adenoma_normal.diff,
                               Paneth_cell_type_1_carcinoma_normal.diff,
                               by = "gene")
##筛选差异显著的基因
Paneth_cell.diff.up <- Paneth_cell.diff %>%
  filter(p_val.x < 0.05) %>%
  filter(p_val.y < 0.05) %>%
  filter(avg_log2FC.x > 0) %>%
  filter(avg_log2FC.y > 0)
  
write.csv(Paneth_cell.diff.up,paste0("./",sam.name,"/Paneth_cell.diff.up.csv"))

# 绘制case vs control的散点图
Idents(obj.combined) <- 'celltype'
theme_set(theme_cowplot())
#提取出Paneth_cell_type_1细胞
Paneth_cell_type_1 <- subset(obj.combined, idents = "Paneth cell type 1")
Idents(Paneth_cell_type_1) <- "orig.ident"
avg.Paneth_cell_type_1 <- log1p(AverageExpression(Paneth_cell_type_1, verbose = FALSE)$RNA)
avg.Paneth_cell_type_1 <- data.frame(avg.Paneth_cell_type_1 ,gene=rownames(avg.Paneth_cell_type_1))


genes.to.label = c("NFKBIA","IGHA1","RPL36A","JCHAIN")
p1 <- ggplot(avg.Paneth_cell_type_1, aes(normal, adenoma)) + geom_point() + ggtitle("Paneth cell type 1")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)


pdf(paste0("./",sam.name,"/Diff_genes_scatter_plot.pdf"),width = 25, height = 8)
p1 
dev.off()



# feature plot
pdf(paste0("./",sam.name,"/feature_plot_orig.ident.pdf"),width = 10, height = 4)
FeaturePlot(obj.combined, features = c("NFKBIE"), split.by = "orig.ident", max.cutoff = 3, 
                 cols = c("grey", "red"))
dev.off()

# violin plot
#split.plot=T 拼成一个小提琴图
pdf(paste0("./",sam.name,"/violin_plot_one.pdf"),width = 10, height = 4)
plots <- VlnPlot(obj.combined, features = c("NFKBIE"), split.by = "orig.ident", group.by = "celltype", 
                 pt.size = 0, combine = FALSE,split.plot=T)
wrap_plots(plots = plots, ncol = 1)
dev.off()

#split.plot=F 分开的小提琴图
pdf(paste0("./",sam.name,"/violin_plot_two.pdf"),width = 12, height = 8)
obj.combined$orig.ident <- factor(obj.combined$orig.ident,levels = c("normal","adenoma","carcinoma"))
plots <- VlnPlot(obj.combined, features = c("NFKBIE", "NFKB1","IL6","ICAM1","REL",
                                            "TNFAIP3", "CXCL10","TLR2","TNF","IL1A"), split.by = "orig.ident", group.by = "celltype", 
                 pt.size = 2, combine = FALSE,split.plot=F,stack = T,flip = T, cols = c("#650d6b","#77d946","#efe400"))
wrap_plots(plots = plots, ncol = 1)
dev.off()

load(file = "./multi/multi_experiment.aggregate.Rdata")
DimPlot(object = experiment.aggregate, 
        split.by ="orig.ident", 
        pt.size=0.5,reduction = "tsne",
        label = T)


