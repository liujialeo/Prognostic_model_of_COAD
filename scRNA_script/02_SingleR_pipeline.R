library(SingleR)
# BiocManager::install("celldex")
library(celldex)

refdata <- HumanPrimaryCellAtlasData()
testdata <- GetAssayData(experiment.aggregate, slot="data")
clusters <- experiment.aggregate@meta.data$seurat_clusters
cellpred <- SingleR(test = testdata, ref = refdata, 
                    labels = refdata$label.main, 
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")
celltype <- data.frame(ClusterID=rownames(cellpred), 
                      celltype=cellpred$labels, stringsAsFactors = F)

experiment.aggregate@meta.data$celltype = "Unknown"
for(i in 1:nrow(celltype)){
  experiment.aggregate@meta.data[which(experiment.aggregate@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]
}

plot_anno <- DimPlot(experiment.aggregate, group.by="celltype", label=F, label.size=5, reduction='tsne') #+ 
  #NoLegend()
ggsave(filename = paste0("./",sam.name,"/CellCluster-TSNEPlot_res0.5_SingleR_",max(dim.use),"PC.pdf"),
       plot = plot_anno) 
