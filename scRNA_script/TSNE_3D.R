
library(RColorBrewer)

set.seed(1)
data.combined <- experiment.aggregate
data.combined <- RunTSNE(data.combined, reduction = "pca", dims = 1:20,dim.embed=3)
dim(data.combined@reductions$tsne)
data.combined <- RenameIdents(
  object = data.combined,
  "0" = "Fibroblasts 1",
  "1" = "Endothelial cells 1",
  "2" = "stem cell 1",
  "3"= "Endothelial cells 2",
  "4" ="Fibroblasts 2",
  "5" = "stem cell 2",
  "6"= "B cells",
  "7" = "Fibroblasts 3",
  "8" = "Neurons 1",
  "9"= "stem cell 3",
  "10" = "Endothelial cells 3",
  "11" = "Neurons 2"
)

tmp.tsne.3<-Embeddings(object = data.combined[["tsne"]])
cb_palette <- c("#ed1299", "#09f9f5", "#246b93", "#cc8e12", "#d561dd", "#c93f00", "#ddd53e","#4aef7b", 
                "#e86502", "#9ed84e", "#39ba30", "#6ad157", "#8249aa", "#99db27", "#e07233", "#ff523f",
                "#ce2523", "#f7aa5d", "#cebb10", "#03827f", "#931635", "#373bbf", "#a1ce4c", "#ef3bb6", 
                "#d66551","#1a918f", "#ff66fc", "#2927c4", "#7149af" ,"#57e559" ,"#8e3af4" ,"#f9a270" ,
                "#22547f", "#db5e92","#edd05e", "#6f25e8", "#0dbc21", "#280f7a", "#6373ed", "#5b910f" ,
                "#7b34c1" ,"#0cf29a","#d80fc1","#dd27ce", "#07a301", "#167275", "#391c82", "#2baeb5",
                "#925bea", "#63ff4f")

cb_palette.use <- cb_palette[1:length(unique(data.combined$seurat_clusters))]
col_match <- data.frame(cluster=unique(data.combined$seurat_clusters),col=cb_palette.use)
col_draw<- col_match[match(data.combined$seurat_clusters,col_match[,1]),2]




library(rgl)
p <- plot3d(
  tmp.tsne.3,
  col = col_draw,
  type = 'p', radius = .001,axes=T,box=F,label = T)

ggsave(p,filename = "./multi/CellCluster_TSNEPlot_3D.focus")

library(plotly)
tmp.tsne.3 <- as.data.frame(tmp.tsne.3)
fig <- plot_ly(tmp.tsne.3, x = ~tSNE_1, y = ~tSNE_2, z = ~tSNE_3, color =data.combined$seurat_clusters, colors = cb_palette.use,size=1)
fig

