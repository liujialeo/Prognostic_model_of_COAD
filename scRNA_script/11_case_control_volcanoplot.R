library(ggplot2)
library(EnhancedVolcano)
#### Volcano plot
DefaultAssay(obj.combined) <- "RNA"
Idents(obj.combined) <- "celltype.orig.ident"
Colonic_stem_cell.diff <- FindMarkers(obj.combined, ident.1 = "Colonic stem cell_adenoma", ident.2 = "Colonic stem cell_normal", verbose = FALSE)

low<-floor(range(Colonic_stem_cell.diff$avg_log2FC)[1])
high<-ceiling(range(Colonic_stem_cell.diff$avg_log2FC)[2])
EnhancedVolcano(Colonic_stem_cell.diff,
                title = 'Colonic_stem_cell.diff adenoma versus normal',
                lab = rownames(Colonic_stem_cell.diff),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                FCcutoff = 1,
                xlim = c(low, high))

#DEA函数
#寻找某一个特定亚群在case和control组中差异表达的基因
#绘制火山图
DEA <- function(x,table_folder,figure_folder){
  diff <- FindMarkers(obj.combined, 
                      ident.1 = paste0(x,"_","adenoma"), 
                      ident.2 = paste0(x,"_","normal"), 
                      verbose = FALSE)
  filename=paste0(table_folder,"/Diff_exp_genes_",x,".csv")
  write.csv(file=filename,diff)
  
  low<-floor(range(diff$avg_log2FC)[1])
  high<-ceiling(range(diff$avg_log2FC)[2])
  pdf(file=paste0(figure_folder,"/",x,"_volcanoplot.pdf"),width=12)
  print(EnhancedVolcano(diff,
                        title = paste0(x,' adenoma versus normal'),
                        lab = rownames(diff),
                        x = 'avg_log2FC',
                        y = 'p_val_adj',
                        FCcutoff = 1,
                        xlim = c(low, high)
  ))
  dev.off()
}


#循环，对每一个亚群做差异表达分析并绘制热图
dir.create("DEG")
dir.create("volcanoplot")


celltypes<-levels(obj.combined$celltype)
for (cluster in celltypes){
  DEA(cluster,"DEG","volcanoplot")
  }


