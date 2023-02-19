library(Seurat)
library(dplyr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)
library(DOSE)

#keytypes(org.Hs.eg.db)


#挑选每个细胞亚群中特意高表达的前200个基因
#sample.markers <- FindAllMarkers(scRNA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top200 <- all.markers %>% group_by(cluster) %>% top_n(n = 200, wt = avg_log2FC) 
top200_gene<-unstack(top200,gene~cluster)
names(top200_gene)=new.cluster.ids

#将gene symbol转换成entriz gene id
#做KEGG富集分析需要
top200_entrez <- lapply(X = top200_gene, FUN = function(x) {
  x <- bitr(x, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  x <- x[,-1]
})

######################################
#对所有亚群做GO富集分析 Gene ontoloty
#####################################
######################################
# CC cellular component
######################################
all_cc <- compareCluster(top200_entrez, 
                         fun='enrichGO',
                         ont= 'CC',
                         OrgDb='org.Hs.eg.db')
write.csv(as.data.frame(all_cc),paste0("./",sam.name,"/all_cc.csv"))

#绘制CC富集图
#includeAll=FALSE 只包含前五个term
pdf(paste0("./",sam.name,"/all_GO_CC.pdf"),width = 15, height = 12)
dotplot(all_cc, includeAll=FALSE)+ theme(axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()

##################################################
# MF molecular function
#################################################
all_mf <- compareCluster(top200_entrez, 
                         fun='enrichGO',
                         ont= 'MF',
                         OrgDb='org.Hs.eg.db')
write.csv(as.data.frame(all_mf),paste0("./",sam.name,"/all_mf.csv"))

#绘制MF富集图
#includeAll=FALSE 只包含前五个term
pdf(paste0("./",sam.name,"/all_GO_MF.pdf"),width = 15, height = 12)
dotplot(all_mf, includeAll=FALSE)+ theme(axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()
#############################################
# BP biological process
#############################################
all_bp <- compareCluster(top200_entrez, 
                         fun='enrichGO',
                         ont= 'BP',
                         OrgDb='org.Hs.eg.db')
write.csv(as.data.frame(all_bp),paste0("./",sam.name,"/all_bp.csv"))

#绘制BP富集图
#includeAll=FALSE 只包含前五个term
pdf(paste0("./",sam.name,"/all_GO_BP.pdf"),width = 15, height = 12)
dotplot(all_bp, includeAll=FALSE)+ theme(axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()



#################################
# 对某一个特定的亚群做GO富集分析
################################
all_go_B  <- enrichGO(gene = top200_entrez$`0`,   #可以修改细胞亚群类型
                          OrgDb = 'org.Hs.eg.db',
                          ont = 'ALL',
                          pvalueCutoff = 0.05,
                          qvalueCutoff =0.2)

all_go_B<-setReadable(all_go_B, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.csv(as.data.frame(all_go_B),paste0("./",sam.name,"/all_go_B.csv"))

pdf(paste0("./",sam.name,"/B_all_GO.pdf"),width = 15, height = 12)
dotplot(all_go_B, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")
dev.off()

###################################
#对所有亚群做KEGG富集分析
###################################
all_kegg <- compareCluster(top200_entrez,
                           fun='enrichKEGG',
                           #pvalueCutoff=0.05, 
                           organism="hsa"
                           )

write.csv(as.data.frame(all_kegg),paste0("./",sam.name,"/all_kegg.csv"))

pdf(paste0("./",sam.name,"/all_KEGG.pdf"),width = 15, height = 12)
dotplot(all_kegg, showCategory=5, includeAll=FALSE)+ theme(axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()

#################################
#对某一个特定的亚群做KEGG富集分析
#################################
kegg_B <- enrichKEGG(gene = top200_entrez$`0`,
                             organism  = 'hsa', 
                             pvalueCutoff = 0.05,
                             qvalueCutoff =0.2)


kegg_B<-setReadable(kegg_B, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.csv(kegg_B,paste0("./",sam.name,"/kegg_B.csv"))
pdf(paste0("./",sam.name,"/KEGG_B.pdf"),width = 12, height = 6)
dotplot(kegg_B)
dev.off()


