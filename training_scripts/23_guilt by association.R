
rm(list = ls())
### 加载表达数据，这是vst标准化后的数据
load(file = "output/exprSet_heatmap.Rdata")
test <- exprSet[1:10,1:10]
##1.设定容器
correlation <- data.frame()
##2.准备数据
data <- as.data.frame(t(exprSet)) 
test <- data[1:10,1:10]
##3.获取基因列表
genelist <- colnames(data)
##4.指定基因
gene <- "UCN"
genedata <- as.numeric(data[,gene])
##5.开始for循环
for(i in 1:length(genelist)){
  ## 1.指示
  print(i)
  ## 2.计算
  dd = cor.test(genedata,as.numeric(data[,i]),method="spearman")
  ## 3.填充
  correlation[i,1] = gene
  correlation[i,2] = genelist[i]
  correlation[i,3] = dd$estimate
  correlation[i,4] = dd$p.value
}

colnames(correlation) <- c("gene1","gene2","cor","p.value")
save(correlation, file = "./output/UCN_association.Rdata")
## 排序
library(dplyr)

genelist = correlation$cor
names(genelist) = correlation$gene2
genelist = sort(genelist,decreasing = T)
head(genelist)
##############################
## 基因集 
library(msigdbr)
h_df <- msigdbr(species = "Homo sapiens") %>% 
  filter(gs_cat == "H") %>% 
  dplyr::select(gs_name,gene_symbol)

#################################
## enrichr
## 选取一定数目上调基因
upgene = head(names(genelist),1000)
## 选取一定数目下调基因
dngene = tail(names(genelist),1000)
## 输入文件
gcSample = list(up = upgene,down = dngene,all = c(upgene,dngene))

## 万能代码！！
library(clusterProfiler)
xx <- compareCluster(gcSample, fun="enricher",TERM2GENE = h_df)
## 作图
dotplot(xx)

#################################
## GSEA
y <- GSEA(genelist,TERM2GENE = h_df)
yd <- as.data.frame(y)
library(ggplot2)
dotplot(y,showCategory=15,split=".sign")+facet_grid(~.sign)
### 自定义画图
ggplot(y, showCategory = 20, aes(NES, forcats::fct_reorder(Description, NES))) + 
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_colour_gradientn(colours=c("#7fc97f","#984ea3")) +
  #scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
  #scale_color_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE))+
  scale_size_continuous(range=c(2, 10)) +
  theme_bw() + 
  xlab("Normalized Enrichment Score") +
  ylab(NULL)

ggsave(filename = "./output/UCN_ssociation_GSEAplot.pdf",width = 8,height = 9)

library(enrichplot)
library(export)
pathway.id = "HALLMARK_MYC_TARGETS_V2"
gseaplot2(y,color = "green",geneSetID = pathway.id,pvalue_table = T)
# ggsave(filename = "./output/UCN_aasociation_GSEA_MYC_TARGETS_V2.pdf")
graph2pdf(file="./output/GABRD_aasociation_GSEA_MYC_TARGETS_V2.pdf")

pathway.id = "HALLMARK_MYC_TARGETS_V1"
gseaplot2(y,color = "green",geneSetID = pathway.id,pvalue_table = T)
# ggsave(filename = "./output/UCN_aasociation_GSEA_MYC_TARGETS_V1.pdf",
#                   width = 12, height = 8)
graph2pdf(file="./output/UCN_aasociation_GSEA_MYC_TARGETS_V1.pdf")

pathway.id = "HALLMARK_DNA_REPAIR"
gseaplot2(y,color = "green",geneSetID = pathway.id,pvalue_table = T)
# ggsave(filename = "./output/UCN_aasociation_GSEA_DNA_REPAIR.pdf",
#                   width = 12, height = 8)
graph2pdf(file="./output/UCN_aasociation_GSEA_DNA_REPAIR.pdf")

gseaplot2(y, geneSetID = c(1:3))
ggsave(filename = "./output/UCN_aasociation_GSEA.pdf")
### cutting edge作图
if(!requireNamespace("ggnewscale",quietly = TRUE)){
  options("repos"=c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
  install.packages("ggnewscale",update = F,ask = F)
}
cnetplot(y,showCategory = 4,foldChange = genelist,colorEdge = T)
ggsave(filename = "./output/UCN_aasociation_GSEA_cnetplot.pdf",width = 15,height = 13)

####提取特定通路的基因
GSEAset_genes <- yd$core_enrichment[3]
GSEAset_genes
GSEAset_genes <- as.data.frame(unlist(strsplit(GSEAset_genes,split = "/"))) 
colnames(GSEAset_genes) <- "gene2"
MYC_target1_genes <- inner_join(correlation,GSEAset_genes, by = "gene2") 
MYC_target1_genes <- arrange(MYC_target1_genes, desc(cor))
###批量输出特定通路的基因列表
dir.create("./output/GSEA_out")
for (i in 1:nrow(yd)) {
  print(i)
  GSEAset = yd$ID[i]
  GSEAset_genes = yd$core_enrichment[i]
  data = as.data.frame(unlist(strsplit(GSEAset_genes,split = "/")))
  write.table(data,file = paste0("output/GSEA_out/",GSEAset,".csv") ,sep = ",",row.names = F)
  
}


## R 包和网站 GTBA
## 如何用?
## 单基因批量相关性分析的妙用
## https://mp.weixin.qq.com/s/TfE2koPhSkFxTWpb7TlGKA
## 单基因批量相关性分析的GSEA
## https://mp.weixin.qq.com/s/sZJPW8OWaLNBiXXrs7UYFw