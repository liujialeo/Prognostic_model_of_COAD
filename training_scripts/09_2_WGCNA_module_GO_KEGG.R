rm(list = ls())
library(clusterProfiler)
library(tidyr)
library(dplyr)

genes=read.table("output/CytoscapeInput-nodes-turquoise.txt",header=T,sep="\t")
genes <- genes%>%
  distinct(nodeName,.keep_all = T)
genes <- as.data.frame(genes[,1])
colnames(genes) <- "gene"
rownames(genes) <- genes[,1]
### 这个分析需要什么数据？
### 获得基因列表
gene <- rownames(genes)
#基因名称转换，返回的是数据框
gene = bitr(gene, 
            fromType="SYMBOL", 
            toType="ENTREZID", 
            OrgDb="org.Hs.eg.db")
head(gene)
#GO分析的细胞组分 CC
ego_CC <- enrichGO(gene = gene$ENTREZID,
                   OrgDb= org.Hs.eg.db,
                   ont = "CC",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.01,
                   readable = TRUE)

#**作图**
#条带图
barplot(ego_CC)
# 点图
dotplot(ego_CC)
# GO 作图
goplot(ego_CC)


#GO分析的生物过程BP
ego_BP <- enrichGO(gene = gene$ENTREZID,
                   OrgDb= org.Hs.eg.db,
                   ont = "BP",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.01,
                   readable = TRUE)

#**作图**
#条带图
barplot(ego_BP)
# 点图
dotplot(ego_BP)
# GO 作图
goplot(ego_BP)

#GO分析分子功能MF：

ego_MF <- enrichGO(gene = gene$ENTREZID,
                   OrgDb= org.Hs.eg.db,
                   ont = "MF",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.01,
                   readable = TRUE)

#**作图**
#条带图
barplot(ego_MF)
# 点图
dotplot(ego_MF)
# GO 作图
goplot(ego_MF)

############################################
### 尝试新的GO聚类的方法
go <- enrichGO(gene = gene$ENTREZID, OrgDb = "org.Hs.eg.db", ont="all")
#save(go,file = "output/go.Rdata")
library(ggplot2)
p <- barplot(go, split="ONTOLOGY") +facet_grid(ONTOLOGY~., scale="free")
p <- dotplot(go, split="ONTOLOGY") +facet_grid(ONTOLOGY~., scale="free")
p

############################################
### 还可以使用enrich函数直接分析
### 使用msigdbr 获取基因集
### https://cran.r-project.org/web/packages/msigdbr/vignettes/msigdbr-intro.html
library(msigdbr)
msigdbr_species()
msigdbr_collections()
library(dplyr)
h_df <- msigdbr(species = "Homo sapiens") %>% 
  ## 筛选条目
  filter(gs_cat == "C5",gs_subcat == "GO:BP")
colnames(h_df)


h_df <- msigdbr(species = "Homo sapiens") %>% 
  ## 筛选条目
  filter(gs_cat == "C5",gs_subcat == "GO:BP") %>% 
  ## 筛选合适的列
  dplyr::select(gs_name,gene_symbol)
go_bp <- enricher(gene = gene$SYMBOL,TERM2GENE = h_df)
dotplot(go_bp)

############################################
### 试试kegg
h_df <- msigdbr(species = "Homo sapiens") %>% 
  ## 筛选条目
  filter(gs_cat == "C2",gs_subcat == "CP:KEGG") %>% 
  ## 筛选合适的列
  dplyr::select(gs_name,gene_symbol)

kegg <- enricher(gene = gene$SYMBOL,TERM2GENE = h_df)
dotplot(kegg)

############################################
### 多组同时做
### 比如: 上调基因，下调基因，所有差异基因
### 函数 compareCluster
rm(list = ls())
library(clusterProfiler)
library(dplyr)
library(tibble)

load(file = "data/allDiff.Rdata")

genes <- allDiff %>% 
  rownames_to_column("gene") %>% 
  filter(adj.P.Val < 0.05) %>% 
  arrange(desc(logFC))

upgeneall <- genes$gene[genes$logFC>0]
dngeneall <- genes$gene[genes$logFC<0]

## 选取一定数目上调基因
upgene = head(upgeneall,500)
## 选取一定数目下调基因
dngene = tail(dngeneall,500)

## 输入文件
gcSample = list(up = upgene,down = dngene,all = c(upgene,dngene))

## 基因集 
h_df <- msigdbr(species = "Homo sapiens") %>% 
  filter(gs_cat == "C2",gs_subcat == "CP:KEGG") %>% 
  dplyr::select(gs_name,gene_symbol)

## 万能代码！！
xx <- compareCluster(gcSample, fun="enricher",TERM2GENE = h_df)
## 作图
dotplot(xx)