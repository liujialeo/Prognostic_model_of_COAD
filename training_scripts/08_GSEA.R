rm(list = ls())

##############################################################
### GSEA 分析
load(file = "output/COAD_allDiff.Rdata")
library(dplyr)
gene_df <- allDiff%>%
  filter(gene_type =="protein_coding" )%>% 
  dplyr::select(gene_id,logFC,gene) %>% 
  ## 去掉NA
  filter(gene!="") %>% 
  ## 去掉重复
  distinct(gene,.keep_all = T)

### 1.获取基因logFC
geneList <- gene_df$logFC
### 2.命名
names(geneList) = gene_df$gene
## 3.排序很重要
geneList = sort(geneList, decreasing = TRUE)

head(geneList)
library(clusterProfiler)

## 读入hallmarks gene set，从哪来？
hallmarks <- read.gmt("resource/h.all.v7.1.symbols.gmt")
#hallmarks <- read.gmt("resource/msigdb.v7.5.1.symbols.gmt")

### GSEA 分析
gseahallmarks <- GSEA(geneList,TERM2GENE =hallmarks)

library(ggplot2)
yd <- as.data.frame(gseahallmarks)

dotplot(gseahallmarks,showCategory=30,split=".sign")+facet_grid(~.sign)

library(enrichplot)
pathway.id = "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
gseaplot2(gseahallmarks, 
          color = "red",
          geneSetID = pathway.id,
          pvalue_table = T)

pathway.id = "HALLMARK_OXIDATIVE_PHOSPHORYLATION"
gseaplot2(gseahallmarks, 
          color = "red",
          geneSetID = pathway.id,
          pvalue_table = T)

##展示多个通路在一张图上
hallmarks_order <- gseahallmarks[order(gseahallmarks$NES,decreasing=T)] #按enrichmentScore或NES降序排列
gseaplot2(gseahallmarks, 
          color = "red",
          geneSetID =rownames(hallmarks_order)[1:5],
          pvalue_table = T)
#将结果画成山脊图
ridgeplot(gseahallmarks,showCategory = 30, fill = "pvalue")
### 转录组数据还可以参考这个帖子
### https://mp.weixin.qq.com/s/7BBJGTlOa5i6YPMlrRqw2w
