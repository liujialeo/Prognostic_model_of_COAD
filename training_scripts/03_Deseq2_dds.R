################################################
################################################
### 作者：果子
### 更新时间：2020-01-04
### 微信公众号:果子学生信
### 私人微信：guotosky
### 个人博客: https://codingsoeasy.com/

rm(list = ls())
load(file = "./output/COAD_RNASEQ_exprdf.Rdata")
class(expr_df)

exprSet <- as.data.frame(expr_df)
### 查看分组
### https://dwz.cn/WVgQUqfw
### 样本名称
TCGA_id <- colnames(exprSet)[-1]
table(substring(TCGA_id,14,15))
### 我们发现了7个转移的样本，本次分析，我们关注的是癌症和癌旁，先把转移的样本去掉
### 原发和转移的对比作为家庭作业

TCGA_id <- TCGA_id[substring(TCGA_id,14,15)!="06"]
TCGA_id <- TCGA_id[substring(TCGA_id,14,15)!="02"]

exprSet <- cbind(exprSet$gene_id,exprSet[,TCGA_id])
TCGA_id <- colnames(exprSet)[-1]
table(substring(TCGA_id,14,15))

### 创建metadata
sample <- ifelse(substring(TCGA_id,14,15)=="01","cancer","normal")
sample <- factor(sample,levels = c("normal","cancer"),ordered = F)
metadata <- data.frame(TCGA_id,sample) 
save(metadata,file = "./output/COAD_metadata.Rdata")
library(DESeq2)
### 第一列有名称，所以tidy=TRUE
dds <-DESeqDataSetFromMatrix(countData=exprSet, 
                             colData=metadata, 
                             design=~sample,
                             tidy=TRUE)
nrow(dds)
### 如果一个基因在所有样本中的counts数小于等于1，我们就把他删掉
dds <- dds[rowSums(counts(dds))>1,]
nrow(dds)

### 数据标准化用于看聚类
### https://dwz.cn/xJTuI4aO
### 很耗时间
if(T){
  vsd <- vst(dds, blind = FALSE)
}
save(vsd,file = "output/COAD_vsd.Rdata")
load(file = "output/COAD_vsd.Rdata")
### PCAf
plotPCA(vsd, "sample")
### 保存数据用于热图
exprSet_vst <- as.data.frame(assay(vsd))
test <- exprSet_vst[1:10,1:10]
save(exprSet_vst,file = "output/COAD_exprSet_vst.Rdata")

### 最困难的一步来了.燃烧你的小电脑。
### 这里还使用了并行化处理来解决速度的问题
### Deseq2 更新后速度大幅度提升
### https://dwz.cn/bb1jDs12

library(BiocParallel)
dds <- DESeq(dds,parallel = T)
save(dds,file="output/dds_very_long.Rdata")

load(file="output/dds_very_long.Rdata")

###################################################################
### 提取标准化后的数据，注意，我们data warngling 用的就是这个数据
normalized_counts <- as.data.frame(counts(dds, normalized=TRUE))
save(normalized_counts,file = "output/COAD_normalized_counts.Rdata")
### 单独作图简易展示
plotCounts(dds, gene = "ENSG00000000003.13", intgroup=c("sample"))
### 还可以把数据返回
plotdata <- plotCounts(dds, gene = "ENSG00000000003.13", intgroup=c("sample","paire_info"),returnData = T)
library(ggplot2)
ggplot(plotdata,aes(x=sample,y=count,col=sample))+
  geom_jitter()+
  theme_bw()

#######################################################################
### logFC矫正,注意顺序哈，?号
contrast <- c("sample","cancer","normal")
dd1 <- results(dds, contrast=contrast, alpha = 0.05)
plotMA(dd1, ylim=c(-5,5))
### 发现样本需要做logFC校正
###dd2 <- lfcShrink(dds, contrast=contrast, res=dd1)
###plotMA(dd2, ylim=c(-5,5))

summary(dd1, alpha = 0.05)
library(dplyr)
library(tibble)
library(tidyr)
### 导出差异分析的结果
res <- dd1 %>% 
  data.frame() %>% 
  rownames_to_column("gene_id") 

### 基因注释
## 准备注释文件
### https://gdc.cancer.gov/about-data/data-harmonization-and-generation/gdc-reference-files
### 下载后先解压gencode.v22.annotation.gtf.gz
### gtf1 <- rtracklayer::import('gencode.v22.annotation.gtf')
### gtf_df <- as.data.frame(gtf1)
### save(gtf_df,file = "gtf_df.Rdata")
library(AnnotationDbi)
library(org.Hs.eg.db)
load(file = "./resource/gtf_df.Rdata")
## 注释差异基因
allDiff <- gtf_df %>% 
  ##筛选gene
  dplyr::filter(type=="gene") %>%
  ## 筛选基因名称和类型
  dplyr::select(c(gene_name,gene_id,gene_type)) %>% 
  ## 合并
  dplyr::inner_join(res,by ="gene_id")

term=as.character(allDiff$gene_id[1]) 
### 1.裂解
term = unlist(strsplit(term, split=".", fixed=T))[1]
term
##批量操作
allDiff$gene_id <- as.character(allDiff$gene_id)
for (i in 1:nrow(allDiff)) {
  print(i)
  term = allDiff$gene_id[i]
  term = unlist(strsplit(term, split=".", fixed=T))[1]
  allDiff$gene_id[i] = term
}
### 增加ENTREZID
allDiff$entrez <- mapIds(org.Hs.eg.db,
                     keys=allDiff$gene_id,
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

### 修改logFC，p值的名称，为的是跟火山图的代码匹配
colnames(allDiff) <- c("gene","gene_id","gene_type","baseMean","logFC","lfcSE","stat","P.Value","adj.P.Val","entrez")


save(allDiff,file = "./output/COAD_allDiff.Rdata")

### 现在有了差异基因列表，热图，火山图，GO，KEGG，GSEA，pathview都可以做了
### 自行完成

