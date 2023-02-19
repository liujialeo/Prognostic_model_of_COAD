################################################
################################################
### 作者：果子
### 更新时间：2020-01-04
### 微信公众号:果子学生信
### 私人微信：guotosky
### 个人博客: https://codingsoeasy.com/

### 分别提取mRNA和lncRNA
rm(list = ls())
load(file = "output/COAD_normalized_counts.Rdata")
expr_df <- cbind(gene_id= rownames(normalized_counts),normalized_counts,stringsAsFactors=F)
## 准备注释文件
### https://gdc.cancer.gov/about-data/data-harmonization-and-generation/gdc-reference-files
### 下载后先解压gencode.v22.annotation.gtf.gz
### gtf1 <- rtracklayer::import('gencode.v22.annotation.gtf')
### gtf_df <- as.data.frame(gtf1)
### save(gtf_df,file = "gtf_df.Rdata")

load(file = "./resource/gtf_df.Rdata")
test <- gtf_df[1:100,]

library(dplyr)
library(tidyr)
## 提取mRNA
mRNA_exprSet <- gtf_df %>% 
  #筛选gene,和编码指标
  dplyr::filter(type=="gene",gene_type=="protein_coding") %>%
  dplyr::select(c(gene_name,gene_id,gene_type)) %>% 
  dplyr::inner_join(expr_df,by ="gene_id") %>% 
  tidyr::unite(gene_id,gene_name,gene_id,gene_type,sep = "|")
test <- mRNA_exprSet[1:10,1:10]


### lncRNA
### 如何定义非编码RNA呢？
### https://www.gencodegenes.org/pages/biotypes.html
ncRNA <- c("3prime_overlapping_ncrna","antisense","lincRNA",
           "macro_lncRNA","non_coding",
           "processed_transcript","sense_intronic",
           "sense_overlapping")

LncRNA_exprSet <- gtf_df %>% 
  dplyr::filter(type=="gene",gene_type %in% ncRNA) %>% 
  dplyr::select(c(gene_name,gene_id,gene_type)) %>% 
  dplyr::inner_join(expr_df,by ="gene_id") %>% 
  tidyr::unite(gene_id,gene_name,gene_id,gene_type,sep = "|")

test <- LncRNA_exprSet[1:10,1:10]


### 可以集体合并，同时提取编码和非编码的信息
mytype<- c("protein_coding","3prime_overlapping_ncrna","antisense","lincRNA",
           "macro_lncRNA","non_coding",
           "processed_transcript","sense_intronic",
           "sense_overlapping")
all_exprSet <- gtf_df %>% 
  dplyr::filter(type=="gene",gene_type %in% mytype) %>% 
  dplyr::select(c(gene_name,gene_id,gene_type)) %>% 
  dplyr::inner_join(expr_df,by ="gene_id") %>% 
  tidyr::unite(gene_id,gene_name,gene_id,gene_type,sep = "|")
save(all_exprSet,file = "output/COAD_all_exprSet.Rdata")


#########################################################################
### 那现在能干什么呢？
### 1.调整数据格式到清洁数据对接ggplot2
### 三大步，基因注释，行列转换，分组信息
##基因注释
###mRNA
View(as.data.frame(mRNA_exprSet$gene_id) )  
###转化为纯基因名
##拿一个举例子
term=as.character(mRNA_exprSet$gene_id[1]) 
### 1.裂解
term = unlist(strsplit(term, split="|", fixed=T))[1]
term
##批量操作
mRNA_exprSet$gene_id <- as.character(mRNA_exprSet$gene_id)
for (i in 1:nrow(mRNA_exprSet)) {
  print(i)
  term = mRNA_exprSet$gene_id[i]
  term = unlist(strsplit(term, split="|", fixed=T))[1]
  mRNA_exprSet$gene_id[i] = term
}
##第一列转行名，按照表达量大小去重
mRNA_exprSet <- mRNA_exprSet %>% 
  mutate(newcolumn = rowMeans(.[,-1])) %>% 
  arrange(desc(newcolumn)) %>% 
  distinct(gene_id,.keep_all = T) %>% 
  dplyr::select(-newcolumn)
rownames(mRNA_exprSet) <- mRNA_exprSet[,1]
mRNA_exprSet <- mRNA_exprSet[,-1]
save(mRNA_exprSet,file = "./output/COAD_mRNA_exprSet.Rdata")

##LINCRNA
View(as.data.frame(LncRNA_exprSet$gene_id) )  
###转化为纯基因名
##拿一个举例子
term=as.character(LncRNA_exprSet$gene_id[1]) 
### 1.裂解
term = unlist(strsplit(term, split="|", fixed=T))[1]
term
##批量操作
LncRNA_exprSet$gene_id <- as.character(LncRNA_exprSet$gene_id)
for (i in 1:nrow(LncRNA_exprSet)) {
  print(i)
  term = LncRNA_exprSet$gene_id[i]
  term = unlist(strsplit(term, split="|", fixed=T))[1]
  LncRNA_exprSet$gene_id[i] = term
}
##第一列转行名，按照表达量大小去重
LncRNA_exprSet <- LncRNA_exprSet %>% 
  mutate(newcolumn = rowMeans(.[,-1])) %>% 
  arrange(desc(newcolumn)) %>% 
  distinct(gene_id,.keep_all = T) %>% 
  dplyr::select(-newcolumn)
rownames(LncRNA_exprSet) <- LncRNA_exprSet[,1]
LncRNA_exprSet <- LncRNA_exprSet[,-1]
save(LncRNA_exprSet,file = "./output/COAD_lincRNA_exprSet.Rdata")

### 2.单基因的批量相关性分析
### 参考a:https://dwz.cn/ebZiHkEK
### 参考b:https://dwz.cn/nQqbZQaH
### 3.单基因GSEA分析，神器
### 可以注释任意基因
### 核心点，在于用单个基因的表达量把样本分成高表达和低表达
### 之后就变成了两组RNA-seq分析，接着进行GSEA分析