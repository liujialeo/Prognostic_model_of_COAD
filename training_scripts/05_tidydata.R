rm(list = ls())
### 本节任务: 转录组数据变成清洁数据
###################################################################
### 一个基因，或者多个基因的差异作图
### 需要清洁数据
load("output/COAD_exprSet_vst.Rdata")
load("output/COAD_allDiff.Rdata")
load("output/COAD_metadata.Rdata")
### 1.获取注释文件
library(dplyr)
ensemble_symbol = allDiff %>% 
  dplyr::select(gene_id,gene) %>% 
  filter(gene !="") 

### 2.准备表达量数据
### 行名变成列
term=as.character(rownames(exprSet_vst)[1])
### 1.裂解
term = unlist(strsplit(term, split=".", fixed=T))[1]
term

rownames(exprSet_vst) <- as.character(rownames(exprSet_vst))
for (i in 1:nrow(exprSet_vst)) {
  print(i)
  term = rownames(exprSet_vst)[i]
  term = unlist(strsplit(term, split=".", fixed=T))[1]
  rownames(exprSet_vst)[i] = term
}

exprSet <- cbind(gene_id=rownames(exprSet_vst),exprSet_vst)

### 3.交叉合并
exprSet <- merge(ensemble_symbol,exprSet,by="gene_id")

### 4.基因名称去重(保留最大值法)
### 列转行名一定要去重，因为行名不支持重复
exprSet <- exprSet %>% 
  dplyr::select(-gene_id) %>% 
  mutate(newcolumn = rowMeans(.[,-1])) %>% 
  arrange(desc(newcolumn)) %>% 
  distinct(gene,.keep_all = T) %>% 
  dplyr::select(-newcolumn)

### 5.列变成行名
rownames(exprSet) <- exprSet[,1]
exprSet <- exprSet[,-1]
### 保存数据，画热图
save(exprSet,file = "output/exprSet_heatmap.Rdata")

### 6.行列转置
exprSet <- t(exprSet)
exprSet <- as.data.frame(exprSet)
test <- exprSet[,1:10]

### 7.添加分组
colnames(metadata) <- c("sample","group")
exprSet <- cbind(group= metadata$group,exprSet)
test <- exprSet[,1:10]
save(exprSet,file = "output/exprSet_tidy.Rdata")

###########################################################

rm(list = ls())
library(RColorBrewer)
display.brewer.all()
###获取自定义颜色
colset <- brewer.pal(10,"Paired") ##brewer.pal选择颜色的最小种类是3
load(file = "output/exprSet_tidy.Rdata")

## steal plot
my_comparisons <- list(
  c("cancer", "normal")
)
library(ggpubr)
ggboxplot(
  exprSet, x = "group", y = "GABRD",
  color = "group", palette = c("#00AFBB", "#E7B800"),
  add = "jitter"
)+
  stat_compare_means(comparisons = my_comparisons, method = "t.test")

## 改写成函数
diffplot <- function(gene){
  my_comparisons <- list(
    c("treat", "con")
  )
  library(ggpubr)
  ggboxplot(
    exprSet, x = "group", y = gene,
    color = "group", palette = c("#00AFBB", "#E7B800"),
    add = "jitter"
  )+
    stat_compare_means(comparisons = my_comparisons, method = "t.test")
}

### AGR3,ESR1,SLC4A10, ALPP,VSIR,PLA2G2F
diffplot("NFKBIE")
diffplot("ALPP")

## 多个基因作图查看
## 先把基因提取出来
genelist <- c("FAM132B","DLX4","UCN","GABRD","MS4A2","DYNC1I1",
              "GRIK3","VWC2","LEP","SNCB","CDH10","GABRG1")
## 再提取表达量，使用名称选取行
data <- exprSet[,c("group",genelist)]
## 用pivot_longer调整数据，数据变长，增加的是行
library(tidyr)
data <- data %>% 
  pivot_longer(cols=-1,
               names_to= "gene",
               values_to = "expression")
## 多基因作图
## 作图
ggplot(data = data,aes(x=group,y=expression,fill=group))+
  geom_boxplot()+
  geom_jitter(size= 1)+
  theme_bw()+
  #facet_grid(.~gene)+
  facet_wrap(~gene,nrow = 1)+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test")+
  scale_fill_manual(values= colset[c(1,2)] )

ggsave(filename = "./output/multi_genes_multivariateCOX.pdf", width = 15, height = 5 )
