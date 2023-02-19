### 本节任务：如何展示一群基因的表达结果
### 热图！
###############
### heatmap热图
#用行名提取数据
rm(list = ls())
### 加载表达数据，这是vst标准化后的数据
load(file = "output/exprSet_heatmap.Rdata")
### 加载差异基因列表
load(file = "output/COAD_allDiff.Rdata")
### 筛选差异基因
library(dplyr)
diffgene <- allDiff %>% 
  filter(gene !="") %>% 
  filter(adj.P.Val < 0.01) %>% 
  filter(abs(logFC) >2)
save(diffgene,file = "output/diffgene.Rdata")

### 加载热图的R包
library(pheatmap)
### 用名称提取部分数据用作热图绘制
heatdata <- exprSet[diffgene$gene,]

### 制作一个分组信息用于注释
load(file = "output/COAD_metadata.Rdata")
colnames(metadata) <- c("sample","group")

annotation_col <- data.frame(group=metadata$group)
rownames(annotation_col) <- metadata$sample

### 直接作图
pheatmap(heatdata)

### 如果注释出界, 可以通过调整格子比例和字体修正
pheatmap(heatdata, #热图的数据
         cluster_rows = TRUE,#行聚类
         cluster_cols = TRUE,#列聚类，可以看出样本之间的区分度
         annotation_col =annotation_col, #标注样本分类
         annotation_legend=TRUE, # 显示注释
         show_rownames = F,# 显示行名
         show_colnames = F,# 显示行名
         scale = "row", #以行来标准化，这个功能很不错
         color =colorRampPalette(c("blue", "white","red"))(100),#调色
         #filename = "heatmap_F.pdf",#是否保存
         cellwidth = 0.2, cellheight = 0.5,# 格子比例
         fontsize = 10)

### 加载R包
library(export)
### 导成PPT可编辑的格式
graph2ppt(file="output/heatmap.pptx")


################################################################
### 如何展示所有的基因变化
### 火山图
rm(list = ls())
### 加载差异基因列表
load(file = "output/COAD_allDiff.Rdata")
### 筛选差异基因
library(dplyr)
library(ggplot2)
library(ggrepel)

dat <- allDiff %>% 
  filter(gene !="") %>% 
    filter(gene_type =="protein_coding" )
#加change列,标记上下调基因
logFC_t=2
P.Value_t = 0.01
k1 = (dat$P.Value < P.Value_t)&(dat$logFC < -logFC_t)
k2 = (dat$P.Value < P.Value_t)&(dat$logFC > logFC_t)
change = ifelse(k1,"down",ifelse(k2,"up","stable"))
dat <- mutate(dat,change)


p <- ggplot(data = dat, 
            aes(x = logFC, 
                y = -log10(P.Value))) +
  geom_point(alpha=0.4, size=3.5, 
             aes(color=change)) +
  ylab("-log10(Pvalue)")+
  scale_color_manual(values=c("deepskyblue", "black","orangered"))+
  geom_vline(xintercept=c(-logFC_t,logFC_t),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(P.Value_t),lty=4,col="black",lwd=0.8) +
  theme_bw()
p

# data$significant <- as.factor(data$P.Value<0.01 & abs(data$logFC) > 2)
# ggplot(data=data, aes(x=logFC, y =-log10(P.Value),color=significant)) +
#   geom_point(alpha=0.8, size=1.2,col="black")+
#   geom_point(data=subset(data, logFC > 2),alpha=2, size=1.4,col="orangered")+
#   geom_point(data=subset(data, logFC < -2),alpha=2, size=1.4,col="deepskyblue")+
#   labs(x="log2 (fold change)",y="-log10 (adj.P.Val)")+
#   theme(plot.title = element_text(hjust = 0.4))+
#   geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.6,alpha=0.8)+
#   geom_vline(xintercept = c(-2,2),lty=4,lwd=0.6,alpha=0.8)+
#   theme_bw()+
#   theme(panel.border = element_blank(),
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),   
#         axis.line = element_line(colour = "black")) #+
#   # geom_point(data=subset(data, abs(logFC) >= 4),alpha=0.8, size=3,col="green")+
#   # geom_text_repel(data=subset(data, abs(logFC) > 4),
#   #                 aes(label=gene),col="black",alpha = 0.8)

ggsave(filename = "./output/volcano.pdf")
