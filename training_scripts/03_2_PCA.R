rm(list = ls())
library(BiocManager)
library(devtools)
#install_github("GuangchuangYu/yyplot")
#install_github('fawda123/ggord')
library(ggplot2)
library(plyr)
library(ggord)
library(yyplot)
library(ggsci)

library(SummarizedExperiment)
load(file = "output/COAD_vsd.Rdata")
exprSet_df <- as.data.frame(assay(vsd))
exprSet_df <- as.data.frame(t(exprSet_df)) 
test <- exprSet_df[1:5,1:5]

### 创建metadata
TCGA_id <- rownames(exprSet_df)
table(substring(TCGA_id,14,15))
sample <- ifelse(substring(TCGA_id,14,15)=="01","cancer","normal")
sample <- factor(sample,levels = c("normal","cancer"),ordered = F)
meta_df1 <- data.frame(TCGA_id,sample) 


rownames(meta_df1) <- meta_df1[,1]
colnames(meta_df1) <- c("sample","group")
meta_df <- as.data.frame(meta_df1[,2]) 
colnames(meta_df) <- "group"
rownames(meta_df) <- rownames(meta_df1)


#用`prcomp`进行PCA分析
pca.results <- prcomp(exprSet_df, center = TRUE, scale. = FALSE)

#定义足够多的颜色，用于展示分组
#mycol <- c("#223D6C","#D20A13","#088247","#FFD121","#11AA4D","#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767")
mycol <- c("deeppink","orangered","deepskyblue")
########经典版

#有可能你的网络在安装yyplot时遇到困难，我把geom_ord_ellipse函数单独下载，通过本地进行加载。
#调用yyplot包里的geom_ord_ellipse函数
source('resource/geom_ord_ellipse.R') #该文件位于当前文件夹

#用ggord画基本PCA图

p <- ggord(pca.results, grp_in = meta_df$group, repel=TRUE,
           ellipse = T, #不显示置信区间背景色
           size =2, #样本的点大小
           alpha=0.3, #设置点为半透明，出现叠加的效果
           poly= F,
           #polylntyp="dashed",
           #xlim= c(-50,50),
           #ylim= c(-30,30),
           #如果用自定义的颜色，就运行下面这行
           cols = mycol[c(3,2)],
           #cols = mycol[1:length(unique(meta_df$group))],
           arrow = NULL,txt = NULL) + #不画箭头和箭头上的文字
  theme(panel.grid =element_blank()) #+ #去除网格线
  
  #用yyplot添加置信区间圆圈
  # geom_ord_ellipse(ellipse_pro = .95, #先画个.95的圆圈
  #                  color='darkgrey', #圈圈的颜色
  #                  size=0.5, lty=1 ) + #画成虚线，可以用1-6的数字设置为其他线型
  # geom_ord_ellipse(ellipse_pro = .98, #再画个.98的圆圈
  #                  #color='grey', #把这行注释掉，就是跟点一样的颜色
  #                  size=0.5, lty=1 ) 
p
#p+scale_color_nejm()  #新英格兰风格
#p+scale_color_lancet() #柳叶刀风格
###在分组后边加上样本数
table(meta_df$group)
meta_df$group <- ifelse(meta_df$group =="cancer","cancer(471)","normal(41)")
p <- ggord(pca.results, grp_in = meta_df$group, repel=TRUE,
           ellipse = T, #不显示置信区间背景色
           size =2, #样本的点大小
           alpha=0.3, #设置点为半透明，出现叠加的效果
           poly= F,
           #polylntyp="dashed",
           #xlim= c(-50,50),
           #ylim= c(-30,30),
           #如果用自定义的颜色，就运行下面这行
           cols = mycol[c(3,2)],
           #cols = mycol[1:length(unique(meta_df$group))],
           arrow = NULL,txt = NULL) + #不画箭头和箭头上的文字
  theme(panel.grid =element_blank()) #+ #去除网格线

#用yyplot添加置信区间圆圈
# geom_ord_ellipse(ellipse_pro = .95, #先画个.95的圆圈
#                  color='darkgrey', #圈圈的颜色
#                  size=0.5, lty=1 ) + #画成虚线，可以用1-6的数字设置为其他线型
# geom_ord_ellipse(ellipse_pro = .98, #再画个.98的圆圈
#                  #color='grey', #把这行注释掉，就是跟点一样的颜色
#                  size=0.5, lty=1 ) 
p
#保存到pdf文件
ggsave("output/PCA.pdf", width = 6, height = 4)
ggsave("output/PCA.tiff", width = 6, height = 4)
