# devtools::install_github("jokergoo/ComplexHeatmap")
# BiocManager::install("ComplexHeatmap")
# library(devtools)
# devtools::install_github("sqjin/CellChat")


#加载需要的R包
rm(list = ls())
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(Seurat)
sam.name <- "cellchat"
if(!dir.exists(sam.name)){
  dir.create(sam.name)
}
##1.加载数据
load(file = "./multi/multi_experiment.rename.Rdata")
scRNA <- experiment.rename

#CellChat需要两个输入
#一个是细胞的基因表达数据，
#另一个是细胞标签

#获取表达矩阵
cellchat <- createCellChat(scRNA@assays$RNA@data)
#获取标签
scRNA$cellType <- Idents(scRNA)
meta <- data.frame(cellType = scRNA$cellType, row.names =  Cells(scRNA))
#把metadata信息加到CellChat对象中，添加细胞标签
cellchat <- addMeta(cellchat, meta = meta, meta.name = "cellType")
#把细胞标签设置成默认的ID
cellchat <- setIdent(cellchat, ident.use = "cellType") 
#统计每个细胞亚群中的细胞数目
groupSize <- as.numeric(table(cellchat@idents)) 


#CellChat提供了人和小鼠的配受体数据库，分别可以用CellChatDB.human,CellChatDB.mouse来导入
CellChatDB <- CellChatDB.human
View(CellChatDB)

#2.择特定的信息描述细胞间的相互作用，比用一个大的配体库更精细
#用Secreted Signaling来分析细胞通信
#查看有哪些可以选择的子数据库
unique(CellChatDB$interaction$annotation)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")


#将使用的数据库信息写入到cellchat对象中
cellchat@DB <- CellChatDB.use

#抽取信号通路相关基因的表达矩阵
cellchat <- subsetData(cellchat) 


#对表达数据进行预处理，用于细胞间的通信分析。
#首先在一个细胞组中识别过表达的配体或受体，然后将基因表达数据投射到蛋白-蛋白相互作用(PPI)网络上。
#如果配体或受体过表达，则识别过表达配体和受体之间的相互作用。
#cellchat@data.signaling
cellchat <- identifyOverExpressedGenes(cellchat) 
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)


#相互作用推断
#为每个相互作用分配一个概率值并进行置换检验来推断生物意义上的细胞-细胞通信
cellchat <- computeCommunProb(cellchat)

#通过计算与每个信号通路相关的所有配体-受体相互作用的通信概率来推断信号通路水平上的通信概率
#通过计算链路的数量或汇总通信概率来计算细胞间的聚合通信网络
cellchat <- computeCommunProbPathway(cellchat)
df.net <- subsetCommunication(cellchat,slot.name = "netP")
write.csv(df.net,paste0("./",sam.name,"/net_pathway.csv"))

##计算整合的细胞通信网络
cellchat <- aggregateNet(cellchat)
#我们还可以可视化整合的细胞通信网络。例如，使用圆图显示任意两个细胞组之间的相互作用次数或总交互强度（比重）。
groupSize <- as.numeric(table(cellchat@idents))

pdf(paste0("./",sam.name,"/netVisual_circle.pdf"),width = 12,height = 6)

netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength")
dev.off()
#由于细胞通信网络复杂，我们可以检查每个细胞组发送的信号。在这里，我们还控制参数edge.weight.max，以便我们可以比较不同网络之间的边缘权重。

pdf(paste0("./",sam.name,"/netVisual_circle_seperate.pdf"),width = 15,height = 15)
mat <- cellchat@net$weight
par(mfrow = c(4,5), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()



#3.显示重要通信的信号路径
cellchat@netP$pathways
levels(cellchat@idents) 
#在层次绘图的时候，第一列显示的细胞类型的数目
vertex.receiver = seq(1,4) 
#显示的信号通路
pathways.show <- "IL1"
#绘制层次图

pdf(paste0("./",sam.name,"/IL16_signaling_pathway_hierarchy.pdf"),width = 12,height = 6)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "hierarchy", vertex.receiver = vertex.receiver, vertex.weight  = groupSize)  
dev.off()

#绘制环状图
pdf(paste0("./",sam.name,"/CXCL_signaling_pathway_circle.pdf"),width = 10,height = 8)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle", vertex.receiver = vertex.receiver, vertex.weight  = groupSize)  
dev.off()

#绘制弦图
pdf(paste0("./",sam.name,"/IL16_signaling_pathway_chord.pdf"),width = 14,height = 10)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord", vertex.receiver = vertex.receiver, vertex.weight  = groupSize)  
dev.off()

pdf(paste0("./",sam.name,"/netVisual_heatmap.pdf"),width = 14,height = 10)
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
dev.off()

#4.识别细胞群的信号转导作用，通过计算每个细胞群的网络中心性指标，
#CellChat允许随时识别细胞间通信网络中的主要发送者、接收者、调解者和影响者。
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

pdf(paste0("./",sam.name,"/IL16_signalingRole_network.pdf"),width = 7,height = 5)
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
dev.off()

#在 2D 空间中可视化占主导地位的发送器（源）和接收器（目标）
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("IL16"))
pdf(paste0("./",sam.name,"/netAnalysis_signalingRole_scatter.pdf"),width = 7,height = 5)
gg1 + gg2
dev.off()

#识别对某些细胞组的传出或传入信号贡献最大的信号
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
pdf(paste0("./",sam.name,"/netAnalysis_signalingRole_heatmap.pdf"),width = 10,height = 8)
ht1 + ht2
dev.off()
#显示所有的显著的（L-R对）
#从'sources.use'定义的细胞亚群到'targets.use'定义的细胞亚群
pdf(paste0("./",sam.name,"/netVisual_bubble.pdf"),width = 7,height = 8)
netVisual_bubble(cellchat, sources.use = 1:5, targets.use = c(1:19), remove.isolate = FALSE)
dev.off()

#显示具体通路的显著的（L-R对）
pdf(paste0("./",sam.name,"/IL16_netVisual_bubble.pdf"),width = 7,height = 8)
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), signaling = c("IL16","MIF"), remove.isolate = FALSE)
dev.off()

#使用小提琴/点图绘制信号基因表达分布
pdf(paste0("./",sam.name,"/plotGeneExpression_CXCL.pdf"),width = 7,height = 8)
plotGeneExpression(cellchat, signaling = "CXCL")
dev.off()

#5.识别和可视化分泌细胞的传出通信模式
##传出模式揭示了发送者细胞（即作为信号源的细胞）如何相互协调，以及它们如何与某些信号通路协调以驱动通信。
library(NMF)
library(ggalluvial)
selectK(cellchat, pattern = "outgoing")
#当传出模式数为 2 时，Cophenetic 和Silhouette值都开始突然下降。
nPatterns = 2 
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
netAnalysis_river(cellchat, pattern = "outgoing")
netAnalysis_dot(cellchat, pattern = "outgoing")

#识别和可视化目标细胞的传入通信模式
selectK(cellchat, pattern = "incoming")
nPatterns = 6 #当传入模式的数量为 4 时，Cophenetic 值开始下降。
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
netAnalysis_river(cellchat, pattern = "incoming")
netAnalysis_dot(cellchat, pattern = "incoming")
##保存结果
save(cellchat, file = paste0("./",sam.name,"/cellchat.Rdata"))

     