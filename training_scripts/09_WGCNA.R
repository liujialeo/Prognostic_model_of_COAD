rm(list = ls())
library(tidyr)
library(dplyr)
library(WGCNA)
load(file = "./output/COAD_allDiff.Rdata")
diffgene <- allDiff %>% 
        filter(gene !="") %>% 
        filter(adj.P.Val < 0.01) %>% 
        filter(abs(logFC) >2)%>%
        filter(gene_type =="protein_coding" )
load("output/COAD_normalized_counts.Rdata")
normalized_counts <- log(normalized_counts+0.1)
range(normalized_counts)
##将行名中点去掉
term=as.character(row.names(normalized_counts)[1]) 
term =  unlist(strsplit(term, split=".", fixed=T))[1]
rownames(normalized_counts) <- as.character(rownames(normalized_counts))
for (i in 1:nrow(normalized_counts)) {
        print(i)
        term =rownames(normalized_counts)[i]
        term= unlist(strsplit(term, split=".", fixed=T))[1]
        rownames(normalized_counts)[i] =term
}
expro <- cbind(gene_id=rownames(normalized_counts),normalized_counts)
diffgene <- diffgene[,1:2]
expro <- inner_join(diffgene, expro, by = "gene_id")
expro <- expro%>%
        mutate(newcolumn = rowMeans(.[,-c(1:2)])) %>% 
        arrange(desc(newcolumn)) %>% 
        distinct(gene,.keep_all = T) 

rownames(expro) <- expro$gene
expro <- expro[,-(1:2)]

dim(expro)
datExpr=as.data.frame(t(expro));
datExpr <- cbind(sample=rownames(datExpr),datExpr)

for (i in 1:nrow(datExpr)) {
        print(i)
        datExpr$sample[i] <- as.character(substring(datExpr$sample[i],1,16)) 
}
test <- datExpr[1:10,1:10]
datExpr <- datExpr%>%
        mutate(newcolumn = rowMeans(.[,-1])) %>% 
        arrange(desc(newcolumn)) %>% 
        distinct(sample,.keep_all = T) 
rownames(datExpr) <- datExpr$sample
datExpr <- datExpr[,-1]

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

##样本聚类检查离群值##
gsg = goodSamplesGenes(datExpr, verbose = 3);
gsg$allOK
sampleTree = hclust(dist(datExpr), method = "average")
pdf("Sample cluster.pdf",height=6,width=8)
plot(sampleTree, main = "Sample clustering to detect outliers"
     , sub="", xlab="")
dev.off()

##没有离群值，不作处理##

##软阈值筛选##
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
pdf("Soft thresholding power.pdf",height=5,width=10)
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.9,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
#可以查看一下系统推荐的阈值，但是可能不准确
power = sft$powerEstimate 
power

##一步法网络构建：One-step network construction and module detection
net = blockwiseModules(datExpr, power = 5, maxBlockSize = 6000,
                       TOMType = "unsigned", minModuleSize = 30, deepSplit = 2,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "FPKM-TOM", loadTOMs = TRUE,
                       verbose = 3)
table(net$colors)

# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
pdf("Dendrogram.pdf",height=6,width=8)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

##结果保存
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
table(moduleColors)
write.table(moduleLabels,file="moduleLabels.txt",quote=F,sep="\t",col.names=F)
write.table(moduleColors,file="moduleColors.txt",quote=F,sep="\t",col.names=F)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];



##表型与模块相关性##
load(file = "./output/COAD_clinical.Rdata")
rownames(clinical) <- clinical[,1]
clinical <- clinical[,-1]
clinical <- clinical[rownames(datExpr),]
dataclinial  <-  clinical
write.table(dataclinial,"./output/dataclinial.txt",sep="\t",quote=F)
##手动更改后重新读入
dataclinical <- data.table::fread(file = "./output/dataclinical.csv",data.table = F) #保存csv格式时第一列有信息丢失
colnames(dataclinical)[1] <- "TCGA_id"
rt=read.table("output/diffgene_risk.txt",header=T,sep="\t") ##读入多因素COX分析结果，将riskScore作为一种临床信息
colnames(rt)[1] <- "TCGA_id"
rt <- rt[,c(1,8)]

dataclinical <- inner_join(rt,dataclinical, by = "TCGA_id")

rownames(dataclinical) <- dataclinical[,1]
dataclinical <- dataclinical[,-1]


MEs0 = moduleEigengenes(datExpr,moduleColors)$eigengenes
MEs0 <- MEs0[rownames(dataclinical),]
MEs = orderMEs(MEs0)
modul_clinical_cor = cor(MEs, dataclinical, use = "p")
write.table(modul_clinical_cor,"output/module-clinial-cor.txt",sep="\t",quote=F)
modul_clinical_p = corPvalueStudent(modul_clinical_cor, nSamples)
write.table(modul_clinical_p,"output/modul-clinical-p.txt",sep="\t",quote=F)
textMatrix = paste(signif(modul_clinical_cor, 2), " (", signif(modul_clinical_p, 1), ")",sep = "")
dim(textMatrix) = dim(modul_clinical_cor)
pdf("output/Module&Clinical.pdf",height=8,width=10)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = modul_clinical_cor, xLabels = names(dataclinical), 
               yLabels = names(MEs),ySymbols = names(MEs), colorLabels = FALSE, 
               colors = blueWhiteRed(50), textMatrix = textMatrix, setStdMargins = FALSE, 
               cex.text = 1, zlim = c(-1,1), main = paste("Module-trait relationships"))
dev.off()
### 如何把基因数目显示出来
dd1 <- data.frame(ME=names(MEs),color=substring(names(MEs),3))
dd2 <- data.frame(table(moduleColors))
colnames(dd2) <- c("color","num")
## 合并
dd3 <- merge(dd1,dd2,by="color")
rownames(dd3) <- dd3$ME
dd3 <- dd3[names(MEs),]

### c(bottom, left, top, right)
pdf("output/Module&Clinical_genenumber.pdf",height=8,width=12)
par(mar = c(6, 10, 3, 3))
ynumbers <- paste0(names(MEs),paste0("(",dd3$num,")"))
ynumbers
labeledHeatmap(Matrix = modul_clinical_cor, xLabels = names(dataclinical), 
               yLabels = names(MEs),ySymbols = ynumbers, colorLabels = FALSE, 
               colors = blueWhiteRed(50), textMatrix = textMatrix, setStdMargins = FALSE, 
               cex.text = 1, zlim = c(-1,1), main = paste("Module-trait relationships"))
dev.off()
###选择感兴趣的模块进行分析 GS-MM
# Firstly, calculate the correlation matrix between module and genes.
# names (colors) of the modules
modNames = substring(names(MEs), 3)
datExpr <- datExpr[rownames(MEs),]
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
# Calculate the Pearson correlation coefficient matrix of each module and its genes
# MEs is the value of each module in each sample
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")


# Secondly, calculate the correlation matrix between conditions and genes. 
# Only continuous properties can be computed. If the variables are discrete, the matrix is converted to 0-1 when the sample table is constructed.
# Here, the variable whether or not it belongs to the P condition is numeralized with 0 and 1.
#选择哪一个临床信息，此处以第7列risk_score为例
P = as.data.frame(dataclinical[,1]) # choose an interested condition!!
names(P) = "riskScore"
geneTraitSignificance = as.data.frame(cor(datExpr, P, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(P), sep="")
names(GSPvalue) = paste("p.GS.", names(P), sep="")

# Then, combine aboved two correlation matrixes
module = "turquoise" # choose interested module
column = match(module, modNames)
# get the genes in the interested module
moduleGenes = moduleColors==module
sizeGrWindow(7, 7)
par(mfrow = c(1,1))
# Genes that are highly correlated with conditions are also highly associated with modules
pdf("output/GS-MM (turquoise-risk_score).pdf")
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = "Gene significance for risk_score",
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module) 
dev.off()

#绘制两两模块间的邻接矩阵
pdf("output/Eigengene adjacency heatmap.pdf",height=8,width=10)
plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap",plotDendrograms = T, xLabelsAngle= 90)
dev.off()

###导出基因模块到Cytoscape###
# Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(datExpr, power = 4);
# Select modules
modules = c("turquoise");
# Select module probes
probes = names(datExpr)
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("output/CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("output/CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               #altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule]);


## 可视化基因网络TOM
# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 4);
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7;
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;
# Call the plot function
pdf("output/TOM.pdf")
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot")
dev.off()
