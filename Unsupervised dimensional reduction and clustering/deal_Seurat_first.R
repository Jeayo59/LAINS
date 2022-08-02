Args <- commandArgs(T)
library(Seurat)
library(dplyr)
npcs <- 30
LUAD <- readRDS(Args[1])
############# Normalization ########################################################################
#LUAD <- NormalizeData(object = LUAD, normalization.method = "LogNormalize", scale.factor = 10000)

######### HVG ######################################################################################
#LUAD <- FindVariableFeatures(object = LUAD, selection.method = "mvp", mean.function = ExpMean, dispersion.function = LogVMR,mean.cutoff = c(0.0125,3), dispersion.cutoff = c(0.5,Inf))
#LUAD <- FindVariableFeatures(object = LUAD,selection.method = "vst", nfeatures = 2000)
###### 
#top10 <- head(x = VariableFeatures(object = LUAD), 10)
#pdf(file="01.featureVar.pdf",width=10,height=6)              #保存基因特征方差图
#plot1 <- VariableFeaturePlot(object = LUAD)
#plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
#CombinePlots(plots = list(plot1, plot2))
#dev.off()
## re integrated ##
DefaultAssay(LUAD) <- "RNA"
LUAD1 <- subset(LUAD,subset=(orig.ident=="LUAD1"))
LUAD2 <- subset(LUAD,subset=(orig.ident=="LUAD2"))

LUADlist = c(LUAD1,LUAD2)
for(i in 1:length(LUADlist)){
    LUADlist[[i]]<- NormalizeData(LUADlist[[i]])
    LUADlist[[i]] <- FindVariableFeatures(LUADlist[[i]], selection.method = "vst")
}
LUAD.anchors <- FindIntegrationAnchors(object.list = LUADlist)
LUAD <- IntegrateData(anchorset = LUAD.anchors)

############# PCA ###################################################################################
LUAD=ScaleData(LUAD)        
LUAD=RunPCA(object= LUAD,npcs = npcs,pc.genes=VariableFeatures(object = LUAD))

#绘制每个PCA成分的相关基因
pdf(file="02.pcaGene.pdf",width=20,height=16)
VizDimLoadings(object = LUAD, dims = 1:npcs, reduction = "pca",nfeatures = 20)
dev.off()

#主成分分析图形
pdf(file="03.PCA.pdf",width=6.5,height=6)
DimPlot(object = LUAD, reduction = "pca")
dev.off()

#主成分分析热图
pdf(file="04.pcaHeatmap.pdf",width=30,height=24)
DimHeatmap(object = LUAD, dims = 1:npcs, cells = 1000, balanced = TRUE,nfeatures = 30,ncol=2)
dev.off()

#每个PC的p值分布和均匀分布
#LUAD <- JackStraw(object = LUAD, dim = npcs, num.replicate = 100) ## 1h20min
#LUAD <- ScoreJackStraw(object = LUAD, dims = 1:npcs)
#pdf(file="05.pcaJackStraw.pdf",width=8,height=6)
#JackStrawPlot(object = LUAD, dims = 1:npcs)
#dev.off()

pdf(file="06.elbowPlot.pdf",width=8,height=6)
ElbowPlot(object = LUAD, ndims = npcs)
dev.off()
print(LUAD)
saveRDS(LUAD,file="LUAD.rds")
