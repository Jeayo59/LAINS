#!/data/sft/software/R-3.6.1/bin/Rscript
Args <- commandArgs(T)
library(Seurat)
library(dplyr)
library(magrittr)
library(ggplot2)
LUAD <- readRDS(Args[1])
npcs <- 30
pcSelect <- as.numeric(Args[2])
resolution <- as.numeric(Args[3])

###################################TSNE聚类分析和marker基因###################################
##TSNE聚类分析
LUAD <- FindNeighbors(object = LUAD, dims = 1:pcSelect)                #计算邻接距离
LUAD <- FindClusters(object = LUAD, resolution = resolution)                  #对细胞分组,优化标准模块化
LUAD <- RunTSNE(object = LUAD, dims = 1:pcSelect)                      #TSNE聚类
write.csv(LUAD$seurat_clusters,file="01.tsneCluster.csv",quote=F)
embed_tsne <- Embeddings(LUAD,'tsne')
write.csv(embed_tsne,"02.embed_tsne.csv",quote=F)

## group by cluster
plot1 = DimPlot(LUAD,reduction="tsne",label=T,raster = F)
ggsave("03.tsne.jpeg",plot=plot1,width=8,height=7)

plot2 = DimPlot(LUAD,reduction="tsne",group.by='orig.ident',raster = F)
ggsave("04.tsne_sample.jpeg",plot=plot2,width=8,height=7)

##UMAP
LUAD <- RunUMAP(object = LUAD, dims = 1:pcSelect)

embed_umap <- Embeddings(LUAD,'umap')
write.csv(embed_umap,"05.embed_umap.csv",quote=F)

plot3 = DimPlot(LUAD,reduction="umap",label=T,raster = F)
ggsave("06.umap.jpeg",plot=plot3,width=8,height=7)

plot4 = DimPlot(LUAD,reduction="umap",group.by='orig.ident',raster = F)
ggsave("07.umap_sample.jpeg",plot=plot4,width=8,height=7)

##寻找差异表达的特征
#DefaultAssay(LUAD) <- "RNA"
logFCfilter=0.5
adjPvalFilter=0.05
LUAD.markers <- FindAllMarkers(object = LUAD,
                               only.pos = FALSE,
                               min.pct = 0.25,
                               logfc.threshold = logFCfilter)
sig.markers=LUAD.markers[(abs(as.numeric(as.vector(LUAD.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(LUAD.markers$p_val_adj))<adjPvalFilter),]
write.table(sig.markers,file="08.markers.csv",sep=",",row.names=F,quote=F)

top10 <- LUAD.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
#绘制marker在各个cluster的热图
#pdf(file="09.tsneHeatmap.pdf",width=12,height=9)
#DoHeatmap(object = LUAD, features = top10$gene) + NoLegend()
#dev.off()
DefaultAssay(LUAD) <- "integrated"
saveRDS(LUAD,file="LUAD_second.rds")
