#!/data/sft/software/R-3.6.1/bin/Rscript
Args <- commandArgs()
library(Seurat)
library(dplyr)
library(magrittr)
library(ggplot2)
LUAD <- readRDS(Args[6])
npcs <- 30
pcSelect <- as.numeric(Args[7])
resolution <- as.numeric(Args[8])

###################################TSNE聚类分析和marker基因###################################
##TSNE聚类分析
LUAD <- FindNeighbors(object = LUAD, dims = 1:pcSelect)                #计算邻接距离
LUAD <- FindClusters(object = LUAD, resolution = resolution)                  #对细胞分组,优化标准模块化
write.table(LUAD$seurat_clusters,file="01.tsneCluster.txt",sep = "\t",quote=F,col.names=F)

