#!/data/user/zengz/Software/anaconda2/envs/py3/anaconda3/envs/r-devel/bin/R
Args <- commandArgs()
library(Seurat)
library(infercnv)
library(plyr)
library(dplyr)

cutoff=Args[7]
outdir=Args[8]

LUAD <- readRDS(Args[6])
DefaultAssay(LUAD) <- "RNA"
counts <- GetAssayData(object = LUAD,slot = "data")
counts <- as.data.frame(counts)
ref_group_names = c("normal")

infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = counts,
                                    annotations_file = "./cell_annotation.txt",
                                    delim = "\t",
                                    gene_order_file = "./human.gene.positions1.txt",
                                    ref_group_names = ref_group_names)

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=as.numeric(cutoff),  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir=outdir,  # 输出文件夹
                             cluster_by_groups=F,   # 聚类
                             denoise=T, #去噪
                             HMM=F)
