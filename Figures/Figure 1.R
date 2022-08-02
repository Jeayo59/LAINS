library(Seurat)
library(dplyr)
library(ggplot2)
library(scales)
library(ggsci)
library(gridExtra)
library(stringi)
library(magrittr)
library(ggsignif)

mypal <- pal_npg("nrc",alpha=0.7)(10)
show_col(mypal)
## Figure 1
## cluster
setwd("I:/LUAD_merge/manuscript_fiture/Figure1/")
LUAD_marker <- readRDS("LUAD_marker.rds")
DefaultAssay(LUAD_marker) <- "RNA"
metadata <- LUAD_marker@meta.data
metadata$Cell <- rownames(metadata)
metadata$Patient_id <- sapply(metadata$Sample,function(x) substr(x,1,3))
Sample_info <- read.table("I:/LUAD_merge/Sample_info.txt",sep="\t",header=T)
Sample_info$SmokingType <- sub("Never-smoker","NS",Sample_info$SmokingType)
Sample_info$SmokingType <- sub("Smoker","S",Sample_info$SmokingType)
head(Sample_info)
metadata <- left_join(metadata,Sample_info[,c("Patient_id","SmokingType")],by="Patient_id",
                      all.x=TRUE)
metadata$TS <- tidyr::unite(metadata,"TS",Type,SmokingType)$TS
metadata$TS <- sub("Normal_","n",metadata$TS)
metadata$TS <- sub("Tumor_","t",metadata$TS)
rownames(metadata) <- metadata$Cell

LUAD_marker <- AddMetaData(LUAD_marker,
                           metadata = metadata[rownames(LUAD_marker@meta.data),"TS"],
                           col.name ="TS")
head(LUAD_marker@meta.data)
## main cell type
cell_color <- data.frame(Cell_type=unique(metadata[order(metadata$Cell_type,decreasing = F),]$Cell_type),
                         color = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A4", "#BDB76B"))

UMAPPlot(LUAD_marker,group="Cell_type",pt.size=0.1,raster=FALSE)+
  theme(legend.position = "none")+
  scale_color_manual(values = as.character(cell_color$color))+
  # geom_text(data=XY,aes(x=x,y=y,label=Cell_type),size=8)+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        plot.title = element_blank())
ggsave("LUAD_main_cell_type.pdf",height=5,width=5.5)
ggsave("LUAD_main_cell_type.jpg",height=5,width=5.5)

## split
UMAPPlot(LUAD_marker,group.by="Cell_type",split.by="TS",raster=FALSE)+
  scale_color_manual(values = as.character(cell_color$color))+
  # geom_text(data=XY,aes(x=x,y=y,label=Cell_type),size=8)+
  theme_bw()+
  theme(legend.position = "none")+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        plot.title = element_blank())
ggsave("LUAD_main_cell_type_split.pdf",height=5,width=22)
ggsave("LUAD_main_cell_type_split.jpg",height=5,width=22)

head(metadata)
# marker gene
marker <- c("KRT19","KRT18","EPCAM","CDH1",# Epithelial
            "CD3D","CD3E","CD3G","CD2", # T
            "GZMA","NKG7","KLRD1","GNLY", # NK
            "CD79A","MS4A1", # B
            "MARCO","LYZ","CD68", # Myeloid
            "KIT","GATA2","TPSAB1","CPA3","HPGDS","RGS13",# MAST
            "THY1","DCN","COL1A2","COL1A1", # Fibroblast
            "RAMP2","FLT1","CLDN5" # Endothelial
) 


DotPlot(LUAD_marker,feature=marker,group.by="Cell_type")+
  scale_colour_gradient2(low="darkslategray",high="chocolate4",mid="ghostwhite")+
  labs(x="",y="")+
  scale_y_discrete(limits=c("Endothelial","Fibrobast","MAST","Myeloid","B","NK","T","Epithelial"),
                   labels = c("Endothelial cells","Fibroblasts","MAST cells","Myeloid cells","B lymphocytes","NK cells",
                              "T lymphocytes","Epithelial cells"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle=30,hjust=1,vjust=1,size=10))
ggsave("main_cellType_maker.pdf",height=3,width=10)