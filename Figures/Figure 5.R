library(Seurat)
library(dplyr)
library(ggplot2)
library(scales)
library(ggsci)
library(gridExtra)
library(stringi)
library(magrittr)
library(ggsignif)
## Figure 5
setwd("I:/LUAD_merge/manuscript_fiture/Figure5/")
T.NK.cell <- readRDS("I:/LUAD_merge/02.Seurat/remove_effect_re_integrate/04.T_NK/LUAD_second.rds")
metadata <- T.NK.cell@meta.data
metadata$Cell <- rownames(metadata)
head(metadata)
DefaultAssay(T.NK.cell) <- "RNA"

umap <- read.csv("I:/LUAD_merge/02.Seurat/remove_effect_re_integrate/04.T_NK/05.embed_umap.csv",row.names = 1)
umap$Cell <- rownames(umap)
head(umap)
metadata_umap <- merge(metadata,umap,by="Cell",all.x=TRUE)
anno <- subset(metadata,select=c(SmokingType,Type))
anno <- tidyr::unite(anno,"TS",Type,SmokingType)
anno$TS <- sub("Normal_","n",anno$TS)
anno$TS <- sub("Tumor_","t",anno$TS)
anno$Cell <- rownames(anno)
metadata_umap <- merge(metadata_umap,anno,by="Cell",all.x=TRUE)
metadata_umap$Patient_id <- substr(metadata_umap$Sample,1,3)
rownames(metadata_umap) <- metadata_umap$Cell
head(metadata_umap)

df <- read.csv("I:/LUAD_merge/02.Seurat/remove_effect_re_integrate/04.T_NK/LUAD_T.NK_metadata.csv",row.names=1)
T.NK.cell$subtype <- df[rownames(metadata_umap),]$subtype
unique(T.NK.cell$subtype)
T.NK.cell$subtype[T.NK.cell$subtype=="iTreg"] <- "Undefined"
T.NK.cell$subtype[T.NK.cell$subtype=="nTreg"] <- "Treg"
DimPlot(T.NK.cell,group.by = "subtype",pt.size=0.2)+
  scale_color_manual(values=c(mypal,"grey"),
                     limits = c("CD4 Active","CD4 Naive",
                                "CD8 GZMB+","CD8 GZMK+","CD8 Exhausted","CD8 undefined",
                                "NK CD16+","NK CD16-",
                                "Treg","gamma-delta T cell",
                                "Undefined"))+
  guides(colour = guide_legend(override.aes = list(size=2.5)))+
  theme_classic()+
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.key.size = unit(2,"pt"))
ggsave("TNK_subtype_2.pdf",height=5,width=5.5)
dev.off()

T_marker <- c("CD3D","CD3E","CD3G")
NK_marker <- c("XCL1","FCGR3A","KLRD1","KLRF1")
CD8_marker <- c("CD8A","CD8B")
CD4_marker <- c("CD4","IL7R")
Naive_marker <- c("TCF7","SELL","LEF1","CCR7")
active_CD4_marker <- c("RBPJ","LGALS3","CD40LG")
Cytotocxic_marker <- c("IL2","GZMA","GNLY","PRF1","GZMB","IFNG","NKG7","FGFBP2","CST7","GZMH","GZMK")
Exhausted_marker <- c("LAG3","TIGIT","PDCD1","HAVCR2","CTLA4")
Treg_marker <- c("IL2RA","FOXP3","IKZF2")
gamma_delta_marker <- c("STMN1","TUBB","MKI67")
marker <- c(T_marker,NK_marker,CD4_marker,CD8_marker,Naive_marker,
            active_CD4_marker,Cytotocxic_marker,
            Exhausted_marker,
            Treg_marker,gamma_delta_marker)

DefaultAssay(T.NK.cell) <- "RNA"


## marker_gene mean expression
marker <- intersect(marker,rownames(T.NK.cell))
data_marker.m <- AverageExpression(T.NK.cell,features = marker,slot = "data")$RNA
data_marker.m
library(pheatmap)
cluster <-c("2","4","14",
            "0","9","5","15",
            "1","7","12","3","8","6",
            "10",
            "13",
            "11","16")

data_marker.m <- data_marker.m[marker,cluster]
pheatmap(data_marker.m,scale = "row",cluster_rows = F,cluster_cols = F,
         # color = colorRampPalette(colors = c("mediumblue","white","brown3"))(100),
         gaps_col = c(3,7,13,14,15),
         gaps_row = c(3,7,9,11,15,18,29,34,37),
         color = c(colorRampPalette(c("darkslategray","ghostwhite"))(100),colorRampPalette(c("ghostwhite","chocolate4"))(100)),
         filename = "T_NK_marker.pdf",height=6.5,width=4.5)

metadata_umap$TS
df_type <- as.data.frame(table(metadata_umap$seurat_clusters,metadata_umap$TS))

ggplot(df_type,aes(x=Var1,y=Freq,fill=Var2))+
  geom_col(position = "fill")+
  scale_y_continuous(expand=c(0,0),labels = scales::percent)+
  scale_x_discrete(limits = cluster)+
  scale_fill_manual(values = c("darkseagreen3","darkseagreen4","salmon3","salmon4"))+
  labs(x="",y="ratio %",fill="")+
  theme_bw()+
  theme(legend.key.size = unit(8,"pt"))

ggsave("Cluster_sample_types_percentage.pdf",height=3,width = 8,units="cm")

dat <- data.frame(Sample = my.df$Sample,
                  Tregper = (my.df$nTreg)/(my.df$CD4_Active+my.df$CD4_Naive+my.df$nTreg)
) %>%
  left_join(unique(metadata_umap[,c("Sample","TS","Type")]),by="Sample")

dat <- reshape2::melt(dat,id.vars = c("Sample","TS","Type"))
head(dat)
colnames(dat)[4] <- "subtype"
compare("Tregper",0.9)
ggsave("Treg_percentage.pdf",height=3,width = 3.2)

## CD8 GZMB
## Non Smoking
CD8_GZMB.NS.signature <- read.csv("I:/LUAD_merge/05.signature/T_cell/T_CD8_GZMB/Non_smoking_Tumor_vs_Normal_signature.csv")
colnames(CD8_GZMB.NS.signature)[1] <- "Gene"

CD8_GZMB.NS.signature$sig[CD8_GZMB.NS.signature$p_val_adj<0.01] <- "Yes"
CD8_GZMB.NS.signature$sig[CD8_GZMB.NS.signature$p_val_adj>=0.01] <- "No"
CD8_GZMB.NS.signature$y <- abs(CD8_GZMB.NS.signature$pct.1 - CD8_GZMB.NS.signature$pct.2)
CD8_GZMB.NS.signature$Type[CD8_GZMB.NS.signature$avg_log2FC>0 &
                             CD8_GZMB.NS.signature$sig=="Yes"] <- "Up"
CD8_GZMB.NS.signature$Type[CD8_GZMB.NS.signature$avg_log2FC<0 &
                             CD8_GZMB.NS.signature$sig=="Yes"] <- "Down"

CD8_GZMB.NS.signature$Type[CD8_GZMB.NS.signature$sig=="No"] <- "not sig."

ggplot(CD8_GZMB.NS.signature,aes(x=avg_log2FC,y=y,color=Type))+
  geom_point(size=0.5)+
  theme_classic()+
  geom_vline(xintercept = c(-0.25,0.25),linetype="dashed")+
  scale_color_manual(values = c("darkseagreen3","grey","salmon3"))+
  theme(legend.key.size = unit(8,"pt"))+
  guides(colour = guide_legend(override.aes = list(size=2)))+
  labs(x="avg_log2FC",y="abs(pct.1-pct.2)")+
  geom_text_repel(data =CD8_GZMB.NS.signature[CD8_GZMB.NS.signature$sig=="Yes",],
                  aes(avg_log2FC, y, label = Gene),size=2.5,max.overlaps = 15,verbose = FALSE)
ggsave("NS_DE.pdf",height=5,width=5)

head(CD8_GZMB.NS.signature)
CD8_GZMB.NS.signature[CD8_GZMB.NS.signature$gene=="GZMA",]

## Smoking
CD8_GZMB.S.signature <- read.csv("I:/LUAD_merge/05.signature/T_cell/T_CD8_GZMB/Smoking_Tumor_vs_Normal_signature.csv")
colnames(CD8_GZMB.S.signature)[1] <- "Gene"
CD8_GZMB.S.signature$sig[CD8_GZMB.S.signature$p_val_adj<0.01] <- "Yes"
CD8_GZMB.S.signature$sig[CD8_GZMB.S.signature$p_val_adj>=0.01] <- "No"
CD8_GZMB.S.signature$y <- abs(CD8_GZMB.S.signature$pct.1 - CD8_GZMB.S.signature$pct.2)
CD8_GZMB.S.signature$Type[CD8_GZMB.S.signature$avg_log2FC>0 &
                            CD8_GZMB.S.signature$sig=="Yes"] <- "Up"
CD8_GZMB.S.signature$Type[CD8_GZMB.S.signature$avg_log2FC<0 &
                            CD8_GZMB.S.signature$sig=="Yes"] <- "Down"

CD8_GZMB.S.signature$Type[CD8_GZMB.S.signature$sig=="No"] <- "not sig."
unique(CD8_GZMB.S.signature$Type)
ggplot(CD8_GZMB.S.signature,aes(x=avg_log2FC,y=y,color=Type))+
  geom_point(size=0.5)+
  theme_classic()+
  geom_vline(xintercept = c(-0.25,0.25),linetype="dashed")+
  scale_color_manual(values = c("darkseagreen4","grey","salmon4"))+
  theme(legend.key.size = unit(8,"pt"))+
  guides(colour = guide_legend(override.aes = list(size=2)))+
  geom_text_repel(data =CD8_GZMB.S.signature[CD8_GZMB.S.signature$sig=="Yes",],
                  aes(avg_log2FC, y, label = Gene),size=2.5,
                  max.overlaps = 15)+
  labs(x="avg_log2FC",y="abs(pct.1-pct.2)")
ggsave("S_DE.pdf",height=5,width=5)

gene <- c("GZMB","GZMA","GZMH","GZMM","GZMK","NKG7","PRF1","KLRF1","GNLY","FGFBP2",
          "TNF","KLRG1","IL7R")
head(T.NK.cell@meta.data)
head(metadata_umap)
T.NK.cell$TS <- metadata_umap[rownames(T.NK.cell@meta.data),]$TS
GZMB <- subset(T.NK.cell,subset=(subtype=="CD8 GZMB+"))
DotPlot(CD8_GZMB,feature=gene,group.by="TS")+
  coord_flip()+
  scale_colour_gradient2(low="darkslategray",high="chocolate4",mid="ghostwhite")+
  scale_x_discrete(limits = rev(gene))+
  labs(x="",y="")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=30,hjust=1,vjust=1,size=12))
ggsave("GZMB_cytotoxic.pdf",height=3,width=4)