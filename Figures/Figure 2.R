library(Seurat)
library(dplyr)
library(ggplot2)
library(scales)
library(ggsci)
library(gridExtra)
library(stringi)
library(magrittr)
library(ggsignif)
## Figure 2
setwd("I:/LUAD_merge/manuscript_fiture/Figure2/")
Epi_n.cell <- readRDS("I:/UAD_merge/02.Seurat/remove_effect_re_integrate/03.Epithelial_normal/LUAD_second.rds")
Epi_n.metadata <- Epi_n.cell@meta.data
Epi_n.metadata$Cell <- rownames(Epi_n.metadata)
Epi_n.umap <- read.csv("I:/LUAD_merge/02.Seurat/remove_effect_re_integrate/03.Epithelial_normal/05.embed_umap.csv")
colnames(Epi_n.umap)[1] <- c("Cell")

Epi_n.metadata_umap <- merge(Epi_n.metadata,Epi_n.umap,by="Cell",all.x=TRUE)
head(Epi_n.metadata_umap)
Epi_n.metadata_umap$Patient_id <- substr(Epi_n.metadata_umap$Sample,1,3)
Epi_n.metadata_umap$subtype[Epi_n.metadata_umap$seurat_clusters==0|
                              Epi_n.metadata_umap$seurat_clusters==2|
                              Epi_n.metadata_umap$seurat_clusters==6] <- "Ciliated"
Epi_n.metadata_umap$subtype[Epi_n.metadata_umap$seurat_clusters==3|
                              Epi_n.metadata_umap$seurat_clusters==9] <- "AT2"
Epi_n.metadata_umap$subtype[Epi_n.metadata_umap$seurat_clusters==4] <- "AT1"
Epi_n.metadata_umap$subtype[Epi_n.metadata_umap$seurat_clusters==1|
                              Epi_n.metadata_umap$seurat_clusters==7] <- "Club"

Epi_n.metadata_umap$subtype[Epi_n.metadata_umap$seurat_clusters==5|
                              Epi_n.metadata_umap$seurat_clusters==8] <- "Lymphocyte infiltrate"

Epi_n.metadata_umap$subtype[Epi_n.metadata_umap$seurat_clusters==10] <- "Basal"

ggplot(Epi_n.metadata_umap,aes(x=UMAP_1,y=UMAP_2,color=subtype))+
  geom_point(size=0.5)+
  labs(color="subtype")+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  scale_color_manual(values = mypal[2:7])+
  theme_classic()+
  theme(legend.position="bottom",
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())
ggsave("Subtype_normal.pdf",height=4,width=4)


Epi.cell <- readRDS("I:/LUAD_merge/02.Seurat/remove_effect_re_integrate/08.Epithelial_all/LUAD_add_info.rds")
Epi.metadata <- Epi.cell@meta.data
Epi.metadata$Cell <- rownames(Epi.metadata)
Epi.metadata$TS <- tidyr::unite(Epi.metadata,"TS",Type,SmokingType)$TS
Epi.metadata$TS <- sub("Normal_","n",Epi.metadata$TS)
Epi.metadata$TS <- sub("Tumor_","t",Epi.metadata$TS)

head(Epi.metadata)
Epi.cell <- AddMetaData(Epi.cell,
                        metadata = Epi.metadata[rownames(Epi.cell@meta.data),"TS"],
                        col.name ="TS")
head(Epi.cell@meta.data)
umap <- read.csv("I:/LUAD_merge/02.Seurat/remove_effect_re_integrate/08.Epithelial_all/05.embed_umap.csv",header=T)
colnames(umap)[1] <- "Cell"
Epi.metadata <- merge(Epi.metadata,umap,by="Cell",all.x=TRUE)
ggplot(Epi.metadata,aes(x=UMAP_1,y=UMAP_2,color=subtype))+
  geom_point(size=0.5)+
  labs(color="subtype")+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  scale_color_manual(values = c(mypal[2:7],mypal[1],"grey"))+
  theme_classic()+
  theme(legend.position="bottom",
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())
ggsave("Subtype.pdf",height=5,width=5)

## split
ggplot(Epi.metadata,aes(x=UMAP_1,y=UMAP_2,color=subtype))+
  geom_point(size=0.5)+
  facet_grid(Type~SmokingType)+
  labs(color="subtype")+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  scale_color_manual(values = c(mypal[2:7],mypal[1],"grey"))+
  theme_classic()+
  theme(legend.position="none",
        axis.title = element_blank(),
        strip.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())
ggsave("subtype_by_Smoking_by_Type.pdf",height=4.5,width=5)

## AT2 signature
AT2_sig <- c("LAMP1","FOXA2","NAPSA","SFTA3","NKX2-1","ABCA3","SCGB3A1","SFTPB","SFTPC","SFTPD")
AT2_sig_lst <- list(AT2_sig)

Epi_sub <- subset(Epi.cell,subset = (subtype=="AT2"|subtype=="malignant"))
Epi_sub <- AddModuleScore(Epi_sub,features = AT2_sig_lst,name="AT2_sig")

AT2_sig_df <- Epi_sub@meta.data %>% select(c("TS","AT2_sig1"))
comparisons <- list(c("tNS","tS"),c("nNS","tNS"),c("nS","tS"),c("tNS","tS"),c("nNS","nS"))
ggplot(AT2_sig_df,aes(x=TS,y=AT2_sig1,fill=TS))+
  geom_boxplot(outlier.size = 0.4,width=0.6)+
  labs(x="",y="AT2 signature score")+
  scale_fill_manual(values = c("darkseagreen3","darkseagreen4","salmon3","salmon4"))+
  theme_bw()+
  theme(legend.position = "none")
ggsave("AT2_signature1.pdf",width=4,height=3)


## signature
tNS.signature <- read.table("I:/LUAD_merge/05.signature/Epithelial_tumor_all_normal/tNS_signature.txt",header=FALSE)$V1 %>%
  as.character()
tNS.signature[5] <- "GPR116"
tS.signature <- read.table("I:/LUAD_merge/05.signature/Epithelial_tumor_all_normal/tS_signature.txt",header=FALSE)$V1 %>%
  as.character()

Epi.cell.fi <- subset(Epi.cell,subset=(subtype!="non_malignant"))
TS <- tidyr::unite(Epi.cell.fi@meta.data,"TS",Type,SmokingType)
TS$TS <- sub("Normal_","n",TS$TS)
TS$TS <- sub("Tumor_","t",TS$TS)

Epi.cell.fi <- AddMetaData(Epi.cell.fi,metadata=TS[rownames(Epi.cell.fi@meta.data),"TS"],col.name = "TS")

mean_exp <- AverageExpression(Epi.cell.fi,assays = "RNA",features=c(tNS.signature,tS.signature),group.by = "TS")$RNA %>%
  as.data.frame()
library(pheatmap)
pheatmap(mean_exp[,c(3,4,1,2)],scale="row",cluster_rows = T,cluster_cols = F,
         treeheight_row = 0,
         color = c(colorRampPalette(c("darkslategray","ghostwhite"))(100),colorRampPalette(c("ghostwhite","chocolate4"))(100)),
         legend = T,filename = "signature.pdf",height = 6,width=4)

## Prognosis
library(UCSCXenaTools)
library(survival)
library(survminer)
suppressMessages(library(dplyr))
luad_cohort = XenaData %>% 
  filter(XenaHostNames == "tcgaHub") %>%
  XenaScan("TCGA Lung Adenocarcinoma")
luad_cohort

cli_query = luad_cohort %>%
  filter(DataSubtype == "phenotype") %>%
  XenaGenerate() %>%
  XenaQuery() %>% 
  XenaDownload()
cli_query

cli = XenaPrepare(cli_query)
cli
clinicalMatrix <- as.data.frame(cli$LUAD_clinicalMatrix)
cli.s = cli$LUAD_survival.txt
ge = luad_cohort %>%
  filter(DataSubtype == "gene expression RNAseq", Label == "IlluminaHiSeq")
## tS
tS = fetch_dense_values(host = ge$XenaHosts, 
                        dataset = ge$XenaDatasets, 
                        identifiers = tS.signature,
                        use_probeMap = FALSE)
tS_mean <- colMeans(tS)

tS_merged_data = tibble(sample = names(tS_mean),
                        tS_mean_expression = as.numeric(tS_mean)) %>% 
  left_join(cli.s, by = "sample") %>%
  filter(OS.time <=1825) %>%
  filter(substr(sample, 14, 15) == "01") %>%   # Keep only 'Primary Tumor'
  # filter(sample %in% TCGA_S_sample) %>%
  select(sample, tS_mean_expression, OS.time, OS) %>% 
  rename(time = OS.time, 
         status = OS)

fit = coxph(Surv(time/30, status) ~ tS_mean_expression, data = tS_merged_data)
fit
summary(fit)$coefficients[,5]
tS_merged_data = tS_merged_data %>% 
  mutate(group = case_when(
    tS_mean_expression > quantile(tS_mean_expression, 0.5) ~ 'tS_High',
    tS_mean_expression < quantile(tS_mean_expression, 0.5) ~ 'tS_Low',
    TRUE ~ NA_character_
  ))
fit = survfit(Surv(time/30, status) ~ group, data = tS_merged_data)
ggsurvplot(fit, pval = TRUE,
           ggtheme = theme_bw(),
           risk.table.col = "strata",
           palette = c("#E7B800", "#2E9FDF"),
           size=0.5)
ggsave("tS_signature_KM.pdf",height=3,width=3.5)

##############cancer cell
setwd("F:/01.scRNA-seq/01.LUAD_merge/14.cancer_cell/")
LUAD <- readRDS("LUAD_second.rds")
head(LUAD@meta.data)
table(LUAD@meta.data$Sample)
DimPlot(LUAD,reduction="umap",group.by='Sample',raster = F)
table(LUAD@meta.data$seurat_clusters,LUAD@meta.data$Sample)
ggsave("umap_sample.pdf",height=7,width=9)
DimPlot(LUAD,reduction="umap",raster = F,label = T)
ggsave("umap.pdf",height=7,width=8)

anno <- read.table("F:/01.scRNA-seq/01.LUAD_merge/All_add_info.txt",sep="\t")
head(anno)
anno <- anno[,c(1,7)]
head(anno)
colnames(anno) <- c("Cell","TS")
metadata <- LUAD@meta.data
head(metadata)
metadata <- merge(metadata,anno,by="Cell",all.x=T)
rownames(metadata) <- metadata$Cell
metadata <- metadata[rownames(LUAD@meta.data),]
LUAD@meta.data <- metadata
DimPlot(LUAD,reduction="umap",group.by='TS',raster = F)
ggsave("umap_TS.pdf",height=7,width=8)


LUAD  <- AddModuleScore(LUAD,features = list(MHCII),name = "MHCII")
head(LUAD@meta.data)
ggplot(LUAD@meta.data,aes(x=seurat_clusters,y=MHCII1,fill= seurat_clusters))+
  geom_boxplot()+
  theme_bw()+
  labs(x="",y="",title="MHC-II")+
  theme(legend.position="none",
        plot.title = element_text(hjust=0.5))
ggsave("MHC-II_signature.jpg",height=3,width= 7)

