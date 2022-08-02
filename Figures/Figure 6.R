library(pheatmap)
library(psych)
library(qgraph)
library(igraph)
library(tidyverse)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(grid)
library(Cairo)
############Figure 6
mynet <- read.delim("I:/LUAD_merge/08.cellphone_main/All_again/count_network.txt", check.names = FALSE,sep='\t')
length(mynet[,1])
mynet <- mynet[!grepl("undefined",mynet$SOURCE) & !grepl("undefined",mynet$TARGET),]
mynet <- mynet[mynet$SOURCE!=mynet$TARGET,]
unique(mynet$SOURCE)
cellType <- c("Epithelial","Fibrobast","Endothelial",
              "CD8_GZMK+","CD8_GZMB+","CD8_Exhausted","CD4_Active","CD4_Naive","NK_CD16+","B",
              "cDC1s","cDC2s","Alveolar_macrophage","Macrophage","Proliferating_macrophages1","Proliferating_macrophages2",
              "nTreg","CD14+_monocyte","CD16+_monocyte")

head(mynet)
mynet1 <- data.frame(matrix(ncol=3,nrow=0))
colnames(mynet1) <- colnames(mynet)
mynet1
for(cell in cellType){
  sub <- mynet[mynet$SOURCE==cell & mynet$TARGET %in% cellType,]
  mynet1 = rbind(mynet1,sub)}
mynet1$SOURCE <- sub("Proliferating_macrophages2","SPPhi pro",mynet1$SOURCE)
mynet1$SOURCE <- sub("Proliferating_macrophages1","FABP4hi pro",mynet1$SOURCE)
mynet1$SOURCE <- sub("Alveolar_macrophage","TRM",mynet1$SOURCE)
mynet1$SOURCE <- sub("Macrophage","MDM",mynet1$SOURCE)
mynet1$SOURCE <- sub("Fibrobast","Fibroblast",mynet1$SOURCE)

mynet1$TARGET <- sub("Proliferating_macrophages2","SPPhi pro",mynet1$TARGET)
mynet1$TARGET <- sub("Proliferating_macrophages1","FABP4hi pro",mynet1$TARGET)
mynet1$TARGET <- sub("Alveolar_macrophage","TRM",mynet1$TARGET)
mynet1$TARGET <- sub("Macrophage","MDM",mynet1$TARGET)
mynet1$TARGET <- sub("Fibrobast","Fibroblast",mynet1$TARGET)
mynet2 <- mynet1[mynet1$count>20,]
net<- graph_from_data_frame(mynet2)
net


net
E(net)$width  <- E(net)$count/25

karate_groups <- cluster_optimal(net)
coords <- layout_in_circle(net, order =order(membership(karate_groups)))
coords

allcolour=c("#20B2AA","#FFA500","#9370DB",
            "#FA8072","#9400D3",
            "#800080","#A0522D","#D2B48C","#D2691E","#87CEEB",
            "#5F9EA0",
            "#FFE4B5","#228B22","#E9967A","#4682B4",
            "#F0E68C","#EE82EE","#FF6347",
            "#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
pdf("All_interaction.pdf",height=8,width=8)
plot(net, edge.arrow.size=0, 
     edge.curved=0, 
     vertex.color=allcolour,
     vertex.frame.color=allcolour,
     vertex.label.color="black",
     layout = coords,
     vertex.label.cex=.7,
     vertex.label.dist=2)
dev.off()

## Immune check point
setwd("I:/LUAD_merge/14.checkpoint/")
library(Seurat)
library(ggplot2)
library(ggridges)
LUAD_checkpoint <- readRDS("LUAD_checkpoint.rds")

anno <- read.table("All_add_info.txt",header=F)
colnames(anno) <- c("Cell","subtype","mainType","Patient","Type","SmokingType","TS")
anno$subtype[anno$subtype=="tS1" |
               anno$subtype=="tS2"|
               anno$subtype=="tS3"] <- "malignant"
anno$subtype[anno$subtype=="Macrophage"] <- "MDM"
anno$subtype[anno$subtype=="Alveolar_macrophage"] <- "TRM"
anno$subtype[anno$subtype=="Proliferating_macrophages1"] <- "FABP4hi_pro"
anno$subtype[anno$subtype=="Proliferating_macrophages2"] <- "SPP1hi_pro"

rownames(anno) <- anno$Cell
head(LUAD_checkpoint@meta.data)
LUAD_checkpoint <- AddMetaData(LUAD_checkpoint,metadata = anno[rownames(LUAD_checkpoint@meta.data),
                                                               c("subtype","Patient","TS")],
                               col.name = c("subtype","Patient","TS"))
gene <- rownames(LUAD_checkpoint)
gene
unique(anno$subtype)

cellType_select <- read.csv("I:/LUAD_merge/08.cellphone_malignant_stage/All/subtype1_select_1.txt",header = F,sep="\t")
colnames(cellType_select) <- c("subtype","mainType")
cellType_select$subtype[cellType_select$subtype=="tS1"] <- "malignant"
cellType_select$subtype[cellType_select$subtype=="tS2"] <- "malignant"
cellType_select$subtype[cellType_select$subtype=="tS3"] <- "malignant"
cellType_select <- cellType_select[cellType_select$subtype!="iTreg",]

LUAD_checkpoint_sub <- subset(LUAD_checkpoint,subset=c(subtype=="malignant"|
                                                         subtype == "Tumor_ECs"|
                                                         subtype == "Lipofibroblast"|
                                                         subtype == "MMP-high"|
                                                         subtype == "Myofibroblast"|
                                                         subtype == "CD14+_monocyte"|
                                                         subtype == "CD16+_monocyte"|
                                                         subtype == "cDC1s"|
                                                         subtype == "cDC2s"|
                                                         subtype == "FABP4hi_pro"|
                                                         subtype == "TRM"|
                                                         subtype == "SPP1hi_pro"|
                                                         subtype == "MDM"|
                                                         subtype == "NK_CD16+"|
                                                         subtype == "nTreg"|
                                                         subtype =="CD4_Active"|
                                                         subtype == "CD4_Naive"|
                                                         subtype == "CD8_Exhausted"|
                                                         subtype == "CD8_GZMK+"|
                                                         subtype =="CD8_GZMB+"))
DotPlot(LUAD_checkpoint_sub,feature=gene,group.by="subtype")+
  scale_colour_gradient2(low="darkslategray",high="chocolate4",mid="ghostwhite")+
  labs(x="",y="")+
  scale_y_discrete(limits=c("malignant","Tumor_ECs","Lipofibroblast","MMP-high","Myofibroblast",
                            "CD14+_monocyte","CD16+_monocyte","cDC1s","cDC2s",
                            "FABP4hi_pro","TRM","SPP1hi_pro","MDM",
                            "NK_CD16+","nTreg","CD4_Active","CD4_Naive",
                            "CD8_Exhausted","CD8_GZMK+","CD8_GZMB+"))+
  theme(axis.text.x = element_text(angle=30,hjust=1,vjust=1,size=15))
ggsave("LUAD_checkpoint_dotplot.pdf",height=6,width=14)

TMB <- data.frame(Patient = c("P05","P07","P09","P11","P17","P18","P21","P22",
                              "P01","P02","P03","P04","P06","P08","P13","P14","P15","P16","P19","P20"),
                  TMB = c(4.8,3.84,5.76,1.92,5.76,11.52,1.92,6.72,
                          1.92,2.88,2.88,3.84,6.72,4.8,3.84,1.92,3.84,2.88,4.8,1.92),
                  TS = c(rep("S",8),rep("NS",12)))
TMB
library(ggplot2)
library(ggsignif)
ggplot(TMB,aes(x=TS,y=TMB,color=TS))+
  # geom_violin()+
  geom_boxplot(outlier.size = 0.5,width=0.5)+
  scale_color_manual(values = c("salmon3","salmon4"))+
  geom_signif(comparisons = list(c("NS","S")),step_increase = 0.05,map_signif_level = F,
              test = wilcox.test,size=0.2,textsize = 4,tip_length=0.01,color="grey44")+
  scale_color_manual(values = c("salmon3","salmon4")) +
  scale_x_discrete(limits = c("NS","S"),labels = c("Non-smoking","Smoking"))+
  theme_classic()+
  labs(x="",y="TMB (Muts/Mb)")+
  theme(legend.position="none")

ggsave("TMB.pdf",height=3.2,width=3)
t.test(TMB[TMB$TS=="NS",]$TMB,TMB[TMB$TS=="S",]$TMB,alternative = "less")

## Immune related gene
dat <- read.csv("LUAD_Immune_related_score.csv",check.names = FALSE,header=T,row.names = 1)
dat$Cell <- rownames(dat)
dim(dat)
dim(df)
dat <- left_join(dat,df[,c("Cell","subtype","TS")],by="Cell",all.x=TRUE)
head(dat)
## PDL1
ggplot(dat[dat$subtype=="malignant",],aes(x=TS,y=PDL1,color=TS))+
  geom_violin(size=1)+
  geom_boxplot(width=0.5,outlier.size = 0.2)+
  geom_signif(comparisons = list(c("tNS","tS")),step_increase = 0.05,map_signif_level = F,
              test = t.test,size=0.2,textsize = 4,tip_length=0.01,color="grey44")+
  scale_color_manual(values = c("salmon3","salmon4"))+
  theme_classic()+
  labs(x="",y="PDL1 score")+
  theme(legend.position="none")
ggsave("PDL1_in_malignant.pdf",height=3.2,width=3)
ggsave("PDL1_in_malignant.jpg",height=3.2,width=3)

## CD47
ggplot(dat[dat$subtype=="malignant" & (dat$TS=="tNS"|dat$TS=="tS"),],aes(x=TS,y=CD47,color=TS))+
  geom_violin(size=1)+
  geom_boxplot(width=0.2,outlier.size = 0.2,size=1)+
  geom_signif(comparisons = list(c("tNS","tS")),step_increase = 0.05,map_signif_level = F,
              test = t.test,size=0.2,textsize = 4,tip_length=0.01,color="grey44")+
  scale_color_manual(values = c("salmon3","salmon4"))+
  theme_classic()+
  labs(x="",y="CD47 score")+
  theme(legend.position="none")
ggsave("CD47_in_malignant.pdf",height=3.2,width=3)
ggsave("CD47_in_malignant.jpg",height=3.2,width=3)
