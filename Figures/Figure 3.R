library(Seurat)
library(dplyr)
library(ggplot2)
library(scales)
library(ggsci)
library(gridExtra)
library(stringi)
library(magrittr)
library(ggsignif)
## Figure 3
setwd("I:/LUAD_merge/manuscript_fiture/Figure3/")
Myeloid <- readRDS("I:/LUAD_merge/02.Seurat/remove_effect_re_integrate/09.Myeloid/LUAD_second.rds")
anno <- read.table("I:/LUAD_merge/02.Seurat/remove_effect_re_integrate/09.Myeloid/Myeloid_anno.txt",header=FALSE)
head(anno)
colnames(anno) <- c("Cell","subtype","Maintype","Patient_id","Type","SmokingType","TS")
rownames(anno) <- anno$Cell
anno <- anno[rownames(Myeloid@meta.data),]
anno$subtype <- sub("Alveolar_macrophage","TRM",anno$subtype)
anno$subtype <- sub("Macrophage","MDM",anno$subtype)
anno$subtype <- sub("Proliferating_macrophages1","FABP4hi_pro",anno$subtype)
anno$subtype <- sub("Proliferating_macrophages2","SPP1hi_pro",anno$subtype)
Myeloid$subtype <- anno$subtype
Myeloid$TS <- anno$TS
Myeloid$Patient_id <- anno$Patient_id
unique(anno$subtype)
Myeloid$subtype1 <- Myeloid$subtype
Myeloid$subtype1 <- sub("FABP4hi_pro","Proliferating mac",Myeloid$subtype1)
Myeloid$subtype1 <- sub("SPP1hi_pro","Proliferating mac",Myeloid$subtype1)

DimPlot(Myeloid,group.by = "subtype1",pt.size=0.2)+
  scale_color_manual(values=c(mypal[c(5,2,8)],"#807DBA","#54278F",mypal[c(1,3,6,9,10)]),
                     limits = c("TRM","MDM","Proliferating mac",
                                "CD14+_monocyte","CD16+_monocyte",
                                "Activated_DCs","cDC1s","cDC2s","pDC",
                                "Neutrophils"))+
  guides(colour = guide_legend(override.aes = list(size=2.5)))+
  theme_classic()+
  theme(legend.position = "bottom",
        plot.title = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.key.size = unit(2,"pt"))
ggsave("Myeloid_subtype.pdf",height=5,width=4.5)

#########markers
Myeloid.marker <- c("MARCO","FABP4", # TRM
                    "C1QA","CD68","APOE","SPP1","CCL2", #MDM
                    "TUBB","BIRC5","CDKN2A","TYMS","MKI67","STMN1", #Proliferating mac
                    "CD14","VCAN","FCN1", # CD14+ monocytes
                    "FCGR3A","CDKN1C", # CD16+ monocytes
                    "CCR7","LAMP3", # Activated DC
                    "CLEC9A","XCR1", # cDC1s
                    "CLEC10A","CD1C", # cDC2s
                    "CLEC4C","PTCRA", # pDC
                    "CSF3R","FCGR3B","ALPL","CXCR2","CMTM2") # Neutrophils

Myeloid.marker.mean <- AverageExpression(Myeloid,features = Myeloid.marker,assays = "RNA",group.by="subtype1")$RNA
pheatmap(Myeloid.marker.mean[,c(10,6,9,2,3,1,4,5,8,7)],scale="row",cluster_rows = F,cluster_cols = F,
         color = c(colorRampPalette(c("darkslategray","ghostwhite"))(100),colorRampPalette(c("ghostwhite","chocolate4"))(100)),
         legend = T,filename = "signature.pdf",height = 5,width=4.5)

##########trajectory Macrophage
library(Seurat)
library(monocle)
LUAD <- readRDS("Macrophage.rds")
cds <- readRDS("monocle.rds")
trajectory_XY <- as.data.frame(t(as.matrix(reducedDimS(dat))))
trajectory_XY$State <- dat$State
head(trajectory_XY)

df_trajactory <- trajectory_XY
head(df_trajactory)
df_trajactory$Cell <- rownames(df_trajactory)
colnames(df_trajactory) <- c("X","Y","State","Cell")
metadata <- read.table("metadata_Myeloid.txt",sep="\t",header=TRUE)
head(metadata)
anno <- read.table("all_anno.txt",sep="\t",header=TRUE)
head(anno)
colnames(metadata)
metadata <- metadata[,-19]
metadata <- merge(metadata,anno,by="Cell",all.x=TRUE)
mac <- metadata[metadata$subtype=="Alveolar_macrophage"|metadata$subtype=="Macrophage"|metadata$subtype=="Proliferating_macrophages1"|metadata$subtype=="Proliferating_macrophages2",]
df_trajectory <- merge(mac,
                       df_trajactory[,c("Cell","X","Y","State")],
                       by="Cell",all.x=TRUE)
head(df_trajectory)

ggplot(df_trajectory,aes(x=X,y=Y,col = subtype))+
  geom_point(size=0.3)+
  labs(x="Component 1",y="Component 2",color="Origin type")+
  scale_color_manual(values = c(mypal[c(5,2)],"orange","chocolate4"))+
  guides(colour = guide_legend(override.aes = list(size=4)))+
  theme_classic()
ggsave("cell_trajectory_by_subtype.pdf",height=5,width=7)

##############DEG
library(ggrepel)
sig <- read.csv("Proliferating_macrophages1_vs_2.csv")
head(sig)
colnames(sig)[1] <- "gene"
sig$gene[sig$change=="Stable"] <- ""
TNM$gene[rownames(TNM)=="IFI27"] <- "IFI27"
TNM$gene[rownames(TNM)=="MDK"] <- "MDK"
sig$change = ifelse(sig$p_val_adj < 0.01 & abs(sig$avg_log2FC) >= 0.25, 
                    ifelse(sig$avg_log2FC> 0.25 ,'Up','Down'),
                    'Stable')
p <- ggplot(data = sig, 
            aes(x = avg_log2FC, 
                y = abs(pct.1-pct.2), 
                colour=change,
                label = gene)) +
  geom_point(alpha=0.4, size=3) +
  scale_color_manual(values=c("chocolate4", "grey","orange"))+
  #xlim(c(-4.5, 4.5)) +
  geom_vline(xintercept=c(-0.25,0.25),lty=4,col="black",lwd=0.8) +
  #geom_hline(yintercept = -log10(0.000001),lty=4,col="black",lwd=0.8) +
  labs(x="avg_log2FC",
       y="abs(pct.1-pct.2)",
       title="")  +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank(),
        axis.text = element_text(size=18),
        axis.text.x = element_text(size=18),
        legend.text=element_text(size=18),
        axis.title = element_text(size=18)
  )
p
p+geom_text_repel(data = sig, aes(x = avg_log2FC, 
                                  y = abs(pct.1-pct.2), 
                                  label = gene),
                  size = 3,box.padding = unit(0.5, "lines"),
                  point.padding = unit(0.8, "lines"), 
                  segment.color = "black", 
                  show.legend = FALSE)
ggsave("sig.pdf",height=10,width=10)

#######signature
Gene <- c("TYMS","MKI67","STMN1","SPP1","MARCO","INHBA","FABP4")
MAC <- LUAD[,LUAD@meta.data$subtype_new=="TRM"|LUAD@meta.data$subtype_new=="MDM"|LUAD@meta.data$subtype_new=="FABP4hi_pro"|LUAD@meta.data$subtype_new=="SPP1hi_pro"]
DotPlot(object =MAC,features =Gene,group.by="subtype_new",dot.scale=8)+
  scale_y_discrete(limits = c("TRM","MDM","FABP4hi_pro","SPP1hi_pro"))+
  coord_flip()+
  scale_colour_gradient2(low="lightgray",high="chocolate4",mid="ghostwhite")+
  labs(x="",y="")+
  theme(axis.text = element_text(size=8))+
  theme(axis.text.x = element_text(angle=30,hjust=1,vjust=1,size=8))+
  theme(legend.text=element_text(size=8),legend.title=element_text(size=8))
ggsave("signature_gene_dotplot.pdf",height=4,width=4)


####proportion
my.df <- read.csv("I:/LUAD_merge/09.tissue_preference/Cell_statistics_updata.csv",header=T,row.names = 1,check.names = FALSE)
dat <- data.frame(Sample = my.df$Sample,
                  PM2per = my.df$PM2/my.df$Total,
                  MDMper = my.df$Macrophage/my.df$Total) %>%
  left_join(unique(LUAD_marker@meta.data[,c("Sample","TS","Type")]),by="Sample")

judgeDegree <- function(value){
  if(value >= 0.05){
    type = "NS."
  }
  else if(value <0.05 & value >= 0.01){
    type = "*"
  }
  else if(value < 0.01 & value >= 0.001){
    type = "**"
  }
  else if(value <= 0.001){
    type = "***"
  }
  return(type)
}
anno <- function(y,height){
  
  annotation_table <- data.frame('xsl' = c(1,3),'xel' = c(1,3), 'ysl' = c(y,y+height),'yel' = c(y+y/30,y+y/30+height),
                                 'xsr' = c(2,4),'xer' = c(2,4), 'ysr' = c(y,y+height),'yer' = c(y+y/30,y+y/30+height),
                                 'xsh' = c(1,3),'xeh' = c(2,4), 'ysh' = c(y+y/30,y+y/30+height),'yeh' = c(y+y/30,y+y/30+height),
                                 'xsu' = c(1.5,3.5),'xeu' = c(1.5,3.5), 'ysu' = c(y+y/30,y+y/30+height),'yeu' = c(y+y/15+abs(height),y+y/15+abs(height)))
  annotation_table2 <- data.frame('xsuh' = 1.5,'xeuh' = 3.5, 'ysuh' =y+y/15+abs(height),'yeuh' = y+y/15+abs(height))
  
  result <- list(annotation_table,annotation_table2)
  return(result)
}
dat <- reshape2::melt(dat,id.vars = c("Sample","TS","Type"))
head(dat)
colnames(dat)[4] <- "subtype"
compare <- function(subtype,h){
  # subtype = "MDMper"
  # h=0.45
  df <- dat[dat$subtype==subtype,]
  p.nNSvsnS <- wilcox.test(df[df$TS=="nNS",]$value,df[df$TS=="nS",]$value)$p.value %>% judgeDegree()
  p.nNSvstNS <- wilcox.test(df[df$TS=="nNS",]$value,df[df$TS=="tNS",]$value)$p.value %>% judgeDegree()
  p.nSvstS <- wilcox.test(df[df$TS=="nS",]$value,df[df$TS=="tS",]$value)$p.value %>% judgeDegree()
  p.tNSvstS <- wilcox.test(df[df$TS=="tNS",]$value,df[df$TS=="tS",]$value)$p.value %>% judgeDegree()
  p.TvsN <- wilcox.test(df[df$Type=="Tumor",]$value,df[df$Type=="Normal",]$value)$p.value %>% judgeDegree()
  if(p.TvsN != "NS."){
    p =  wilcox.test(df[df$Type=="Tumor",]$value,df[df$Type=="Normal",]$value,alternative="greater")$p.value
    print(p)
    if(p<0.1){height=h/10}
    else{height=-1*h/10}
  }
  else{height=0}
  
  annotation_table <- anno(h,height)[[1]]
  annotation_table1 <- anno(h,height)[[2]]
  loc <- data.frame(x=c(2,3,4,4,2.5),
                    y=c(h-h/25,h-h/25+height,h-h/25+height,h+(h/20+height),h+(h/10+abs(height))),
                    value=c(p.nNSvsnS,p.nNSvstNS,p.nSvstS,p.tNSvstS,p.TvsN),
                    type=c(1,2,2,3,4))
  
  ggplot(df,aes(x=TS,y=value,color=TS))+
    # geom_violin(size=0.8)+
    geom_boxplot(outlier.size = 0.3,size=0.5)+
    geom_segment(inherit.aes = F,data= annotation_table,aes(x=xsl,xend=xel,y=ysl,yend=yel),color="grey40")+
    geom_segment(inherit.aes = F,data= annotation_table,aes(x=xsr,xend=xer,y=ysr,yend=yer),color="grey40")+
    geom_segment(inherit.aes = F,data= annotation_table,aes(x=xsh,xend=xeh,y=ysh,yend=yeh),color="grey40")+
    geom_segment(inherit.aes = F,data= annotation_table,aes(x=xsu,xend=xeu,y=ysu,yend=yeu),color="grey40")+
    geom_segment(inherit.aes = F,data= annotation_table1,aes(x=xsuh,xend=xeuh,y=ysuh,yend=yeuh),color="grey40")+
    geom_text(data=loc,aes(x=x,y=y,label=value,color=as.factor(type)))+
    scale_color_manual(values = c("black","orange","red","brown","black","black","black","black"))+
    scale_y_continuous(labels = scales::percent)+
    labs(x="",y="Percentage")+
    theme_classic()+
    theme(plot.title = element_text(size=14,hjust=0.5,vjust=0.5),
          axis.text.x = element_text(size=13),
          axis.text.y = element_text(size=12),
          axis.title.y= element_text(size=14),
          legend.position = "none")
}

compare("MDMper",0.5)
ggsave("MDM_percentage.pdf",width=2.5,height=3.5)
compare("PM2per",0.12)
ggsave("SPP1hipro_percentage_new.pdf",width=2.5,height=3.5) 

##########trajectory MDM
library(monocle)
setwd("I:/LUAD_merge/06.monocle/MDM_pro2_mono/")
cds <- readRDS("monocle.rds") 
ordergene <- rownames(cds)
meta <- pData(cds)
head(meta)
unique(meta$subtype)
plot_cell_trajectory(cds,color_by = "subtype",cell_size=0.2,alpha=0.8,show_branch_points = F)+
  guides(colour = guide_legend(override.aes = list(size=2),ncol=2))+
  labs(color="")+
  scale_color_manual(values=c("#807DBA","#54278F",mypal[2],"chocolate4"),
                     limits = c("CD14+_monocyte","CD16+_monocyte","Macrophage","Proliferating_macrophages2"),
                     labels = c("CD14+ mono","CD16+ mono","MDM","SPPhi pro"))+
  theme_bw()+
  theme(legend.position=c(0.7,0.2))+
  theme(legend.key.size = unit(0.8,"pt"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        strip.text = element_blank())
ggsave("trajectory_MDM.pdf",height=3.5,width=3.5)


