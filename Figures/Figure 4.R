library(Seurat)
library(dplyr)
library(ggplot2)
library(scales)
library(ggsci)
library(gridExtra)
library(stringi)
library(magrittr)
library(ggsignif)
## Figure 4
##########################M1/M2
setwd("F:/01.scRNA-seq/01.LUAD_merge/Fig3/M1M2/")
M1 <- c("CCR7","KYNU","IDO1","CD40","IRF1","IRF5","CCL5","IL6","IL1B","IL1A","CD86","CXCL11","CXCL10","CXCL9","TNF")
M2 <- c("IRF4","FN1","MSR1","VTCN1","CD276","TNFSF8","TNFSF12","WNT7B","CLEC7A","MMP9","MMP14","TGFB3","TGFB2","TGFB1","CTSD","CTSB","CTSA","EGF","VEGFC","VEGFB","VEGFA","LYVE1","CCL24","CCL22","CCL18","CCL17","CCL20","CCL13","CCL4")

Mac <- subset(LUAD,subset=(subtype=="MDM" |
                             subtype=="TRM"|
                             subtype=="SPP1hi_pro"|
                             subtype=="FABP4hi_pro"))

Mac <- AddModuleScore(Mac,features = list(M1),name = "M1")
Mac <- AddModuleScore(Mac,features = list(M1),name = "M2")

metadata_M1M2 <- read.table("metadata_mac.txt",header=TRUE)
head(metadata_M1M2)
head(anno2)
colnames(metadata_M1M2) 
metadata_M1M2 <- metadata_M1M2[,-18]
metadata_M1M2$Cell <- rownames(metadata_M1M2)
metadata_M1M2 <- merge(metadata_M1M2,anno2,by="Cell",all.x=TRUE)


metadata <- Mac@meta.data
metadata$Cell <- rownames(metadata)
metadata_M1M2 <- merge(metadata_M1M2,metadata,by="Cell",all.x=T)
ggplot(metadata_M1M2[metadata_M1M2$TS=="tNS",])+
  geom_point(aes(x=M21,y=M11,color=subtype),size=0.5)+
  scale_color_manual(values = c(mypal[c(5,2)],"#DC0101"))+
  theme_classic()+
  theme(legend.position="null")
ggsave("tNS_M1M2_3.pdf",height=4.5,width=3)
ggplot(metadata_M1M2[metadata_M1M2$TS=="tS",])+
  geom_point(aes(x=M21,y=M11,color=subtype),size=0.5)+
  scale_color_manual(values = c(mypal[c(5,2)],"#DC0101"))+
  theme_classic()+
  theme(legend.position="null")
ggsave("tS_M1M2_3.pdf",height=4.5,width=3)

#########score
Mac_metadata <- Mac@meta.data
M1_marker <- c("IL23","TNF","CXCL9","CXCL10","CXCL11","CD86","IL1A","IL1B","IL6",
               "CCL5","IRF5","IRF1","CD40","IDO1","KYNU","CCR7")
M2_marker <- c("L4R","CCL4","CCL13","CCL20","CCL17","CCL18","CCL22","CCL24","LYVE1",
               "VEGFA","VEGFB","VEGFC","VEGFD","EGF","CTSA","CTSB","CTSC","CTSD",
               "TGFB1","TGFB2","TGFB3","MMP14","MMP19","MMP9","CLEC7A","WNT7B",
               "FASL","TNFSF12","TNFSF8","CD276","VTCN1","MSR1","FN1","IRF4")

anti_inflammatory_marker <- c("IL4","IL10","IL13","IFNA1","IFNA2","TGFB1","FCN1",
                              "CD163","CD36","S100A9","S100A8","VCAN","CD14","F13A1",
                              "APOE","SEPP1","C1QB","C1QA","C1QC")

AP_marker <- c("B2M","CIITA","CTSS","HLA-A","HLA-B","HLA-C","HLA-DMA","HLA-DMB",
               "HLA-DOA","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQA2","HLA-DQB1",
               "HLA-DRA","HLA-DRB1","HLA-DRB5","HLA-E","HLA-F","PDIA3","PSME2",
               "TAP1","TAP2","TAPBP")
M1_marker <- intersect(M1_marker,rownames(Mac))
M2_marker <- intersect(M2_marker,rownames(Mac))
anti_inflammatory_marker <- intersect(anti_inflammatory_marker,rownames(Mac))
AP_marker_marker <- intersect(AP_marker,rownames(Mac))

M1_marker_lst <- list(M1_marker)
M2_marker_lst <- list(M2_marker)
anti_inflammatory_marker_lst <- list(anti_inflammatory_marker)
AP_marker_marker_lst <- list(AP_marker_marker)

Mac <- AddModuleScore(Mac,features = M1_marker_lst,name="M1")
Mac <- AddModuleScore(Mac,features = M2_marker_lst,name="M2")
Mac <- AddModuleScore(Mac,features = anti_inflammatory_marker_lst,name="anti")
Mac <- AddModuleScore(Mac,features = AP_marker_marker_lst,name="AP")

mydata <- FetchData(Mac,vars = c("UMAP_1","UMAP_2","M11","M21","Type","SmokingType","subtype","anti1","AP1"))
colnames(mydata) <- c("UMAP_1","UMAP_2","M1","M2","Type","SmokingType","subtype","anti","AP")
head(mydata)
anno <- tidyr::unite(mydata,"TS",Type,SmokingType) %>% select(TS)
head(anno)
head(mydata)
anno$TS <- sub("Normal_","n",anno$TS)
anno$TS <- sub("Tumor_","t",anno$TS)
mydata$Cell <- rownames(mydata)
anno$Cell <- rownames(anno)
mydata <- merge(mydata,anno,by="Cell",all.x=TRUE)
head(mydata)

## compare 
judgeDegree <- function(value){
  if(value >= 0.1){
    type = "NS."
  }
  else if(value <0.1 & value >= 0.05){
    type = "*"
  }
  else if(value < 0.05 & value >= 0.01){
    type = "**"
  }
  else if(value <= 0.01){
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

compare <- function(cellType,name,geneSet,cat,h){
  df <- mydata[mydata$subtype==cellType,]
  p.nNSvsnS <- wilcox.test(df[df$TS=="nNS",][[geneSet]],df[df$TS=="nS",][[geneSet]])$p.value %>% judgeDegree()
  p.nNSvstNS <- wilcox.test(df[df$TS=="nNS",][[geneSet]],df[df$TS=="tNS",][[geneSet]])$p.value %>% judgeDegree()
  p.nSvstS <- wilcox.test(df[df$TS=="nS",][[geneSet]],df[df$TS=="tS",][[geneSet]])$p.value %>% judgeDegree()
  p.tNSvstS <- wilcox.test(df[df$TS=="tNS",][[geneSet]],df[df$TS=="tS",][[geneSet]])$p.value %>% judgeDegree()
  p.TvsN <- wilcox.test(df[df$Type=="Tumor",][[geneSet]],df[df$Type=="Normal",][[geneSet]])$p.value %>% judgeDegree()
  if(p.TvsN != "NS."){
    p =  wilcox.test(df[df$Type=="Tumor",][[geneSet]],df[df$Type=="Normal",][[geneSet]],alternative="greater")$p.value
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
  ggplot(df,aes_string(x=cat,y=geneSet,color=cat))+
    geom_violin(size=0.8)+
    geom_boxplot(width=0.2,outlier.size = 0.5,size=0.8)+
    geom_segment(inherit.aes = F,data= annotation_table,aes(x=xsl,xend=xel,y=ysl,yend=yel),color="grey40")+
    geom_segment(inherit.aes = F,data= annotation_table,aes(x=xsr,xend=xer,y=ysr,yend=yer),color="grey40")+
    geom_segment(inherit.aes = F,data= annotation_table,aes(x=xsh,xend=xeh,y=ysh,yend=yeh),color="grey40")+
    geom_segment(inherit.aes = F,data= annotation_table,aes(x=xsu,xend=xeu,y=ysu,yend=yeu),color="grey40")+
    geom_segment(inherit.aes = F,data= annotation_table1,aes(x=xsuh,xend=xeuh,y=ysuh,yend=yeuh),color="grey40")+
    geom_text(data=loc,aes(x=x,y=y,label=value,color=as.factor(type)))+
    scale_color_manual(values = c("black","orange","red","brown",
                                  "darkseagreen3","darkseagreen4","salmon3","salmon4"))+
    labs(x="",y=paste0(geneSet," score"),title=name)+
    theme_classic()+
    theme(plot.title = element_text(size=14,hjust=0.5,vjust=0.5),
          axis.text.x = element_text(size=13),
          axis.text.y = element_text(size=12),
          axis.title.y= element_text(size=14),
          legend.position = "none")
}

AM_mydata <- mydata[mydata$subtype=="TRM",]
Mac_mydata <- mydata[mydata$subtype=="MDM",]


compare("TRM","TRM","M1","TS",2)
ggsave("AM_M1_score.jpeg",height=4.5,width=3)


compare("MDM","MDM","M1","TS",1.65)
ggsave("Mac_M1_score.jpeg",height=4.5,width=3)


## M2
compare("TRM","TRM","M2","TS",1)
ggsave("AM_M2_score.jpeg",height=4.5,width=3)


compare("MDM","MDM","M2","TS",1.1)
ggsave("Mac_M2_score.jpeg",height=4.5,width=3)

## anti
compare("TRM","TRM","anti","TS",1.6)
ggsave("AM_anti_score.jpeg",height=4.5,width=3)


compare("MDM","MDM","anti","TS",1.7)
ggsave("Mac_anti_score.jpeg",height=4.5,width=3)


## AP
compare("TRM","TRM","AP","TS",2.4)
ggsave("AM_AP_score.jpeg",height=4.5,width=3)


compare("MDM","MDM","AP","TS",2.7)
ggsave("Mac_AP_score.jpeg",height=4.5,width=3)



## SPP1hi Mac
SPP1hi_Mac <- subset(Mac,subset=(seurat_clusters=="2" | seurat_clusters=="10"))
SPP1hi_Mac_tumor <- subset(SPP1hi_Mac,subset=(Type=="Tumor"))
SPP1hi_Mac_tumor <- AddModuleScore(SPP1hi_Mac_tumor,features = c("SPP1"),name="SPP1")

unique(SPP1hi_Mac_tumor@meta.data$seurat_clusters)
df <- FetchData(SPP1hi_Mac_tumor,vars = c("TS","M11","M21","TAM1","Angiogenesis1",
                                          "Phagocytosis1","anti1","pro1","AP1","IFN1"))
colnames(df) <- c("TS","M1","M2","TAM","Angiogenesis",
                  "Phagocytosis","anti","pro","AP","IFN")
head(df)
df.m <- reshape2::melt(df,id.vars = c("TS"))
head(df.m)

ggplot(df,aes(x=TS,y=M1,color=TS))+
  geom_violin(size=0.9)+
  geom_boxplot(width=0.4,outlier.size = 0.4,size=0.9)+
  theme_classic()+
  labs(x="",y="M1 score")+
  geom_signif(comparisons = list(c("tNS","tS")),step_increase = 0.05,map_signif_level = T,
              test = t.test,size=0.2,textsize = 4,tip_length=0.01,color="grey44")+
  scale_color_manual(values = c("salmon3","salmon4"))+
  theme(legend.position = "none",
        axis.text.x = element_text(size=13),
        axis.text.y = element_text(size=12),
        axis.title.y= element_text(size=14))
ggsave("SPP1hi_Mac_M1_tNS_vs_tS.pdf",height=2.5,width=2.5)


ggplot(df,aes(x=TS,y=M2,color=TS))+
  geom_violin(size=0.8)+
  geom_boxplot(width=0.4,outlier.size = 0.4,size=0.8)+
  theme_classic()+
  labs(x="",y="M2 score")+
  geom_signif(comparisons = list(c("tNS","tS")),step_increase = 0.05,map_signif_level = T,
              test = t.test,size=0.2,textsize = 4,tip_length=0.01,color="grey44")+
  scale_color_manual(values = c("salmon3","salmon4"))+
  theme(legend.position = "none",
        axis.text.x = element_text(size=13),
        axis.text.y = element_text(size=12),
        axis.title.y= element_text(size=14))
ggsave("SPP1hi_Mac_M2_tNS_vs_tS.pdf",height=2.5,width=2.5)

ggplot(df,aes(x=TS,y=anti,color=TS))+
  geom_violin(size=0.8)+
  geom_boxplot(width=0.4,outlier.size = 0.4,size=0.8)+
  theme_classic()+
  labs(x="",y="Anti-inflammatory score")+
  geom_signif(comparisons = list(c("tNS","tS")),step_increase = 0.05,map_signif_level = T,
              test = t.test,size=0.2,textsize = 4,tip_length=0.01,color="grey44")+
  scale_color_manual(values = c("salmon3","salmon4"))+
  theme(legend.position = "none",
        axis.text.x = element_text(size=13),
        axis.text.y = element_text(size=12),
        axis.title.y= element_text(size=14))
ggsave("SPP1hi_Mac_anti_tNS_vs_tS.jpeg",height=4,width=2.5)

ggplot(df,aes(x=TS,y=AP,color=TS))+
  geom_violin(size=0.8)+
  geom_boxplot(width=0.4,outlier.size = 0.4,size=0.8)+
  theme_classic()+
  labs(x="",y="Antigene presentation score")+
  geom_signif(comparisons = list(c("tNS","tS")),step_increase = 0.05,map_signif_level = T,
              test = t.test,size=0.2,textsize = 4,tip_length=0.01,color="grey44")+
  scale_color_manual(values = c("salmon3","salmon4"))+
  theme(legend.position = "none",
        axis.text.x = element_text(size=13),
        axis.text.y = element_text(size=12),
        axis.title.y= element_text(size=14))
ggsave("SPP1hi_Mac_AP_tNS_vs_tS.pdf",height=2.5,width=2.5)
