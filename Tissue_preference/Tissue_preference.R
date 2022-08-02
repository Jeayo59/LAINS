library(colorspace)
## tissue preference
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
compare <- function(data,cellType,h){

  df <- data[data$CellType==cellType,]
  df
   
  p.nNSvstNS <- wilcox.test(df[df$TS=="nNS",]$value,df[df$TS=="tNS",]$value,paired = T)$p.value %>% judgeDegree()
  p.nSvstS <- wilcox.test(df[df$TS=="nS",]$value,df[df$TS=="tS",]$value)$p.value %>% judgeDegree()
  p.TvsN <- wilcox.test(df[df$Type=="Tumor",]$value,df[df$Type=="Normal",]$value,paired = T)$p.value %>% judgeDegree()
  
  if(p.TvsN != "NS."){
    p =  wilcox.test(df[df$Type=="Tumor",]$value,df[df$Type=="Normal",]$value,alternative="greater")$p.value
    if(p<0.1){height=h/10}
    else{height=-1*h/10}
  }
  else{height=0}
  
  annotation_table <- anno(h,height)[[1]]
  annotation_table1 <- anno(h,height)[[2]]
  
  loc <- data.frame(x=c(3,4,2.5),
                    y=c(h-h/25+height,h-h/25+height,h+(h/10+abs(height))),
                    value=c(p.nNSvstNS,p.nSvstS,p.TvsN),
                    type=c(2,2,4))
  df
  loc
  ggplot(df)+
    geom_boxplot(aes(x=TS,y=value,color=TS),outlier.size = 0.5)+
    geom_point(aes(x=TS,y=value,color=Patient),size=1)+
    geom_segment(inherit.aes = F,data= annotation_table,aes(x=xsl,xend=xel,y=ysl,yend=yel),color="grey40")+
    geom_segment(inherit.aes = F,data= annotation_table,aes(x=xsr,xend=xer,y=ysr,yend=yer),color="grey40")+
    geom_segment(inherit.aes = F,data= annotation_table,aes(x=xsh,xend=xeh,y=ysh,yend=yeh),color="grey40")+
    geom_segment(inherit.aes = F,data= annotation_table,aes(x=xsu,xend=xeu,y=ysu,yend=yeu),color="grey40")+
    geom_segment(inherit.aes = F,data= annotation_table1,aes(x=xsuh,xend=xeuh,y=ysuh,yend=yeuh),color="grey40")+
    geom_text(data=loc,aes(x=x,y=y,label=value,color=as.factor(type)))+
    scale_color_manual(values = c("green","red",
                                  "darkseagreen3","darkseagreen4",
                                  Patient_color[unique(df$Patient),]$color,
                                  "salmon3","salmon4"))+
    labs(x="",y="o/e",title=cellType)+
    theme_classic()+
    theme(plot.title = element_text(size=14,hjust=0.5,vjust=0.5),
          legend.position = "none",
          axis.text.x = element_text(size=13),
          axis.text.y = element_text(size=12),
          axis.title.y= element_text(size=14))
}


setwd("I:/交接/LUAD_merge/09.tissue_preference/")
df_main <- read.table("Main_type_statistics.txt",header=FALSE,check.names = FALSE)
colnames(df_main) <- c("CellType","SmokingType","Patient","tumorNumber","normalNumber")
library(dplyr)
total <- ddply(df_main,"Patient",summarise,tumorTotal=sum(tumorNumber),normalTotal=sum(normalNumber))
df_main <- merge(df_main,total,by="Patient",all.x=TRUE)
df_main$tumorPer <- df_main$tumorNumber/df_main$tumorTotal
df_main$normalPer <- df_main$normalNumber/df_main$normalTotal
## observe/expect
tumorE <- c()
normalE <- c()
for(i in 1:length(df_main[,1])){
  tumorE <- c(tumorE,chisq.test(df_main[i,c(8,9)])$expect[1])
  normalE <- c(normalE,chisq.test(df_main[i,c(8,9)])$expect[2])
}
df_main$tumorE <- tumorE
df_main$normalE <- normalE
head(df_main)
df_main$tumor_oe <- df_main$tumorPer/df_main$tumorE
df_main$normal_oe <- df_main$normalPer/df_main$normalE
head(df_main)
df_plot <- df_main[,c("Patient","CellType","SmokingType","tumor_oe","normal_oe")]
head(df_plot)
df_plot.melt <- reshape2::melt(df_plot,id.vars=c("Patient","CellType","SmokingType"))
head(df_plot.melt)
df_plot.melt <- tidyr::unite(df_plot.melt,"TS",variable,SmokingType)
df_plot.melt$TS <- sub("tumor_oe_NS","tNS",df_plot.melt$TS)
df_plot.melt$TS <- sub("tumor_oe_S","tS",df_plot.melt$TS)
df_plot.melt$TS <- sub("normal_oe_NS","nNS",df_plot.melt$TS)
df_plot.melt$TS <- sub("normal_oe_S","nS",df_plot.melt$TS)
df_plot.melt$Type[df_plot.melt$TS=="tNS" |
                    df_plot.melt$TS=="tS"] <- "Tumor"

df_plot.melt$Type[df_plot.melt$TS=="nNS" |
                    df_plot.melt$TS=="nS"] <- "Normal"

#Patient_color <- data.frame(Patient=unique(df_plot.melt$Patient),
#                            color=diverging_hcl(17, "Tropic"))
Patient_color <- data.frame(Patient=unique(df_plot.melt$Patient),
                            color=rep("black", 17))
rownames(Patient_color) <- Patient_color$Patient

head(df_plot.melt)

setwd("./Main/")
for(cell in unique(df_plot.melt$CellType)){
  compare(df_plot.melt,cell,2.5)
  ggsave(paste0(cell,"_tissue_preference.pdf"),height=4,width=3)
}

library(reshape2)
setwd("../")
df_sub <- read.table("Sub_type_statistics.txt",header=FALSE,check.names = FALSE)
head(df_sub)
colnames(df_sub) <- c("CellType","Patient","SmokingType","tumorNumber","normalNumber")
total <- ddply(df_sub,"Patient",summarise,tumorTotal=sum(tumorNumber),normalTotal=sum(normalNumber))
df_sub <- merge(df_sub,total,by="Patient",all.x=TRUE)

df_sub$tumorPer <- df_sub$tumorNumber/df_sub$tumorTotal
df_sub$normalPer <- df_sub$normalNumber/df_sub$normalTotal

## observe/expect
df_sub <- df_sub[df_sub$tumorNumber!=0 & df_sub$normalNumber!=0,]

tumorE <- c()
normalE <- c()
for(i in 1:length(df_sub[,1])){
  tumorE <- c(tumorE,chisq.test(df_sub[i,c(8,9)])$expect[1])
  normalE <- c(normalE,chisq.test(df_sub[i,c(8,9)])$expect[2])
}
df_sub$tumorE <- tumorE
df_sub$normalE <- normalE
head(df_sub)
df_sub$tumor_oe <- df_sub$tumorPer/df_sub$tumorE
df_sub$normal_oe <- df_sub$normalPer/df_sub$normalE
head(df_sub)
df_plot <- df_sub[,c("Patient","CellType","SmokingType","tumor_oe","normal_oe")]
head(df_plot)
df_plot.melt <- reshape2::melt(df_plot,id.vars=c("Patient","CellType","SmokingType"))
head(df_plot.melt)
df_plot.melt <- tidyr::unite(df_plot.melt,"TS",variable,SmokingType)
df_plot.melt$TS <- sub("tumor_oe_NS","tNS",df_plot.melt$TS)
df_plot.melt$TS <- sub("tumor_oe_S","tS",df_plot.melt$TS)
df_plot.melt$TS <- sub("normal_oe_NS","nNS",df_plot.melt$TS)
df_plot.melt$TS <- sub("normal_oe_S","nS",df_plot.melt$TS)
df_plot.melt$Type[df_plot.melt$TS=="tNS" |
                    df_plot.melt$TS=="tS"] <- "Tumor"

df_plot.melt$Type[df_plot.melt$TS=="nNS" |
                    df_plot.melt$TS=="nS"] <- "Normal"

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
unique(df_plot.melt$Patient)
Patient_color <- data.frame(Patient=unique(df_plot.melt$Patient),
                            color=diverging_hcl(17, "Tropic"))
rownames(Patient_color) <- Patient_color$Patient



setwd("I:/交接/LUAD_merge/09.tissue_preference/Sub/")
for(cell in unique(df_plot.melt$CellType)[-c(36,40)]){
  compare(df_plot.melt,cell,2.5)
  ggsave(paste0(cell,"_tissue_preference.pdf"),height=4,width=3)
}




