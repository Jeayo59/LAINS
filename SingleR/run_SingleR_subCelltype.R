Args <- commandArgs()
library(SingleR)
library(scRNAseq)
LUAD <- readRDS(Args[6])
Cell_type <- Args[7]

Count <- LUAD@assays$RNA@counts
hpca.se <- celldex::HumanPrimaryCellAtlasData()
bp.se <- BlueprintEncodeData()
## run SingleR
pred.hesc <- SingleR(test=Count, ref=list(BP=bp.se,HPCA=hpca.se),
                     labels=list(bp.se$label.fine,hpca.se$label.fine),
                     clusters = LUAD@meta.data$seurat_clusters)

## heatmap
pdf.file <- paste0(Cell_type,"_SingleR.pdf")
pdf(pdf.file,height=10,width=25)
plotScoreHeatmap(pred.hesc)
dev.off()

dat <- as.data.frame(pred.hesc[,c(1:5)])
dat <- dat[,(ncol(dat)-3):ncol(dat)]
dat$cluster <- rownames(dat)
dat <- dat[,c(5,1,2,3,4)]
file.name <- paste0(Cell_type,"_SingleR_result.txt")
write.table(dat,file.name,sep="\t",quote = FALSE,row.names = FALSE)

