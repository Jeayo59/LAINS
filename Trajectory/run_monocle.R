Args <- commandArgs(T)
library(Seurat)
library(monocle)
LUAD <- readRDS(Args[1])

data <- as(as.matrix(LUAD@assays$RNA@counts), 'sparseMatrix')

pd <- LUAD@meta.data

pd <- new('AnnotatedDataFrame', data = pd)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

monocle_cds <- newCellDataSet(as.matrix(data),
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())

dat <- monocle_cds
dat <- estimateSizeFactors(dat)
dat <- estimateDispersions(dat)

disp_table <- dispersionTable(dat)

unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
unsup_clustering_genes

dat <- setOrderingFilter(dat, unsup_clustering_genes$gene_id)

diff_genes <- differentialGeneTest(dat, fullModelFormulaStr = "~seurat_clusters",cores = 10)
head(diff_genes)
ordering_genes <- row.names(subset(diff_genes, qval < 0.01))
dat <- setOrderingFilter(dat, ordering_genes)
dat <- reduceDimension(dat, max_components = 2, method = 'DDRTree')
dat <- orderCells(dat)
pdf("trajectory_plot_State.pdf")
plot_cell_trajectory(dat,color_by = "State")
dev.off()

pdf("trajectory_plot_Pseudotime.pdf")
plot_cell_trajectory(dat,color_by = "Pseudotime")
dev.off()

pdf("trajectory_plot_seurat_cluster.pdf")
plot_cell_trajectory(dat, color_by = "seurat_clusters", cell_size = 1) +
  scale_color_brewer(name = "sub type", type = "qual", palette = "Set2")
dev.off()

saveRDS(dat,"monocle.rds")

