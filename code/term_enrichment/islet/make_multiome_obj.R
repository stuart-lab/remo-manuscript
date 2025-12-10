library(Signac)
library(Seurat)
library(Matrix)
library(ggplot2)
library(Azimuth)
library(SeuratData)
source("code/utilities.R")

gex <- snakemake@input[['gex']]
genes <- snakemake@input[['genes']]
rna_cells <- snakemake@input[['allcells']]
qc_cells <- readLines(snakemake@input[['cells']])

gex_mtx <- readMM(gex)
colnames(gex_mtx) <- readLines(rna_cells)
rownames(gex_mtx) <- readLines(genes)
rownames(gex_mtx) <- make.unique(rownames(gex_mtx))

obj <- CreateSeuratObject(counts = gex_mtx[, qc_cells], assay = 'RNA')

# annotate cell types
DefaultAssay(obj) <- 'RNA'
obj <- SCTransform(obj)
obj <- RunPCA(obj)
obj <- RunAzimuth(obj, reference = "pancreasref")

saveRDS(obj, snakemake@output[['obj']])

# multiome umap
obj <- NormalizeData(obj, scale.factor = 10000)
obj <- FindVariableFeatures(obj, nfeatures = 2000)
obj <- ScaleData(obj)
obj <- RunPCA(obj, verbose = FALSE)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, verbose = FALSE)

p1 <- DimPlot(obj, group.by = 'predicted.annotation.l1', reduction = 'umap', label = TRUE, repel = TRUE) + NoLegend() + ggtitle("Label transfer annotations")
ggsave(file = "plots/islet_labeltransfer_umap.png", p1, width = 5, height = 5, dpi = 500)

p2 <- ggplot(data.frame(score = obj$predicted.annotation.l1.score),
             aes(x = score)) +
  geom_histogram(bins = 100) +
  theme_minimal()
ggsave(file = "plots/islet_labeltransfer_predscore.png", p2, width = 5, height = 5, dpi = 500)