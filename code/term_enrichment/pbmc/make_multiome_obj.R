library(Signac)
library(Seurat)
library(Matrix)
library(ggplot2)
library(Azimuth)
library(SeuratData)
source("code/utilities.R")

frags <- snakemake@input[['frags']]
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
ref <- readRDS(snakemake@input[['ref']])
ref <- UpdateSeuratObject(ref)

query <- obj

DefaultAssay(query) <- "RNA"
query <- NormalizeData(query)
query <- FindVariableFeatures(query, nfeatures = 2000)
query <- ScaleData(query)
query <- RunPCA(query, verbose=FALSE)

transfer.anchors <- FindTransferAnchors(
  reference = ref,
  query = query,
  reduction = 'pcaproject'
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = ref$celltype,
  weight.reduction = query[['pca']],
  dims = 1:20
)

obj <- AddMetaData(obj, predicted.labels)

saveRDS(obj, snakemake@output[['obj']])


# multiome umap
obj <- NormalizeData(obj, scale.factor = 10000)
obj <- FindVariableFeatures(obj, nfeatures = 2000)
obj <- ScaleData(obj)
obj <- RunPCA(obj, verbose = FALSE)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, verbose = FALSE)

p1 <- DimPlot(obj, group.by = 'predicted.id', reduction = 'umap', label = TRUE, repel = TRUE) + NoLegend() + ggtitle("Label transfer annotations")
ggsave(file = "plots/pbmc_labeltransfer_umap.png", p1, width = 5, height = 5, dpi = 500)

p2 <- ggplot(data.frame(score = obj$prediction.score.max),
             aes(x = score)) +
  geom_histogram(bins = 100) +
  theme_minimal()
ggsave(file = "plots/pbmc_labeltransfer_predscore.png", p2, width = 5, height = 5, dpi = 500)