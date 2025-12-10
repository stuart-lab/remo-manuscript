library(Signac)
library(Seurat)
library(Azimuth)
library(SeuratData)
source("code/utilities.R")


gex_obj <- readRDS(snakemake@input[['obj']])
DefaultAssay(gex_obj) <- 'RNA'
gex_obj <- RunAzimuth(gex_obj, reference = "humancortexref")

saveRDS(gex_obj, snakemake@output[['obj']])

# multiome umap
library(ggplot2)
obj <- gex_obj
obj <- NormalizeData(obj, scale.factor = 10000)
obj <- FindVariableFeatures(obj, nfeatures = 2000)
obj <- ScaleData(obj)
obj <- RunPCA(obj, verbose = FALSE)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, verbose = FALSE)

p1 <- DimPlot(obj, group.by = 'predicted.subclass', reduction = 'umap', label = TRUE, repel = TRUE) + NoLegend() + ggtitle("Label transfer annotations")
ggsave(file = "plots/brain_labeltransfer_umap.png", p1, width = 5, height = 5, dpi = 500)

p2 <- ggplot(data.frame(score = obj$predicted.subclass.score),
             aes(x = score)) +
  geom_histogram(bins = 100) +
  theme_minimal()
ggsave(file = "plots/brain_labeltransfer_predscore.png", p2, width = 5, height = 5, dpi = 500)