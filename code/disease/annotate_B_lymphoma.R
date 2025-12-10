# annotate celltypes
library(ggplot2)
library(patchwork)
library(Signac)
library(Seurat)
source("code/utilities.R")

### 10x B lymphoma ###
# load data
gex <- readRDS(snakemake@input[['gex']])
remo <- readRDS(snakemake@input[['remo']])
peaks <- readRDS(snakemake@input[['peaks']])
loupe_celltypes <- read.table(snakemake@input[['celltype_annot']], sep = ",", header = TRUE)

# annotate cell types
rownames(loupe_celltypes) <- loupe_celltypes$Barcode

source_celltypes <- loupe_celltypes[colnames(remo), "Cell.Types", drop = FALSE]
remo <- AddMetaData(remo, metadata = source_celltypes)

source_celltypes <- loupe_celltypes[colnames(gex), "Cell.Types", drop = FALSE]
gex <- AddMetaData(gex, metadata = source_celltypes)

source_celltypes <- loupe_celltypes[colnames(peaks), "Cell.Types", drop = FALSE]
peaks <- AddMetaData(peaks, metadata = source_celltypes)

saveRDS(gex, "objects/10x_B_lymphoma_multiome_gex.rds")
saveRDS(remo, "objects/10x_B_lymphoma_multiome_remo.rds")
saveRDS(peaks, "objects/10x_B_lymphoma_multiome_peaks.rds")

# dimplots
umap_gex <- DimPlot(gex, group.by = 'Cell.Types') +
  ggtitle("Gene expression") + theme_dimplot(legend = TRUE)

umap_peaks <- DimPlot(peaks, group.by = 'Cell.Types') +
  ggtitle("Peaks") + theme_dimplot(legend = TRUE)

umap_remo <- DimPlot(remo, group.by = 'Cell.Types') +
  ggtitle("") + theme_dimplot(legend = FALSE)

ggsave(file = snakemake@output[['umap_remo']], umap_remo, width = 4, height = 4, dpi = 500)
ggsave(file = snakemake@output[['umap_peaks']], umap_peaks, width = 5, height = 5, dpi = 500)
ggsave(file = snakemake@output[['umap_gex']], umap_gex, width = 5, height = 5, dpi = 500)