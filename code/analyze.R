library(Signac)
library(Seurat)
library(Matrix)
library(patchwork)
library(ggplot2)
source("code/utilities.R")

gex <- snakemake@input[['gex']]
genes <- snakemake@input[['genes']]
rna_cells <- snakemake@input[['rna_cells']]
qc_cells <- snakemake@input[['qc_cells']]
peak_counts_dir <- paste0(snakemake@input[['peak_counts']], "/")
remo_counts_dir <- paste0(snakemake@input[['remo_counts']], "/")
annotations <- snakemake@input[['annotations']]
remo <- snakemake@input[['remo']]

dims <- as.numeric(snakemake@params[["dims"]])
resolution <- 0.5

grouped <- Read10X(remo_counts_dir, gene.column=1)
peak_counts <- Read10X(peak_counts_dir, gene.column=1)
gex_counts <- as(readMM(gex), "CsparseMatrix")
rownames(gex_counts) <- readLines(genes)
colnames(gex_counts) <- readLines(rna_cells)

cells.use <- readLines(qc_cells)

grouped <- grouped[, cells.use]
peak_counts <- peak_counts[, cells.use]
gex_counts <- gex_counts[, cells.use]

# make gex rownames unique
rownames(gex_counts) <- make.unique(rownames(gex_counts))

remo_obj <- process_remo_obj(grouped, dims = 1:dims, resolution = resolution)
peak_obj <- process_atac_obj(peak_counts, dims = 2:dims, resolution = resolution)
gex_obj <- process_rna_obj(gex_counts, dims = 1:dims, resolution = resolution)

remo_obj <- AddMetaData(remo_obj, gex_obj$seurat_clusters, col.name = "cluster")
peak_obj <- AddMetaData(peak_obj, gex_obj$seurat_clusters, col.name = "cluster")

saveRDS(remo_obj, snakemake@output[['remo']])
saveRDS(peak_obj, snakemake@output[['peaks']])
saveRDS(gex_obj, snakemake@output[['gex']])