library(Signac)
library(Seurat)
library(Matrix)
library(patchwork)
library(ggplot2)
library(scDblFinder)
source("code/utilities.R")

# inputs
frags <- snakemake@input[['frags']]
gex <- snakemake@input[['gex']]
genes <- snakemake@input[['genes']]
rna_cells <- snakemake@input[['rna_cells']]
peak_counts_dir <- snakemake@input[['peak_counts']]
annotations <- readRDS(snakemake@input[['annotations']])

# QC params
nCount_ATAC_above <- as.numeric(snakemake@params[["nCount_ATAC_above"]])
nCount_ATAC_below <- as.numeric(snakemake@params[["nCount_ATAC_below"]])
TSS_above <- as.numeric(snakemake@params[["TSS_above"]])
nCount_RNA_above <- as.numeric(snakemake@params[["nCount_RNA_above"]])
nCount_RNA_below <- as.numeric(snakemake@params[["nCount_RNA_below"]])
percent_mt_below <- as.numeric(snakemake@params[["percent_mt_below"]])

gex_mtx <- readMM(gex)
rownames(gex_mtx) <- readLines(genes)
colnames(gex_mtx) <- readLines(rna_cells)
rownames(gex_mtx) <- make.unique(rownames(gex_mtx))

peakcounts <- Read10X(peak_counts_dir, gene.column=1)

common_cells <- intersect(colnames(peakcounts), colnames(gex_mtx))

obj <- CreateSeuratObject(counts = gex_mtx[, common_cells], assay = 'RNA')

obj[['ATAC']] <- CreateChromatinAssay(
    counts = peakcounts[, common_cells],
    fragments = frags,
    annotation = annotations
)

DefaultAssay(obj) <- 'ATAC'

obj <- ATACqc(obj)

DefaultAssay(obj) <- 'RNA'

obj <- NormalizeData(obj)
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

obj_filt <- subset(obj, subset = nCount_ATAC > nCount_ATAC_above &
                       nCount_ATAC < nCount_ATAC_below &
                       TSS_enrichment > TSS_above &
                       nCount_RNA > nCount_RNA_above &
                       nCount_RNA < nCount_RNA_below &
                       percent.mt < percent_mt_below)

# doublet score
DefaultAssay(obj_filt) <- "RNA"
scDbl <- as.SingleCellExperiment(obj_filt)
scDbl <- scDblFinder(scDbl)

obj_filt$dbl.class <- scDbl$scDblFinder.class
obj_filt <- obj_filt[, obj_filt$dbl.class == "singlet"]

writeLines(con = snakemake@output[['qc_cells']], text = colnames(obj_filt))