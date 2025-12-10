# QC islet
library(Signac)
library(Seurat)
library(Matrix)
library(ggplot2)
library(scDblFinder)
source("code/utilities.R")

frags <- snakemake@input[['frags']]
gex <- snakemake@input[['gex']]
genes <- snakemake@input[['genes']]
rna_cells <- snakemake@input[['cells']]
annotations <- readRDS(snakemake@input[['annotations']])

gex_mtx <- readMM(gex)
rownames(gex_mtx) <- readLines(genes)
colnames(gex_mtx) <- readLines(rna_cells)
rownames(gex_mtx) <- make.unique(rownames(gex_mtx))

md <- ATACqc(object = frags, annotations = annotations) 
obj <- CreateSeuratObject(counts = gex_mtx, assay = 'RNA', meta.data = md)

DefaultAssay(obj) <- 'RNA'
obj <- NormalizeData(obj)
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

obj_filt <- subset(obj, subset = total_fragments > 5000 &
                       total_fragments < 100000 &
                       TSS_enrichment > 9 &
                       nCount_RNA > 1000 &
                       percent.mt < 30)

# doublet score
DefaultAssay(obj_filt) <- "RNA"
scDbl <- as.SingleCellExperiment(obj_filt)
scDbl <- scDblFinder(scDbl)

obj_filt$dbl.class <- scDbl$scDblFinder.class
obj_filt <- obj_filt[, obj_filt$dbl.class == "singlet"]

writeLines(con = snakemake@output[['cells']], text = colnames(obj_filt))
