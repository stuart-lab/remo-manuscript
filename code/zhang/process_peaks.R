library(Signac)
library(Seurat)
source("code/utilities.R")
options(future.globals.maxSize = 100 * 1024 ^ 3)

peaks_fetal <- Read10X(snakemake@input[[1]], gene.column = 1)
peaks_adult <- Read10X(snakemake@input[[2]], gene.column = 1)

peaks <- cbind(peaks_adult, peaks_fetal)
rm(peaks_fetal, peaks_adult)
gc()

peaks <- RunTFIDF(peaks)
emb <- RunSVD(peaks)
saveRDS(emb, snakemake@output[['object']])