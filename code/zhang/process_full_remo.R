library(Seurat)
library(ggplot2)
library(RcppHNSW)
library(uwot)
source("code/utilities.R")

# PCA embeddings
adult <- Read10X(snakemake@input[['adult']], gene.column = 1)
fetal <- Read10X(snakemake@input[['fetal']], gene.column = 1)

counts <- cbind(adult, fetal)
rm(adult, fetal)
gc()

counts <- NormalizeData(counts)
emb <- RunSVD(counts, pca = TRUE)
saveRDS(emb, snakemake@output[['object']])