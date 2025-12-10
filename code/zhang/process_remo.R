library(Signac)
library(Seurat)
source("code/utilities.R")
options(future.globals.maxSize = 100 * 1024 ^ 3)

remo_fetal <- Read10X(snakemake@input[[1]], gene.column = 1)
remo_adult <- Read10X(snakemake@input[[2]], gene.column = 1)

remo <- cbind(remo_adult, remo_fetal)
rm(remo_fetal, remo_adult)
gc()

remo <- NormalizeData(remo)
emb <- RunSVD(remo, pca = TRUE)
saveRDS(emb, snakemake@output[['object']])