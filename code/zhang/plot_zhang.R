library(Seurat)
library(ggplot2)
library(RcppHNSW)
library(uwot)
source("code/utilities.R")

# PCA embeddings
atlas <- readRDS(snakemake@input[['emb']])

# run UMAP
umap.emb <- umap2(Embeddings(atlas), nn_method = "hnsw")
umap.emb <- as.data.frame(umap.emb)

umap.emb$Dataset <- c(rep("Adult", 600000), rep("Fetal", 600000))

p <- ggplot(umap.emb, aes(V1, V2, color = Dataset)) + geom_point(size=0.1) + theme_dimplot()
ggsave(snakemake@output[['plot']], p, height = 25, width = 25)