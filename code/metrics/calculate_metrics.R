library(Signac)
library(Seurat)
library(tidyr)
library(dplyr)
library(ggplot2)
source("code/utilities.R")

remo_dir <- snakemake@input[['remo']]
peaks_dir <- snakemake@input[['peaks']]

dims <- as.numeric(snakemake@params[["dims"]])
sample_name <- as.character(snakemake@params[['tissue_name']])

remo <- readRDS(remo_dir)
peaks <- readRDS(peaks_dir)

metrics <- get_metrics(
    remo_obj = remo,
    peaks_obj = peaks,
    remo_emb = Embeddings(remo[['pca']])[, 1:dims],
    peaks_emb = Embeddings(peaks[['lsi']])[, 2:dims],
    tissue_name = sample_name
)

saveRDS(metrics, snakemake@output[['metrics']])