# disease metrics
# calculates: silhouette, kNN purity
library(ggplot2)
library(patchwork)
library(Signac)
library(Seurat)
library(cluster)
library(dplyr)
source("code/utilities.R")


### B lymphoma ###
remo <- readRDS(snakemake@input[['B_lymphoma_remo']])
peaks <- readRDS(snakemake@input[['B_lymphoma_peaks']])
loupe_celltypes <- read.table(snakemake@input[['celltype_annot']], sep = ",", header = TRUE)
rownames(loupe_celltypes) <- loupe_celltypes$Barcode

remo <- remo[, colnames(remo) %in% rownames(loupe_celltypes)]
peaks <- peaks[, colnames(peaks) %in% rownames(loupe_celltypes)]

source_celltypes <- loupe_celltypes[colnames(remo), "Cell.Types", drop = FALSE]
remo <- AddMetaData(remo, metadata = source_celltypes)
peaks <- AddMetaData(peaks, metadata = source_celltypes)

# remo
remo_emb <- Embeddings(remo[['pca']])[, 1:30]
idents <- remo$Cell.Types

remo_metric <- list()

# get celltype kNN purity
remo_metric[['knn']] <- data.frame(
        "Score" = knn_purity(remo_emb, idents),
        "Celltype" = idents,
        "Metric" = "kNN purity",
        "Method" = "REMO"
    )
remo_metric[['knn']] <- remo_metric[['knn']][remo_metric[['knn']]$Celltype == "Tumor B",]

# get celltype silhouette score
dist.matrix <- dist(remo_emb)
sil_remo <- cluster::silhouette(x = as.numeric(x = as.factor(x = idents)), dist = dist.matrix)
sil_remo_df <- as.data.frame(sil_remo)

n_cells <- length(which(remo$Cell.Types == 'Tumor B'))
clust_num <- which(summary(sil_remo)$clus.sizes == n_cells)
Tumor_B_sil <- sil_remo_df[sil_remo_df$cluster == clust_num, ]

Tumor_B_sil$Celltype <- "Tumor B"
Tumor_B_sil$Method <- "REMO"
Tumor_B_sil$Metric <- "Silhouette"

remo_metric[['sil']] <- Tumor_B_sil

# get CH index
ch_remo <- fpc::calinhara(remo_emb, as.numeric(as.factor(idents)))

# peaks
peaks_emb <- Embeddings(peaks[['lsi']])[, 2:30]
idents_peaks <- peaks$Cell.Types

peaks_metric <- list()

# get celltype kNN purity
peaks_metric[['knn']] <- data.frame(
        "Score" = knn_purity(peaks_emb, idents_peaks),
        "Celltype" = idents_peaks,
        "Metric" = "kNN purity",
        "Method" = "Peaks"
    )
peaks_metric[['knn']] <- peaks_metric[['knn']][peaks_metric[['knn']]$Celltype == "Tumor B",]

# get celltype silhouette score
dist.matrix <- dist(peaks_emb)
sil_peaks <- cluster::silhouette(x = as.numeric(x = as.factor(x = idents_peaks)), dist = dist.matrix)
sil_peaks_df <- as.data.frame(sil_peaks)

n_cells <- length(which(peaks$Cell.Types == 'Tumor B'))
clust_num <- which(summary(sil_peaks)$clus.sizes == n_cells)
Tumor_B_sil <- sil_peaks_df[sil_peaks_df$cluster == clust_num, ]

Tumor_B_sil$Celltype <- "Tumor B"
Tumor_B_sil$Method <- "Peaks"
Tumor_B_sil$Metric <- "Silhouette"

peaks_metric[['sil']] <- Tumor_B_sil
knn_metrics <- rbind(remo_metric[['knn']], peaks_metric[['knn']])
sil_metrics <- rbind(remo_metric[['sil']], peaks_metric[['sil']])

# get CH
ch_peaks <- fpc::calinhara(peaks_emb, as.numeric(as.factor(idents_peaks)))

ch <- data.frame("Method" = c("REMO", "Peaks"), "Metric" = c("CH", "CH"), "Score" = c(ch_remo, ch_peaks))

B_lymphoma_metrics <- list()
B_lymphoma_metrics[['knn']] <- knn_metrics
B_lymphoma_metrics[['sil']] <- sil_metrics
B_lymphoma_metrics[['ch']] <- ch

saveRDS(B_lymphoma_metrics, snakemake@output[['B_lymphoma_metrics']])

### brain cancer ###
remo <- readRDS(snakemake@input[['brain_cancer_remo']])
peaks <- readRDS(snakemake@input[['brain_cancer_peaks']])

# remo
remo_emb <- Embeddings(remo[['pca']])[, 1:30]
idents <- remo$celltype
idents2 <- remo$orig.ident

sil_list_remo <- list()
knn_list_remo <- list()
batches <- unique(idents2)
for (batch in batches) {
    idx <- which(idents2 == batch)
    batch_emb <- remo_emb[idx, ]
    batch_idents <- idents[idx]
  
    dist.matrix <- dist(batch_emb)
    sil_batch <- cluster::silhouette(x = as.numeric(x = as.factor(x = batch_idents)), dist = dist.matrix)
    sil_batch_df <- as.data.frame(sil_batch)

    n_cells <- sum(batch_idents == 'MLCs')
    clust_num <- which(summary(sil_batch)$clus.sizes == n_cells)
    MLCs_batch_sil <- sil_batch_df[sil_batch_df$cluster == clust_num, ]

    MLCs_batch_sil$Method <- "REMO"
    MLCs_batch_sil$Metric <- "Silhouette"
    MLCs_batch_sil$Celltype <- "MLCs"
    sil_list_remo[[batch]] <- MLCs_batch_sil

    knn_batch <- data.frame(
        "Score" = knn_purity(batch_emb, batch_idents),
        "Celltype" = batch_idents,
        "Metric" = "kNN purity",
        "Method" = "REMO"
    )
    knn_list_remo[[batch]] <- knn_batch
}

remo_sil_I2 <- sil_list_remo[["I2"]]
remo_knn_I2 <- knn_list_remo[["I2"]][knn_list_remo[["I2"]]$Celltype == "MLCs", ]
remo_sil_M7 <- sil_list_remo[["M7"]]
remo_knn_M7 <- knn_list_remo[["M7"]][knn_list_remo[["M7"]]$Celltype == "MLCs", ]


remo$cluster <- paste0(remo$celltype, remo$orig.ident)
ch_remo <- fpc::calinhara(remo_emb, as.numeric(as.factor(remo$cluster)))

# peaks
peaks_emb <- Embeddings(peaks[['lsi']])[, 2:30]
idents_peaks <- peaks$celltype
idents2_peaks <- peaks$orig.ident
peaks$cluster <- paste0(peaks$celltype, peaks$orig.ident)

sil_list_peaks <- list()
knn_list_peaks <- list()
batches <- unique(idents2_peaks)
for (batch in batches) {
    idx <- which(idents2_peaks == batch)
    batch_emb <- peaks_emb[idx, ]
    batch_idents <- idents[idx]
  
    dist.matrix <- dist(batch_emb)
    sil_batch <- cluster::silhouette(x = as.numeric(x = as.factor(x = batch_idents)), dist = dist.matrix)
    sil_batch_df <- as.data.frame(sil_batch)

    n_cells <- sum(batch_idents == 'MLCs')
    clust_num <- which(summary(sil_batch)$clus.sizes == n_cells)
    MLCs_batch_sil <- sil_batch_df[sil_batch_df$cluster == clust_num, ]

    MLCs_batch_sil$Method <- "Peaks"
    MLCs_batch_sil$Metric <- "Silhouette"
    MLCs_batch_sil$Celltype <- "MLCs"
    sil_list_peaks[[batch]] <- MLCs_batch_sil

    knn_batch <- data.frame(
        "Score" = knn_purity(batch_emb, batch_idents),
        "Celltype" = batch_idents,
        "Metric" = "kNN purity",
        "Method" = "Peaks"
    )
    knn_list_peaks[[batch]] <- knn_batch
}

peaks_sil_I2 <- sil_list_peaks[["I2"]]
peaks_knn_I2 <- knn_list_peaks[["I2"]][knn_list_peaks[["I2"]]$Celltype == "MLCs", ]
peaks_sil_M7 <- sil_list_peaks[["M7"]]
peaks_knn_M7 <- knn_list_peaks[["M7"]][knn_list_peaks[["M7"]]$Celltype == "MLCs", ]

ch_peaks <- fpc::calinhara(peaks_emb, as.numeric(as.factor(peaks$cluster)))

ch <- data.frame("Method" = c("REMO", "Peaks"), "Metric" = c("CH", "CH"), "Score" = c(ch_remo, ch_peaks))

brain_cancer_metrics <- list()
brain_cancer_metrics[['knn']] <- rbind(remo_knn_I2, remo_knn_M7, peaks_knn_I2, peaks_knn_M7)
brain_cancer_metrics[['sil']] <- rbind(remo_sil_I2, remo_sil_M7, peaks_sil_I2, peaks_sil_M7)
brain_cancer_metrics[['ch']] <- ch

saveRDS(brain_cancer_metrics, snakemake@output[['brain_cancer_metrics']])


### plot metrics ###
all_metrics_knn <- rbind(B_lymphoma_metrics[['knn']], brain_cancer_metrics[['knn']])
all_metrics_sil <- rbind(B_lymphoma_metrics[['sil']], brain_cancer_metrics[['sil']])

knn_metrics <- all_metrics_knn %>%
  group_by(Celltype, Method, Metric)

knn_metrics$Method <- factor(knn_metrics$Method, levels = c("REMO", "Peaks"))

knn_plot <- ggplot(knn_metrics, aes(Celltype, Score, fill = Method)) +
  geom_boxplot() +
  theme_minimal() +
  ylab("kNN purity score") +
  xlab(" ") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none') +
  scale_fill_manual(values = c(REMO_COL, PEAKS_COL))

sil_metrics <- all_metrics_sil %>%
  group_by(Celltype, Method, Metric)

sil_metrics$Method <- factor(sil_metrics$Method, levels = c("REMO", "Peaks"))

sil_plot <- ggplot(sil_metrics, aes(Celltype, sil_width, fill = Method)) +
  geom_boxplot(linewidth=0.5) +
  theme_minimal() +
  ylab("Silhouette score") +
  xlab(" ") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'right') +
  scale_fill_manual(values = c(REMO_COL, PEAKS_COL))

ch_b <- B_lymphoma_metrics[['ch']]
ch_b$Dataset <- "B-cell lymphoma"
ch_brain <- brain_cancer_metrics[['ch']]
ch_brain$Dataset <- "Brain cancer"

ch <- rbind(ch_b, ch_brain)
ch$Method <- factor(ch$Method, levels = c("REMO", "Peaks"))
ch$Dataset <- factor(ch$Dataset, levels = c("Brain cancer", "B-cell lymphoma"))

ch_plot <- ggplot(ch, aes(Dataset, Score, fill = Method)) +
  geom_bar(stat = "identity", position = 'dodge') +
  xlab("") + 
  ylab("Calinski-Harabasz index") +
  theme_minimal() +
  scale_fill_manual(values = c(REMO_COL, PEAKS_COL)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none')

disease_metrics_plot <- ch_plot | knn_plot | sil_plot

ggsave(filename = snakemake@output[['disease_metrics_plot']], plot = disease_metrics_plot, height = 4, width = 8, dpi = 500)
