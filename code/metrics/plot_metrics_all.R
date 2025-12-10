library(Signac)
library(Seurat)
library(tidyr)
library(dplyr)
library(ggplot2)
source("code/utilities.R")

liver_metrics <- readRDS(snakemake@input[['liver']])
brain_metrics <- readRDS(snakemake@input[['brain']])
pbmc_metrics <- readRDS(snakemake@input[['pbmc']])
jejunum_metrics <- readRDS(snakemake@input[['jejunum']])
muscle_metrics <- readRDS(snakemake@input[['muscle']])
pancreas_metrics <- readRDS(snakemake@input[['pancreas']])
lung_metrics <- readRDS(snakemake@input[['lung']])
heart_metrics <- readRDS(snakemake@input[['heart']])
left_colon_metrics <- readRDS(snakemake@input[['left_colon']])
bile_duct_metrics <- readRDS(snakemake@input[['bile_duct']])
fallopian_tube_metrics <- readRDS(snakemake@input[['fallopian_tube']])
grey_matter_metrics <- readRDS(snakemake@input[['grey_matter']])
kidney_cortex_metrics <- readRDS(snakemake@input[['kidney_cortex']])
pituitary_metrics <- readRDS(snakemake@input[['pituitary']])
placenta_metrics <- readRDS(snakemake@input[['placenta']])
retina_metrics <- readRDS(snakemake@input[['retina']])

datasets <- c("Brain", "PBMC", "Jejunum", "Heart", "Left colon", "Liver",
                "Lung", "Muscle", "Pancreas", "Bile duct", "Fallopian tube", "Grey matter",
                "Kidney", "Pituitary gland", "Placenta", "Retina", "B lymphoma")

names(datasets) <- c("brain", "pbmc", "jejunum", "heart", "left_colon", "liver",
                "lung", "muscle", "pancreas", "bile_duct", "fallopian_tube", "grey_matter",
                "kidney_cortex", "pituitary", "placenta", "retina", "10x_B_lymphoma")


# combine
all_metrics <- rbind(
    liver_metrics$AN,
    brain_metrics$AN,
    pbmc_metrics$AN,
    jejunum_metrics$AN, 
    muscle_metrics$AN,
    pancreas_metrics$AN,
    lung_metrics$AN,
    heart_metrics$AN,
    left_colon_metrics$AN,
    bile_duct_metrics$AN,
    fallopian_tube_metrics$AN,
    grey_matter_metrics$AN,
    kidney_cortex_metrics$AN,
    pituitary_metrics$AN,
    placenta_metrics$AN,
    retina_metrics$AN
)

all_sil <- rbind(
    liver_metrics$SIL,
    brain_metrics$SIL,
    pbmc_metrics$SIL,
    jejunum_metrics$SIL, 
    muscle_metrics$SIL,
    pancreas_metrics$SIL,
    lung_metrics$SIL,
    heart_metrics$SIL,
    left_colon_metrics$SIL,
    bile_duct_metrics$SIL,
    fallopian_tube_metrics$SIL,
    grey_matter_metrics$SIL,
    kidney_cortex_metrics$SIL,
    pituitary_metrics$SIL,
    placenta_metrics$SIL,
    retina_metrics$SIL
)

all_knn <- rbind(
    liver_metrics$KNN,
    brain_metrics$KNN,
    pbmc_metrics$KNN,
    jejunum_metrics$KNN, 
    muscle_metrics$KNN,
    pancreas_metrics$KNN,
    lung_metrics$KNN,
    heart_metrics$KNN,
    left_colon_metrics$KNN,
    bile_duct_metrics$KNN,
    fallopian_tube_metrics$KNN,
    grey_matter_metrics$KNN,
    kidney_cortex_metrics$KNN,
    pituitary_metrics$KNN,
    placenta_metrics$KNN,
    retina_metrics$KNN
)

all_metrics$Dataset <- datasets[all_metrics$Dataset]
all_sil$Dataset <- datasets[all_sil$Dataset]
all_knn$Dataset <- datasets[all_knn$Dataset]

saveRDS(all_metrics, snakemake@output[['all_metrics']])
saveRDS(all_sil, snakemake@output[['all_sil']])
saveRDS(all_knn, snakemake@output[['all_knn']])

metrics <- all_metrics
metrics_sil <- all_sil
metrics_knn <- all_knn

metrics$Method <- factor(metrics$Method, levels = c("REMO", "Peaks"))
metrics$Dataset <- factor(metrics$Dataset, levels = datasets)

metrics_sil$Method <- factor(metrics_sil$Method, levels = c("REMO", "Peaks"))
metrics_knn$Method <- factor(metrics_knn$Method, levels = c("REMO", "Peaks"))
metrics_sil$Dataset <- factor(metrics_sil$Dataset, levels = datasets)
metrics_knn$Dataset <- factor(metrics_knn$Dataset, levels = datasets)

sil_plot <- ggplot(metrics_sil, aes(Dataset, mn, fill = Method)) +
  geom_boxplot() +
  theme_minimal() +
  ylab("Mean Silhouette") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none') +
  scale_fill_manual(values = c(REMO_COL, PEAKS_COL))

knn_plot <- ggplot(metrics_knn, aes(Dataset, mn, fill = Method)) +
  geom_boxplot() +
  theme_minimal() +
  ylab("Mean kNN purity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none') +
  scale_fill_manual(values = c(REMO_COL, PEAKS_COL))

ch_plot <- ggplot(metrics[metrics$Metric == 'CH', ], aes(x = Dataset, y = Score, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  scale_y_continuous(name = NULL, sec.axis = dup_axis(name = "Calinski–Harabas\nIndex")) +
  theme(
    axis.title.y.right = element_text(angle = -90, vjust = 0.5),
    axis.text.y.right  = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "none"
  ) +
  scale_fill_manual(values = c(REMO_COL, PEAKS_COL))


metrics_knn$Metric <- "Mean\nkNN purity"
metrics_sil$Metric <- "Mean\nSilhouette"

combined_metrics <- rbind(metrics_knn, metrics_sil)

pp <- ggplot(combined_metrics, aes(Dataset, mn, fill = Method)) +
  geom_boxplot(linewidth=0.2, outlier.size=0.5) +
  theme_minimal() +
  ylab("") + xlab("") +
  facet_wrap(~Metric, ncol=1, strip.position = 'right', scales = "free_y") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 12),
        strip.text = element_text(size = 12),
        legend.position = 'none') +
  scale_fill_manual(values = c(REMO_COL, PEAKS_COL))

ggsave(filename = snakemake@output[['plots_combined']], plot = pp, height = 3, width = 14)
ggsave(filename = snakemake@output[['plots_ch']], plot = ch_plot, height = 2, width = 14)
ggsave(filename = snakemake@output[['plots_sil']], plot = sil_plot, height = 3, width = 18, dpi = 500)
ggsave(filename = snakemake@output[['plots_knn']], plot = knn_plot, height = 3, width = 18, dpi = 500)

# PBMC only
pbmc_metrics <- metrics[metrics$Dataset == "PBMC", ]
combined_metrics <- combined_metrics[combined_metrics$Dataset == "PBMC", ]

knn <- ggplot(combined_metrics[combined_metrics$Metric == "Mean kNN purity", ], aes(Metric, mn, fill = Method)) +
  geom_boxplot(linewidth=0.2, outlier.size=0.5) +
  theme_minimal() +
  ylab("Score") + xlab("") +
  scale_fill_manual(values = c(REMO_COL, PEAKS_COL))

sil <- ggplot(combined_metrics[combined_metrics$Metric == "Mean Silhouette", ], aes(Metric, mn, fill = Method)) +
  geom_boxplot(linewidth=0.2, outlier.size=0.5) +
  theme_minimal() +
  ylab("Score") + xlab("") +
  scale_fill_manual(values = c(REMO_COL, PEAKS_COL))

ggsave(filename = "plots/sil_pbmc.pdf", plot = sil, height = 1.8, width = 3, dpi = 1000)
ggsave(filename = "plots/knn_pbmc.pdf", plot = knn, height = 1.8, width = 3, dpi = 1000)