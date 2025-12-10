library(Signac)
library(Seurat)
library(ggplot2)
library(Matrix)
library(dplyr)
library(sparseMatrixStats)
library(REMO.v1.GRCh38)
source('code/utilities.R')
options(scipen = 999)
set.seed(1234)


remo <- Read10X('data/pbmc/remo_multiome/', gene.column=1)
peaks <- Read10X('data/pbmc/peaks_multiome/', gene.column=1)
gex <- as(readMM('data/pbmc/gex.mtx'), "CsparseMatrix")
rownames(gex) <- readLines('data/pbmc/genes.tsv')
colnames(gex) <- readLines('data/pbmc/rna_cells.txt')
qc_cells <- readLines("data/pbmc/cells_pbmc_multiome.txt")

remo <- remo[, qc_cells]
peaks <- peaks[, qc_cells]
gex <- gex[, qc_cells]

# number of CREs per module
# distance between elements within module
remo_regions <- as.data.frame(REMO.v1.GRCh38)
remo_regions$position <- (remo_regions$start + remo_regions$end) / 2

mean_pairwise_distance <- function(midpoints) {
  if (length(midpoints) < 2) return(NA)
  dists <- combn(midpoints, 2, function(x) abs(diff(x)))
  mean(dists)
}
                 
module_mean_distance <- remo_regions |>
  group_by(REMO) |>
  summarise(mean_distance = mean_pairwise_distance(position)) |>
  ungroup()

mn_dist <- module_mean_distance$mean_distance
names(mn_dist) <- module_mean_distance$REMO
                 
REMO.v1.GRCh38.metadata$mean_distance <- mn_dist[REMO.v1.GRCh38.metadata$REMO]

p1 <- ggplot(REMO.v1.GRCh38.metadata, aes(CREs)) +
  geom_histogram(binwidth=1, fill = REMO_COL) +
  theme_bw() +
  xlab("Number of CREs in module") + ggtitle("Module size")

p2 <- ggplot(REMO.v1.GRCh38.metadata, aes(mean_distance/1000)) +
  geom_histogram(fill = REMO_COL) +
  theme_bw() +
  xlab("Mean distance between CREs (kb)") + ggtitle("Genomic distance")

mean(REMO.v1.GRCh38.metadata$mean_distance, na.rm = TRUE) / 1000
mean(REMO.v1.GRCh38.metadata$CREs, na.rm = TRUE)
sd(REMO.v1.GRCh38.metadata$CREs, na.rm = TRUE)

# coaccessibility within module in single-cell data
ccre <- Read10X('data/pbmc/ccre_multiome/', gene.column=1)
remo_regions$pk <- paste(remo_regions$seqnames, remo_regions$start, remo_regions$end, sep = "-")

# get chr1 ccres with >5 modules
ncre <- REMO.v1.GRCh38.metadata$CREs
names(ncre) <- REMO.v1.GRCh38.metadata$REMO
remo_regions$CREs <- ncre[remo_regions$REMO]
chr1_modules <- remo_regions[remo_regions$seqnames == 'chr1' & remo_regions$CREs > 5, ]
ccre <- ccre[chr1_modules$pk, qc_cells]

# subset to peaks with decent number of counts in PBMC dataset
ccre <- ccre[rowSums(ccre) > 10, ]

unique_modules <- unique(chr1_modules$REMO)

# find Pearson correlation between peaks that are in the same cluster
pearson_within <- vector(mode = "numeric", length = length(unique_modules))
for (i in seq_along(unique_modules)) {
    cluster_name <- unique_modules[i]
    pks_use <- chr1_modules[chr1_modules$REMO == cluster_name, 'pk']
    pks_keep <- intersect(pks_use, rownames(ccre))
    if (length(pks_keep) > 5) {
        d_use <- t(ccre[pks_keep, ])
        r <- cor(as.matrix(d_use))
        diag(r) <- NA
        mean_r <- mean(r, na.rm = TRUE)
        pearson_within[i] <- mean_r
    } else {
        pearson_within[i] <- NA
    }
}

# get correlation between peaks that are in different clusters
pearson_outside <- vector(mode = "numeric", length = length(unique_modules))
for (i in seq_along(unique_modules)) {
    cluster_name <- unique_modules[i]
    pks_use <- chr1_modules[chr1_modules$REMO == cluster_name, 'pk']
    diff_peaks <- sample(chr1_modules[chr1_modules$REMO != cluster_name, 'pk'], 5000)
    pks_keep <- intersect(pks_use, rownames(ccre))
    diff_pks_keep <- intersect(diff_peaks, rownames(ccre))
    if (length(pks_keep) > 5) {
        d_use <- t(ccre[pks_keep, ])
        diff_use <- t(ccre[diff_pks_keep, ])
        r <- cor(as.matrix(d_use), as.matrix(diff_use))
        diag(r) <- NA
        mean_r <- mean(r, na.rm = TRUE)
        pearson_outside[i] <- mean_r
    } else {
        pearson_outside[i] <- NA
    }
}

df <- data.frame(
    "R" = c(pearson_within, pearson_outside),
    "Module" = c(rep("Same", length(pearson_within)), rep("Different", length(pearson_outside)))
)

p3 <- ggplot(df, aes(x = R, color = Module, fill = Module)) +
  geom_density(alpha = 0.5, adjust = 2) +
  theme_bw() +
  xlab("Pearson correlation") +
  ggtitle("CRE coaccessibility") +
  ylab("Density") +
  scale_fill_manual(values = c("#FF9149", "#9FB3DF")) +
  scale_color_manual(values = c("#FF9149", "#9FB3DF")) +
  theme(legend.position = c(.95, .95),
          legend.justification = c("right", "top"),
          legend.box.just = "right",
          legend.margin = margin(6, 6, 6, 6))


# library size per cell
libsize <- data.frame(
    libsize = c(colSums(remo), colSums(peaks), colSums(gex)),
    Method = c(rep("REMO", ncol(remo)), rep("Peaks", ncol(peaks)), rep("GEX", ncol(gex)))
)

libsize$Method <- factor(libsize$Method, levels = c("REMO", "Peaks", "GEX"))

p4 <- ggplot(libsize, aes(libsize, fill = Method)) +
  geom_density(alpha=0.5) + theme_minimal() + scale_x_log10() + xlab("Total cell counts") +
  scale_fill_manual(values = c(REMO_COL, PEAKS_COL, GENE_COL)) +
  theme(legend.position = 'none')

# compare count distributions
countdist_remo <- data.frame(table(sample(remo@x, 1000000, replace = FALSE)))
countdist_remo$Method <- "REMO"

countdist_pk <- data.frame(table(sample(peaks@x, 1000000, replace = FALSE)))
countdist_pk$Method <- "Peaks"

countdist_gex <- data.frame(table(sample(gex@x, 1000000, replace = FALSE)))
countdist_gex$Method <- "GEX"

countdist <- rbind(countdist_remo, countdist_pk, countdist_gex)
countdist$Method <- factor(countdist$Method, levels = c("REMO", "Peaks", "GEX"))
countdist$Var1 <- as.numeric(countdist$Var1)

p5 <- ggplot(countdist, aes(Var1, Freq, color = Method)) +
  geom_line() + geom_point() + theme_bw() + xlim(c(0, 50)) +
  scale_color_manual(values = c(REMO_COL, PEAKS_COL, GENE_COL)) + scale_y_log10() +
  xlab("Count") + ylab("Frequency") + NoLegend() + ggtitle("Matrix count value")

# number of zeros per cell for each modality
get_prop_nonzero <- function(x) {
    nz <- length(x@x)
    total_positions <- as.numeric(ncol(x)) * as.numeric(nrow(x))
    return(nz/total_positions * 100)
}

# distribution of nonzero values
get_mean_total <- function(x) {
    totals <- colSums(x)
    nz_mean <- totals / colSums(x > 0)
    return(list(totals, nz_mean))
}
remo_stat <- get_mean_total(remo)
peak_stat <- get_mean_total(peaks)
gex_stat <- get_mean_total(gex)

nrep <- length(remo_stat[[1]])

df <- data.frame(total_counts = c(remo_stat[[1]], peak_stat[[1]], gex_stat[[1]]),
                 nonzero_mean = c(remo_stat[[2]], peak_stat[[2]], gex_stat[[2]]),
                 Method = c(rep("REMO", nrep), rep("Peaks", nrep), rep("Gene expression", nrep)))

df$Method <- factor(df$Method, levels = c("Gene expression", "REMO", "Peaks"))

p6 <- ggplot(df, aes(x = total_counts, y = nonzero_mean, color = Method)) +
  geom_point(size=0.1) +
  geom_hline(yintercept=1, color = 'grey') +
  theme_classic() +
  scale_color_manual(values = c(GENE_COL, REMO_COL, PEAKS_COL)) +
  scale_x_log10() +
  ylab("Mean of non-zero counts") +
  xlab("Total cell counts") + NoLegend()

# mean-variance relationship
meanvar <- data.frame(
    mean = c(rowMeans(remo), rowMeans(peaks), rowMeans(gex)),
    var  = c(rowVars(remo), rowVars(peaks), rowVars(gex)),
    Method = c(rep("REMO", nrow(remo)), rep("Peaks", nrow(peaks)), rep("Gene expression", nrow(gex)))
)
meanvar$Method <- factor(meanvar$Method, levels = c("Gene expression", "REMO", "Peaks"))

custom_labels <- function(x) {
  ifelse(abs(x) >= 1, 
         as.character(round(x)),
         as.character(x))
}

p7 <- ggplot(meanvar, aes(x = mean, y = var, color = Method)) +
  geom_point(size=0.1, alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
  theme_bw() +
  scale_color_manual(values = c(GENE_COL, REMO_COL, PEAKS_COL)) +
  scale_x_log10(limits = c(0.001, 100), labels = custom_labels) +
  scale_y_log10(limits = c(0.001, 100), labels = custom_labels) + 
  ylab("Variance") + xlab("Mean") + facet_wrap(~Method, nrow=1) + theme(legend.position = 'none')

# save plots
ggsave("plots/fig1/module_size_distance.pdf", p1 | p2, height = 3, width = 6, create.dir = TRUE)
ggsave("plots/fig1/module_distance.pdf", p2, height = 3, width = 3, create.dir = TRUE)
ggsave("plots/fig1/module_coaccess.pdf", p3, height = 3, width = 3, create.dir = TRUE)
ggsave("plots/fig2/ncount_density.pdf", p4, height = 3, width = 3, create.dir = TRUE)
ggsave("plots/fig2/matrix_count_var.pdf", p5 | p6, height = 3, width = 6, create.dir = TRUE)
ggsave("plots/fig2/nonzero_mean.png", p6 + NoLegend(), height = 3, width = 3, dpi=500, create.dir = TRUE)
ggsave("plots/fig2/mean_variance.png", p7, height = 2.5, width = 6, dpi=500, create.dir = TRUE)
