library(ggplot2)
library(patchwork)
library(Signac)
library(Seurat)
source("code/utilities.R")

tissue_vector <- unlist(snakemake@params[['tissues']])

umap_peaks_list <- list()
umap_remo_list <- list()
umap_gex_list <- list()

fig_titles <- c("Brain", "PBMC", "Jejunum", "Heart", "Left colon", "Liver",
                "Lung", "Muscle", "Pancreas", "Bile duct", "Fallopian tube", "Grey matter",
                "Kidney", "Pituitary gland", "Placenta", "Retina", "B lymphoma")

for (i in 1:length(tissue_vector)) {
    tissue <- tissue_vector[i]  
    print(tissue)
    print(fig_titles[i])

    remo <- readRDS(paste0("objects/", tissue, "_multiome_remo.rds"))
    peaks <- readRDS(paste0("objects/", tissue, "_multiome_peaks.rds"))
    gex <- readRDS(paste0("objects/", tissue, "_multiome_gex.rds"))

    # Generate UMAP plots
    umap_peaks <- DimPlot(peaks, group.by = "cluster", label = TRUE, label.size = 2.5, pt.size=0.1) +
      theme_dimplot(legend=FALSE) +
      ggtitle(fig_titles[i])

    umap_remo <- DimPlot(remo, group.by = "cluster", label = TRUE, label.size = 2.5, pt.size=0.1) +
      theme_dimplot(legend=FALSE) +
      ggtitle(fig_titles[i])

    umap_gex <- DimPlot(gex, label = TRUE, label.size = 2.5, pt.size=0.1) +
      theme_dimplot(legend=FALSE) +
      ggtitle(fig_titles[i])

    # save plots to list
    umap_peaks_list[[tissue]] <- umap_peaks
    umap_remo_list[[tissue]] <- umap_remo
    umap_gex_list[[tissue]] <- umap_gex

    # PBMC UMAP
    if (tissue == "pbmc") {
        annot <- readRDS("data/pbmc/multiome_labels.rds")
        remo <- AddMetaData(remo, annot, col.name = "predicted.id")
        peaks <- AddMetaData(peaks, annot, col.name = "predicted.id")
        gex <- AddMetaData(gex, annot, col.name = "predicted.id")
        
        umap_peaks <- DimPlot(peaks, group.by = "predicted.id", label = TRUE, repel = TRUE) +
          theme_dimplot(legend=FALSE) +
          ggtitle("Peaks", subtitle = paste0(format(nrow(peaks), big.mark=","), " peaks"))
        
        umap_remo <- DimPlot(remo, group.by = "predicted.id", label = TRUE, repel = TRUE) +
          theme_dimplot(legend=FALSE) +
          ggtitle("REMO", subtitle = "20,000 variable modules")

        umap_gex <- DimPlot(gex, group.by = "predicted.id", label = TRUE, repel = TRUE) +
          theme_dimplot(legend=FALSE) +
          ggtitle("Gene expression", subtitle = "2,000 variable genes")

        ggsave(file = snakemake@output[['pbmc_umap_peaks']], umap_peaks, width = 5, height = 5, dpi = 500)
        ggsave(file = snakemake@output[['pbmc_umap_remo']], umap_remo, width = 5, height = 5, dpi = 500)
        ggsave(file = snakemake@output[['pbmc_umap_gex']], umap_gex, width = 5, height = 5, dpi = 500)
    }
}

# combine umaps 
umap_list_1 <- list(
    umap_gex_list[['brain']],
    umap_gex_list[['pbmc']],
    umap_gex_list[['jejunum']],
    umap_gex_list[['heart']],
    umap_gex_list[['left_colon']],
    umap_gex_list[['liver']],
    umap_gex_list[['muscle']],
    umap_gex_list[['lung']],
    umap_remo_list[['brain']],
    umap_remo_list[['pbmc']],
    umap_remo_list[['jejunum']],
    umap_remo_list[['heart']],
    umap_remo_list[['left_colon']],
    umap_remo_list[['liver']],
    umap_remo_list[['muscle']],
    umap_remo_list[['lung']],
    umap_peaks_list[['brain']],
    umap_peaks_list[['pbmc']],
    umap_peaks_list[['jejunum']],
    umap_peaks_list[['heart']],
    umap_peaks_list[['left_colon']],
    umap_peaks_list[['liver']],
    umap_peaks_list[['muscle']],
    umap_peaks_list[['lung']]
)

umap_list_2 <- list(
    umap_gex_list[['pancreas']],
    umap_gex_list[['bile_duct']],
    umap_gex_list[['fallopian_tube']],
    umap_gex_list[['grey_matter']],
    umap_gex_list[['kidney_cortex']],
    umap_gex_list[['pituitary']],
    umap_gex_list[['placenta']],
    umap_gex_list[['retina']],
    umap_remo_list[['pancreas']],
    umap_remo_list[['bile_duct']],
    umap_remo_list[['fallopian_tube']],
    umap_remo_list[['grey_matter']],
    umap_remo_list[['kidney_cortex']],
    umap_remo_list[['pituitary']],
    umap_remo_list[['placenta']],
    umap_remo_list[['retina']],
    umap_peaks_list[['pancreas']],
    umap_peaks_list[['bile_duct']],
    umap_peaks_list[['fallopian_tube']],
    umap_peaks_list[['grey_matter']],
    umap_peaks_list[['kidney_cortex']],
    umap_peaks_list[['pituitary']],
    umap_peaks_list[['placenta']],
    umap_peaks_list[['retina']]
)

# make row labels
row_labels <- list(
    ggplot() + 
        annotate("text", x = 1, y = 1, label = "RNA", angle = 90, size = 10, hjust = 0.5, vjust = 0.5) + 
        theme_void(),
    ggplot() + 
        annotate("text", x = 1, y = 1, label = "REMO", angle = 90, size = 10, hjust = 0.5, vjust = 0.5) + 
        theme_void(),
    ggplot() + 
        annotate("text", x = 1, y = 1, label = "Peaks", angle = 90, size = 10, hjust = 0.5, vjust = 0.5) + 
        theme_void()
)

# dimplot 1 (brain, fetal_cerebellum, fetal_brain, pbmc, jejunum, heart, left colon, liver, lung)
umap_rows <- list(
    wrap_plots(umap_list_1[1:8], ncol = 8),
    wrap_plots(umap_list_1[9:16], ncol = 8),
    wrap_plots(umap_list_1[17:24], ncol = 8)
)
dimplots_1 <- wrap_plots(
    row_labels[[1]], umap_rows[[1]],
    row_labels[[2]], umap_rows[[2]],
    row_labels[[3]], umap_rows[[3]],
    ncol = 2, widths = c(0.5, 8)
)

# dimplot 2 (muscle, pancreas, bile duct, fallopian tube, grey matter, kidney, pituitary, placenta, retina)
umap_rows <- list(
    wrap_plots(umap_list_2[1:8], ncol = 8),
    wrap_plots(umap_list_2[9:16], ncol = 8),
    wrap_plots(umap_list_2[17:24], ncol = 8)
)
dimplots_2 <- wrap_plots(
    row_labels[[1]], umap_rows[[1]],
    row_labels[[2]], umap_rows[[2]],
    row_labels[[3]], umap_rows[[3]],
    ncol = 2, widths = c(0.5, 8)
)

ggsave(file = snakemake@output[['umap_1']], dimplots_1, width = 24, height = 9, dpi = 500)
ggsave(file = snakemake@output[['umap_2']], dimplots_2, width = 24, height = 9, dpi = 500)