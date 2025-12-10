# analyze brain cancer
# combine + annotate celltype + plot umap
library(ggplot2)
library(patchwork)
library(Signac)
library(Seurat)
source("code/utilities.R")

### REMO ###
I2_remo <- Read10X("data/brain_cancer/I2/remo/", gene.column = 1)
I3_remo <- Read10X("data/brain_cancer/I3/remo/", gene.column = 1) 
I4_remo <- Read10X("data/brain_cancer/I4/remo/", gene.column = 1)
I5_remo <- Read10X("data/brain_cancer/I5/remo/", gene.column = 1)
M7_remo <- Read10X("data/brain_cancer/M7/remo/", gene.column = 1)
M8_remo <- Read10X("data/brain_cancer/M8/remo/", gene.column = 1)

colnames(I2_remo) <- paste0("I2_", colnames(I2_remo))
colnames(I3_remo) <- paste0("I3_", colnames(I3_remo))
colnames(I4_remo) <- paste0("I4_", colnames(I4_remo))
colnames(I5_remo) <- paste0("I5_", colnames(I5_remo))
colnames(M7_remo) <- paste0("M7_", colnames(M7_remo))
colnames(M8_remo) <- paste0("M8_", colnames(M8_remo))

combined_remo <- cbind(I2_remo, I3_remo, I4_remo, I5_remo, M7_remo, M8_remo)

# process obj
remo_obj <- process_remo_obj(combined_remo, dims = 1:30, resolution = 0.5)

### Peaks ###
I2_peaks <- Read10X("data/brain_cancer/I2/peaks/", gene.column = 1)
I3_peaks <- Read10X("data/brain_cancer/I3/peaks/", gene.column = 1)
I4_peaks <- Read10X("data/brain_cancer/I4/peaks/", gene.column = 1)
I5_peaks <- Read10X("data/brain_cancer/I5/peaks/", gene.column = 1)
M7_peaks <- Read10X("data/brain_cancer/M7/peaks/", gene.column = 1)
M8_peaks <- Read10X("data/brain_cancer/M8/peaks/", gene.column = 1)

colnames(I2_peaks) <- paste0("I2_", colnames(I2_peaks))
colnames(I3_peaks) <- paste0("I3_", colnames(I3_peaks))
colnames(I4_peaks) <- paste0("I4_", colnames(I4_peaks))
colnames(I5_peaks) <- paste0("I5_", colnames(I5_peaks))
colnames(M7_peaks) <- paste0("M7_", colnames(M7_peaks))
colnames(M8_peaks) <- paste0("M8_", colnames(M8_peaks))

combined_peaks <- cbind(I2_peaks, I3_peaks, I4_peaks, I5_peaks, M7_peaks, M8_peaks)

# process obj
peaks_obj <- process_atac_obj(combined_peaks, dims = 2:30, features = rownames(combined_peaks), resolution = 0.5)

### annotate celltypes ###
brain_cancer <- readRDS(snakemake@input[['source_obj']])
colnames(brain_cancer) <- paste0(colnames(brain_cancer), "-1")
celltypes <- Idents(brain_cancer)

remo_obj$celltype <- celltypes
peaks_obj$celltype <- celltypes

saveRDS(remo_obj, snakemake@output[['remo_obj']])
saveRDS(peaks_obj, snakemake@output[['peaks_obj']])

# dimplots
umap_remo <- DimPlot(remo_obj, group.by = 'celltype', pt.size = 0.1) + ggtitle("") + theme_dimplot(legend=FALSE)
umap_remo2 <- DimPlot(remo_obj, group.by = 'celltype', pt.size = 0.1) + ggtitle("REMO") + theme_dimplot()
umap_peaks <- DimPlot(peaks_obj, group.by = 'celltype', pt.size = 0.1) + ggtitle("Peaks") + theme_dimplot()

ggsave(file = snakemake@output[['umap_peaks']], umap_peaks, width = 5.5, height = 4, dpi = 500)
ggsave(file = snakemake@output[['umap_remo']], umap_remo, width = 5, height = 5, dpi = 500)
ggsave(file = "plots/umap_remo_brain_cancer_presentation.png", umap_remo2, width = 5.5, height = 4, dpi = 500)
