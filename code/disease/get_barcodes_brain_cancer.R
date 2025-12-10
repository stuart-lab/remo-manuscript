library(Seurat)
library(Signac)
options(scipen = 999)

brain_cancer <- readRDS(snakemake@input[['source_obj']])

# get peak set
bed_data <- data.frame(do.call(rbind, strsplit(rownames(brain_cancer), "-")))
colnames(bed_data) <- c("chr", "start", "end")
bed_data$start <- as.numeric(bed_data$start)
bed_data$end <- as.numeric(bed_data$end)

write.table(bed_data, file = snakemake@output[['peaks_bed']], sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# get cell barcodes
barcodes <- colnames(brain_cancer)
sample_name <- brain_cancer$biological_replicates

df <- data.frame(barcode = barcodes, sample = sample_name, stringsAsFactors = FALSE)
samples <- sort(unique(sample_name))

for (sample in samples) {
    barcodes_for_sample <- df$barcode[df$sample == sample]

    cleaned_barcodes <- sub("^[^_]+_", "", barcodes_for_sample)
    cleaned_barcodes <- paste0(cleaned_barcodes, "-1")
    
    out_path <- paste0("data/brain_cancer/", sample, "/", sample, "_barcodes.txt")
    write.table(cleaned_barcodes, file = out_path, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}