library(GenomicRanges)

input_files <- snakemake@input
output_file <- snakemake@output[[1]]

# filter those that exist
existing <- sapply(input_files, file.exists)
if (!all(existing)) {
    warning("Some peaks missing")
    input_files <- input_files[existing]
}

import_peaks <- function(x) {
    df <- read.table(x, sep = "\t")
    gr <- makeGRangesFromDataFrame(df, seqnames.field = "V1", start.field = "V2", end.field = "V3", keep.extra.columns = FALSE)
    return(gr)
}

# load each peak file
all_peaks <- c()
for (i in input_files) {
    all_peaks <- c(all_peaks, import_peaks(i))
}

combined_peaks <- Reduce(all_peaks, f = c)
combined_peaks <- reduce(combined_peaks)

# filter width and blacklist
combined_peaks <- combined_peaks[width(combined_peaks) < 5000 & width(combined_peaks) > 300]
combined_peaks <- keepStandardChromosomes(combined_peaks, pruning.mode = "coarse")
combined_peaks <- subsetByOverlaps(combined_peaks, Signac::blacklist_hg38_unified, invert = TRUE)

# write to file
peak_df <- as.data.frame(combined_peaks)
write.table(peak_df, output_file, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)