extract_count_matrix <- function(h5, output_dir) {
    counts <- Seurat::Read10X_h5(h5)
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }
    if ("Gene Expression" %in% names(counts)) {
        rna_matrix <- counts[["Gene Expression"]]
        Matrix::writeMM(rna_matrix, file.path(output_dir, "gex.mtx"))
        genes <- rownames(rna_matrix)
        write.table(genes, file.path(output_dir, "genes.tsv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
        cells <- colnames(rna_matrix)
        write.table(cells, file.path(output_dir, "rna_cells.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
    } else {
        stop("Gene Expression matrix not found in H5 file.")
    }
}
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
    stop("Usage: Rscript run_extraction.R <input_h5> <output_dir>")
}

h5_file <- args[1]
output_dir <- args[2]
extract_count_matrix(h5_file, output_dir)
