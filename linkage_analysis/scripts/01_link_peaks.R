#!/usr/bin/env Rscript

# 01_link_peaks.R
#
# Per-chromosome Signac peak–gene linking using a custom LinkPeaks() wrapper.
#
# This script:
#   1) Loads a preprocessed Seurat object
#   2) Loads REMO peak definitions (BED)
#   3) Loads GENCODE transcripts (GTF)
#   4) Runs a modified LinkPeaks per chromosome
#   5) Writes per-chromosome TSV files with columns: peak, gene, score, pvalue


## ---- renv bootstrap (must run before library()) ----
if (file.exists("renv/activate.R")) {
  source("renv/activate.R")
} else {
  stop("renv/activate.R not found. Run this script from the project root.", call. = FALSE)
}

suppressPackageStartupMessages({
  library(argparse)
  library(Signac)
  library(Seurat)
  library(GenomicRanges)
  library(Matrix)
  library(data.table)
  library(rtracklayer)
  library(future)
  library(future.apply)
  library(pbapply)
  library(REMO.v1.GRCh38)
})

# -------------------------- #
#  Argument parsing
# -------------------------- #

parser <- ArgumentParser(
  description = "Per-chromosome Signac peak–gene linking with REMO peaks and GENCODE transcripts"
)

parser$add_argument(
  "--seurat_rds", required = TRUE,
  help = "Path to Seurat object preprocessed for Signac peak-gene linking (.rds)"
)
parser$add_argument(
  "--remo_bed", required = TRUE,
  help = "Path to REMO peak BED file"
)
parser$add_argument(
  "--gencode_gtf", required = TRUE,
  help = "Path to GENCODE GTF file (gz or plain)"
)
parser$add_argument(
  "--output_dir", required = TRUE,
  help = "Directory for per-chromosome link TSVs"
)

parser$add_argument(
  "--peak_assay", default = "peaks",
  help = "Name of the assay containing peak accessibility data [default: peaks]"
)
parser$add_argument(
  "--expression_assay", default = "RNA",
  help = "Name of the assay containing gene expression data [default: RNA]"
)
parser$add_argument(
  "--peak_layer", default = "counts",
  help = "Layer (slot) of the peak assay to use [default: counts]"
)
parser$add_argument(
  "--expression_layer", default = "data",
  help = "Layer (slot) of the expression assay to use [default: data]"
)

parser$add_argument(
  "--distance", type = "double", default = 5e5,
  help = "Window size (+/-) around TSS for candidate peaks [default: 5e5]"
)
parser$add_argument(
  "--min_distance", type = "double", default = 500,
  help = "Minimum distance to TSS to exclude promoter-proximal peaks [default: 500]"
)
parser$add_argument(
  "--min_cells", type = "integer", default = 10,
  help = "Minimum number of cells with non-zero peak/gene counts [default: 10]"
)
parser$add_argument(
  "--n_sample", type = "integer", default = 200,
  help = "Number of background peaks sampled per test peak [default: 200]"
)
parser$add_argument(
  "--pvalue_cutoff", type = "double", default = 1,
  help = "P-value cutoff for keeping links [default: 1 (no filtering)]"
)
parser$add_argument(
  "--score_cutoff", type = "double", default = 0,
  help = "Absolute correlation cutoff for peaks before background z-scoring [default: 0]"
)

parser$add_argument(
  "--workers", type = "integer", default = 1,
  help = "Number of parallel workers (future) [default: 1]"
)

args <- parser$parse_args()

seurat_rds       <- args$seurat_rds
remo_bed         <- args$remo_bed
gencode_gtf      <- args$gencode_gtf
output_dir       <- args$output_dir
peak_assay       <- args$peak_assay
expression_assay <- args$expression_assay
peak_layer       <- args$peak_layer
expression_layer <- args$expression_layer

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

message("Input Seurat RDS:    ", seurat_rds)
message("REMO BED:            ", remo_bed)
message("GENCODE GTF:         ", gencode_gtf)
message("Output directory:    ", output_dir)
message("Peak assay:          ", peak_assay,       " | layer: ", peak_layer)
message("Expression assay:    ", expression_assay, " | layer: ", expression_layer)

# -------------------------- #
#  Optional renv activation
# -------------------------- #
if (requireNamespace("renv", quietly = TRUE)) {
  try(renv::load(), silent = TRUE)
}

# -------------------------- #
#  Parallel setup
# -------------------------- #

workers <- args$workers

if (workers > 1) {
  message("Configuring future with ", workers, " worker(s)")
  if (.Platform$OS.type == "unix") {
    future::plan("multicore", workers = workers)
  } else {
    future::plan("multisession", workers = workers)
  }
  options(future.globals.maxSize = 50 * 1024 ^ 3) # 50 GB
} else {
  message("Using sequential execution (workers = 1)")
  future::plan("sequential")
}

# -------------------------- #
#  Helper: CollapseDuplicateColumns
# -------------------------- #

#' Collapse duplicate columns of a dgCMatrix using max aggregation.
#'
#' @param mat A dgCMatrix with possibly duplicated column names.
#' @return A dgCMatrix with unique column names, where each entry is the max
#'         across all original duplicate columns.
CollapseDuplicateColumns <- function(mat) {
  if (!inherits(mat, "dgCMatrix")) {
    stop("Input must be of class dgCMatrix as sparse matrix operations are used")
  }
  columns <- colnames(mat)
  if (is.null(columns)) {
    stop("Matrix must have column names.")
  }

  u_colnames <- unique(columns)

  # Coerce to triplet form for efficient grouped aggregation
  T <- as(mat, "dgTMatrix")

  # If there are no non-zero entries, we still want to
  # collapse to one column per unique name, but everything is zero.
  if (length(T@x) == 0L) {
    out <- Matrix::sparseMatrix(
      i    = integer(0),
      j    = integer(0),
      x    = numeric(0),
      dims = c(nrow(mat), length(u_colnames))
    )
    rownames(out) <- rownames(mat)
    colnames(out) <- u_colnames
    return(out)
  }

  # Map each non-zero entry's column to its unique-name index
  group <- match(columns[T@j + 1L], u_colnames)

  df <- data.frame(
    i = T@i + 1L,  # 1-based row index
    j = group,     # 1-based unique column index
    x = T@x        # values
  )

  agg <- stats::aggregate(x ~ i + j, data = df, FUN = max)

  out <- Matrix::sparseMatrix(
    i    = agg$i,
    j    = agg$j,
    x    = agg$x,
    dims = c(nrow(mat), length(u_colnames))
  )
  colnames(out) <- u_colnames
  rownames(out) <- rownames(mat)
  out
}

# -------------------------- #
#  Helper: LinkPeaksCustom
# -------------------------- #

LinkPeaksCustom <- function(
  object,
  peak.assay,
  expression.assay,
  peak.slot       = "counts",
  expression.slot = "data",
  method          = "pearson",
  gene.coords     = NULL,
  distance        = 5e+05,
  min.distance    = NULL,
  min.cells       = 10,
  genes.use       = NULL,
  n_sample        = 200,
  pvalue_cutoff   = 0.05,
  score_cutoff    = 0.05,
  gene.id         = FALSE,
  verbose         = TRUE
) {

  if (!inherits(x = object[[peak.assay]], what = "ChromatinAssay")) {
    stop("The requested assay is not a ChromatinAssay")
  }

  if (!is.null(x = min.distance)) {
    if (!is.numeric(x = min.distance)) stop("min.distance must be numeric")
    if (min.distance <= 0) min.distance <- NULL
  }

  features.match <- c("GC.percent", "count", "sequence.length")
  cor_method <- switch(
    method,
    pearson  = corSparse,
    spearman = SparseSpearmanCor,
    stop("method must be 'pearson' or 'spearman'")
  )

  if (is.null(x = gene.coords)) {
    annot <- Annotation(object = object[[peak.assay]])
    if (is.null(x = annot)) stop("Gene annotations not found")
    gene.coords <- CollapseToLongestTranscript(annot)
  }

  meta.features <- GetAssayData(
    object = object, assay = peak.assay, layer = "meta.features"
  )
  if (!(all(c("GC.percent", "sequence.length") %in% colnames(meta.features)))) {
    stop("Run RegionsStats before calling this function.")
  }

  if (!("count" %in% colnames(meta.features))) {
    data.use <- GetAssayData(object = object[[peak.assay]], layer = "counts")
    hvf.info <- FindTopFeatures(object = data.use, verbose = FALSE)
    hvf.info <- hvf.info[rownames(meta.features), , drop = FALSE]
    meta.features <- cbind(meta.features, hvf.info)
  }

  peak.data <- GetAssayData(
    object = object, assay = peak.assay, layer = peak.slot
  )

  if (!(expression.slot %in% Layers(object[[expression.assay]]))) {
    stop("Requested expression layer not found in assay '", expression.assay, "'")
  }

  expression.data <- GetAssayData(
    object = object, assay = expression.assay, layer = expression.slot
  )

  peakcounts <- Matrix::rowSums(peak.data > 0)
  genecounts <- Matrix::rowSums(expression.data > 0)

  peaks.keep <- peakcounts > min.cells
  genes.keep <- genecounts > min.cells
  peak.data  <- peak.data[peaks.keep, ]

  if (!is.null(genes.use)) {
    genes.keep <- intersect(names(genes.keep[genes.keep]), genes.use)
  }

  ## ---- CASE 1: no genes left for this chromosome ----
  if (length(genes.keep) == 0L) {
    if (verbose) message("No genes remain after filtering; skipping.")
    Links(object[[peak.assay]]) <- GenomicRanges::GRanges()
    return(object)
  }

  expression.data <- expression.data[genes.keep, , drop = FALSE]
  genes <- rownames(expression.data)

  if (gene.id) {
    gene.coords.use <- gene.coords[gene.coords$gene_id %in% genes, ]
    gene.coords.use$gene_name <- gene.coords.use$gene_id
  } else {
    gene.coords.use <- gene.coords[gene.coords$gene_name %in% genes, ]
  }

  ## ---- CASE 2: no coordinates for these genes ----
  if (length(gene.coords.use) == 0L) {
    if (verbose) message("No gene coordinates found; skipping chromosome.")
    Links(object[[peak.assay]]) <- GenomicRanges::GRanges()
    return(object)
  }

  if (length(gene.coords.use) < nrow(expression.data)) {
    message("Found coordinates for ", length(gene.coords.use),
            " of ", nrow(expression.data), " genes.")
  }

  peaks <- granges(object[[peak.assay]])[peaks.keep]

  peak_distance_matrix <- Signac:::DistanceToTSS(
    peaks = peaks,
    genes = gene.coords.use,
    distance = distance
  )
  peak_distance_matrix <- CollapseDuplicateColumns(peak_distance_matrix)

  if (!is.null(min.distance)) {
    peak_distance_matrix_min <- Signac:::DistanceToTSS(
      peaks = peaks,
      genes = gene.coords.use,
      distance = min.distance
    )
    peak_distance_matrix_min <- CollapseDuplicateColumns(peak_distance_matrix_min)
    peak_distance_matrix <- peak_distance_matrix - peak_distance_matrix_min
  }

  ## ---- CASE 3: no peaks within distance window ----
  if (sum(peak_distance_matrix) == 0) {
    if (verbose) message("No peaks within distance; skipping chromosome.")
    Links(object[[peak.assay]]) <- GenomicRanges::GRanges()
    return(object)
  }

  if (verbose) {
    message("Testing ", nrow(expression.data), " genes and ",
            sum(Matrix::rowSums(peak_distance_matrix) > 0), " peaks")
  }

  genes.use <- colnames(peak_distance_matrix)
  all.peaks <- rownames(peak.data)

  peak.data <- Matrix::t(peak.data)

  coef.vec <- gene.vec <- zscore.vec <- c()

  if (future::nbrOfWorkers() > 1) {
    mylapply <- future.apply::future_lapply
  } else {
    mylapply <- if (verbose) pbapply::pblapply else lapply
  }

  res <- mylapply(seq_along(genes.use), function(i) {
    peak.use <- as.logical(peak_distance_matrix[, genes.use[[i]]])
    gene.expression <- Matrix::t(expression.data[genes.use[[i]], , drop = FALSE])
    gene.chrom <- as.character(seqnames(gene.coords.use[i]))

    if (verbose) message(sum(peak.use), " peaks found for ", genes.use[[i]])
    if (sum(peak.use) < 2)
      return(list(gene=NULL, coef=NULL, zscore=NULL))

    peak.access <- peak.data[, peak.use, drop = FALSE]
    coef.result <- cor_method(X = peak.access, Y = gene.expression)
    rownames(coef.result) <- colnames(peak.access)
    coef.result <- coef.result[abs(coef.result) > score_cutoff, , drop=FALSE]

    if (nrow(coef.result) == 0)
      return(list(gene=NULL, coef=NULL, zscore=NULL))

    peaks.test <- rownames(coef.result)
    trans.peaks <- all.peaks[!grepl(paste0("^", gene.chrom, "-"), all.peaks)]
    meta.use <- meta.features[trans.peaks, ]
    pk.use <- meta.features[peaks.test, ]

    bg.peaks <- lapply(seq_len(nrow(pk.use)), function(x) {
      MatchRegionStats(
        meta.feature  = meta.use,
        query.feature = pk.use[x, , drop = FALSE],
        features.match = features.match,
        n = n_sample,
        verbose = FALSE
      )
    })

    bg.access <- peak.data[, unlist(bg.peaks), drop = FALSE]
    bg.coef   <- cor_method(X = bg.access, Y = gene.expression)
    rownames(bg.coef) <- colnames(bg.access)

    zscores <- numeric(length(peaks.test))
    for (j in seq_along(peaks.test)) {
      coef.use <- bg.coef[((j-1)*n_sample+1):(j*n_sample), ]
      zscores[j] <- (coef.result[j] - mean(coef.use)) / sd(coef.use)
    }

    names(coef.result) <- peaks.test
    names(zscores)     <- peaks.test

    list(
      gene   = rep(i, length(coef.result)),
      coef   = coef.result,
      zscore = zscores
    )
  })

  ## Combine outputs
  gene.vec   <- unlist(lapply(res, `[[`, 1))
  coef.vec   <- unlist(lapply(res, `[[`, 2))
  zscore.vec <- unlist(lapply(res, `[[`, 3))

  ## ---- CASE 4: no significant links ----
  if (length(coef.vec) == 0) {
    if (verbose) message("No significant links; skipping chromosome.")
    Links(object[[peak.assay]]) <- GenomicRanges::GRanges()
    return(object)
  }

  peak.key <- seq_along(unique(names(coef.vec)))
  names(peak.key) <- unique(names(coef.vec))

  coef.matrix <- Matrix::sparseMatrix(
    i = gene.vec,
    j = peak.key[names(coef.vec)],
    x = coef.vec,
    dims = c(length(genes.use), max(peak.key))
  )
  rownames(coef.matrix) <- genes.use
  colnames(coef.matrix) <- names(peak.key)

  links <- Signac:::LinksToGRanges(coef.matrix, gene.coords = gene.coords.use)

  z.matrix <- Matrix::sparseMatrix(
    i = gene.vec,
    j = peak.key[names(zscore.vec)],
    x = zscore.vec,
    dims = c(length(genes.use), max(peak.key))
  )
  rownames(z.matrix) <- genes.use
  colnames(z.matrix) <- names(peak.key)

  z.lnk <- Signac:::LinksToGRanges(z.matrix, gene.coords = gene.coords.use)

  links$zscore <- z.lnk$score
  links$pvalue <- pnorm(-abs(links$zscore))
  links <- links[links$pvalue < pvalue_cutoff]

  Links(object[[peak.assay]]) <- links
  return(object)
}

# -------------------------- #
#  Main: load data
# -------------------------- #

message("Loading Seurat object...")
data <- readRDS(seurat_rds)

message("Loading REMO peaks (BED)...")
remo_peaks <- import(remo_bed, format = "BED")
rel_levels <- levels(seqnames(remo_peaks))

message("Loading GENCODE annotations (GTF)...")
gencode <- import(gencode_gtf)
gencode_transcripts <- gencode[
  gencode$type == "transcript" &
    seqnames(gencode) %in% rel_levels
]

message("Chromosomes to process: ", paste(rel_levels, collapse = ", "))

DefaultAssay(data) <- peak_assay

# -------------------------- #
#  Per-chromosome linking
# -------------------------- #

for (level in rel_levels) {
  links_outfile <- file.path(output_dir, sprintf("%s.tsv", level))

  message("Processing chromosome: ", level)

  gene_coords_chr <- gencode_transcripts[seqnames(gencode_transcripts) == level]
  gene_set_chr    <- gene_coords_chr$gene_name

  message("  Genes on this chromosome: ", length(gene_set_chr))

  data <- LinkPeaksCustom(
    object           = data,
    peak.assay       = peak_assay,
    expression.assay = expression_assay,
    peak.slot        = peak_layer,
    expression.slot  = expression_layer,
    gene.coords      = gene_coords_chr,
    genes.use        = gene_set_chr,
    distance         = args$distance,
    min.distance     = args$min_distance,
    min.cells        = args$min_cells,
    n_sample         = args$n_sample,
    pvalue_cutoff    = args$pvalue_cutoff,
    score_cutoff     = args$score_cutoff,
    verbose          = TRUE
  )

  links_gr <- Links(data[[peak_assay]])
  links_df <- if (length(links_gr) > 0) {
    as.data.frame(links_gr)[, c("peak", "gene", "score", "pvalue")]
  } else {
    data.frame(peak = character(), gene = character(),
               score = numeric(), pvalue = numeric())
  }

  data.table::fwrite(
    links_df,
    file      = links_outfile,
    sep       = "\t",
    quote     = FALSE
  )

  message("  Wrote links to: ", links_outfile, " (", nrow(links_df), " rows)")
}

message("All chromosomes processed.")
