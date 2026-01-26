suppressPackageStartupMessages({
  library(argparse)
  library(Signac)
  library(Seurat)
  library(GenomicRanges)
  library(Matrix)
  library(data.table)
  library(rtracklayer)
})

# -----------------------------
# Args
# -----------------------------
parser <- ArgumentParser()

parser$add_argument("--links_dir_csv", required = TRUE,
  help = "CSV with column: links_dir")
parser$add_argument("--filter_conds_csv", required = TRUE,
  help = "CSV with columns: condition,pvalue,score")
parser$add_argument("--combined_outdir", required = TRUE,
  help = "Output dir for combined per-chr metadata links")
parser$add_argument("--filtered_parent_dir", required = TRUE,
  help = "Parent dir for per-condition outputs")
parser$add_argument("--gencode_gtf", required = TRUE,
  help = "GENCODE GTF")
parser$add_argument("--remo_bed", required = TRUE,
  help = "REMO BED")
parser$add_argument("--distance_to_tss", default = 500,
  help = "Distance used in Signac:::DistanceToTSS (proximal genes)")
parser$add_argument("--overwrite", action = "store_true",
  help = "Overwrite per-chr combined + per-condition outputs if they exist")

args <- parser$parse_args()

links_dir_csv <- args$links_dir_csv
filter_conds_csv <- args$filter_conds_csv
combined_outdir <- args$combined_outdir
filtered_parent_dir <- args$filtered_parent_dir
gencode_gtf <- args$gencode_gtf
remo_bed <- args$remo_bed
distance_to_tss <- as.numeric(args$distance_to_tss)
overwrite <- isTRUE(args$overwrite)

dir.create(combined_outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(filtered_parent_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Helpers
# -----------------------------
CollapseDuplicateColumns <- function(mat) {
  if (!inherits(mat, "dgCMatrix")) stop("Input must be dgCMatrix.")
  columns <- colnames(mat)
  if (is.null(columns)) stop("Matrix must have column names.")
  u_colnames <- unique(columns)

  T <- as(mat, "TsparseMatrix")
  if (length(T@x) == 0L) {
    out <- Matrix::sparseMatrix(
      i = integer(0), j = integer(0), x = numeric(0),
      dims = c(nrow(mat), length(u_colnames))
    )
    rownames(out) <- rownames(mat)
    colnames(out) <- u_colnames
    return(out)
  }

  group <- match(columns[T@j + 1L], u_colnames)
  df <- data.frame(i = T@i + 1L, j = group, x = T@x)
  agg <- stats::aggregate(x ~ i + j, data = df, FUN = max)

  out <- Matrix::sparseMatrix(
    i = agg$i, j = agg$j, x = agg$x,
    dims = c(nrow(mat), length(u_colnames))
  )
  colnames(out) <- u_colnames
  rownames(out) <- rownames(mat)
  out
}

collapse_vec <- function(x) {
  x <- as.character(unlist(x))
  if (length(x) == 0L) "" else paste(x, collapse = ",")
}

# -----------------------------
# Load configs
# -----------------------------
dt_obj_links <- fread(links_dir_csv)
stopifnot("links_dir" %in% names(dt_obj_links))
link_dirs <- dt_obj_links$links_dir

config <- fread(filter_conds_csv)
stopifnot(all(c("condition", "pvalue", "score") %in% names(config)))
conds <- config$condition

# -----------------------------
# Load REMO + GENCODE, compute pCRE
# -----------------------------
summary_fpath <- file.path(combined_outdir, "summary.txt")
if (overwrite && file.exists(summary_fpath)) file.remove(summary_fpath)

remo <- import(remo_bed, format = "BED")
start(remo) <- start(remo) - 1
rel_levels <- levels(seqnames(remo))

peak_strs <- GRangesToString(remo)
peak_grps <- remo$name
names(peak_grps) <- peak_strs

gencode <- import(gencode_gtf)
gencode_transcripts <- gencode[gencode$type == "transcript" & seqnames(gencode) %in% rel_levels]

dt_gc <- as.data.table(gencode_transcripts)[, .(gene_name, seqnames = as.character(seqnames))]
multi_chr_genes <- dt_gc[, .(n_chr = uniqueN(seqnames)), by = gene_name][n_chr > 1]

pdm <- Signac:::DistanceToTSS(peaks = remo, genes = gencode_transcripts, distance = distance_to_tss)
pdm <- CollapseDuplicateColumns(pdm)

T <- as(pdm, "TsparseMatrix")
dt_long <- data.table(
  peak = rownames(pdm)[T@i + 1L],
  gene = colnames(pdm)[T@j + 1L]
)
dt_hits <- dt_long[, .(genes = list(unique(gene))), by = peak]

pcre_by_peak <- dt_hits[, .(proximal_genes = list(unique(genes))), by = peak]
setkey(pcre_by_peak, peak)

# -----------------------------
# Combine PGLs per chr across datasets + add metadata
# -----------------------------
chr_links_list <- list()

for (chr in rel_levels) {
  links_list <- vector("list", length(link_dirs))
  for (i in seq_along(link_dirs)) {
    link_dir <- link_dirs[[i]]
    link_fpath <- file.path(link_dir, paste0(chr, ".tsv"))
    links <- fread(link_fpath)

    links <- links[!(gene %in% multi_chr_genes$gene_name), ]
    links[, remo_name := peak_grps[peak]]
    links[, dataset := basename(link_dir)]

    links_list[[i]] <- links
  }

  chr_links <- rbindlist(links_list, use.names = TRUE, fill = TRUE)[order(peak, gene)]
  chr_links_list[[chr]] <- chr_links

  chr_all_peaks <- remo[seqnames(remo) == chr]
  cat(sprintf("%s: %d linked peaks out of %d peaks\n",
              chr, length(unique(chr_links$peak)), length(chr_all_peaks)),
      file = summary_fpath, append = TRUE)
  cat(sprintf("%s: %d linked modules out of %d modules\n\n",
              chr, length(unique(chr_links$remo_name)), length(unique(chr_all_peaks$name))),
      file = summary_fpath, append = TRUE)

  outlink_fpath <- file.path(combined_outdir, paste0(chr, ".tsv"))
  if (overwrite || !file.exists(outlink_fpath)) {
    write.table(chr_links, outlink_fpath, sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

# -----------------------------
# For each condition: filter + collect genes_per_peak
# -----------------------------
for (cond in conds) {
  cnd <- config[condition == cond]
  p_value <- cnd$pvalue
  corr_coeff <- cnd$score
  if (is.null(p_value) || is.na(p_value)) p_value <- 1
  if (is.null(corr_coeff) || is.na(corr_coeff)) corr_coeff <- 0

  output_dir <- file.path(filtered_parent_dir, cond)
  pgl_metadata_dir <- file.path(output_dir, "pgl_metadata")
  genes_per_peak_dir <- file.path(output_dir, "genes_per_peak")
  dir.create(pgl_metadata_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(genes_per_peak_dir, recursive = TRUE, showWarnings = FALSE)

  all_pgl_metadata_list <- list()
  all_genesPerPeak_list <- list()

  for (chr in rel_levels) {
    links <- chr_links_list[[chr]]
    sig_links <- links[pvalue < p_value & abs(score) >= corr_coeff]

    chr_pgl_fpath <- file.path(pgl_metadata_dir, paste0(chr, ".tsv"))
    if (overwrite || !file.exists(chr_pgl_fpath)) {
      write.table(sig_links, chr_pgl_fpath, sep = "\t", quote = FALSE, row.names = FALSE)
    }
    all_pgl_metadata_list[[chr]] <- sig_links

    chr_peak_strs <- GRangesToString(remo[seqnames(remo) == chr])
    chr_remos <- remo[seqnames(remo) == chr]$name
    peaks_dt <- data.table(peak = chr_peak_strs, remo_name = chr_remos)

    chr_genesPerPeak <- pcre_by_peak[peaks_dt, on = "peak"]

    dcre_by_peak <- sig_links[
      peak %in% chr_peak_strs,
      .(distal_genes = list(unique(gene))),
      by = peak
    ]
    setkey(dcre_by_peak, peak)
    chr_genesPerPeak <- dcre_by_peak[chr_genesPerPeak, on = "peak"]

    chr_genesPerPeak[, all_genes := Map(function(p, d) unique(c(p, d)), proximal_genes, distal_genes)]

    cols <- c("proximal_genes", "distal_genes", "all_genes")
    chr_genesPerPeak[, (cols) := lapply(.SD, function(col) vapply(col, collapse_vec, character(1L))),
                     .SDcols = cols]

    setcolorder(chr_genesPerPeak, c("peak", "remo_name", "proximal_genes", "distal_genes", "all_genes"))

    chr_gpp_fpath <- file.path(genes_per_peak_dir, paste0(chr, ".tsv"))
    if (overwrite || !file.exists(chr_gpp_fpath)) {
      write.table(chr_genesPerPeak, chr_gpp_fpath, sep = "\t", quote = FALSE, row.names = FALSE)
    }
    all_genesPerPeak_list[[chr]] <- chr_genesPerPeak
  }

  all_flinks <- rbindlist(all_pgl_metadata_list, use.names = TRUE, fill = TRUE)
  all_flinks_fpath <- file.path(pgl_metadata_dir, "pgls_metadata_all_chr.tsv")
  if (overwrite || !file.exists(all_flinks_fpath)) {
    write.table(all_flinks, all_flinks_fpath, sep = "\t", quote = FALSE, row.names = FALSE)
  }

  all_gpp <- rbindlist(all_genesPerPeak_list, use.names = TRUE, fill = TRUE)
  all_gpp_fpath <- file.path(genes_per_peak_dir, "genes_per_peak_all_chr.tsv")
  if (overwrite || !file.exists(all_gpp_fpath)) {
    write.table(all_gpp, all_gpp_fpath, sep = "\t", quote = FALSE, row.names = FALSE)
  }
}
