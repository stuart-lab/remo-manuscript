#!/usr/bin/env Rscript

suppressPackageStartupMessages({
	library(data.table)
	library(argparse)
})

parser <- ArgumentParser()

parser$add_argument("--renv_dir", required = TRUE,
	help = "Path to renv project root")
parser$add_argument("glm_dir",
	help = "Directory containing per-chromosome GLM result TSV files")
parser$add_argument("--outcome", required = TRUE,
	help = "Outcome column value to filter on (e.g. same_remo)")
parser$add_argument("--term", required = TRUE,
	help = "Term column value to filter on (e.g. same_gene)")
parser$add_argument("--output_name", default = "glm_results_all.tsv",
	help = "Output TSV filename")
parser$add_argument("--overwrite", action = "store_true",
	help = "Overwrite output file if exists")

args <- parser$parse_args()

renv::load(args$renv_dir)

glm_dir <- args$glm_dir
outcome_target <- args$outcome
term_target <- args$term
out_fpath <- file.path(glm_dir, args$output_name)
overwrite <- isTRUE(args$overwrite)

if (!dir.exists(glm_dir)) {
	stop(sprintf("Directory does not exist: %s", glm_dir), call. = FALSE)
}

if (file.exists(out_fpath) && !overwrite) {
	message(sprintf("Output exists and --overwrite not set, skipping: %s", out_fpath))
	quit(save = "no", status = 0)
}

tsv_files <- list.files(glm_dir, pattern = "\\.tsv$", full.names = TRUE, recursive = FALSE)
if (length(tsv_files) == 0L) {
	stop(sprintf("No .tsv files found in: %s", glm_dir), call. = FALSE)
}

# -----------------------------
# Find TSV files
# -----------------------------
tsv_files <- list.files(
	glm_dir,
	pattern = "\\.tsv$",
	full.names = TRUE,
	recursive = FALSE
)

if (length(tsv_files) == 0L) {
	stop(sprintf("No .tsv files found in: %s", glm_dir), call. = FALSE)
}

# -----------------------------
# Helpers
# -----------------------------
extract_chr <- function(path) {
	f <- basename(path)
	m <- regexec("_(chr[^._]+)\\.tsv$", f, perl = TRUE)
	g <- regmatches(f, m)[[1L]]
	if (length(g) >= 2L) return(g[[2L]])

	m2 <- regexec("(chr[0-9]+|chrX|chrY)", f, perl = TRUE)
	g2 <- regmatches(f, m2)[[1L]]
	if (length(g2) >= 2L) return(g2[[2L]])

	NA_character_
}

# -----------------------------
# Extract results
# -----------------------------
out_rows <- list()

for (f in tsv_files) {
	chr <- sub("^chr", "", extract_chr(f))

	dt <- tryCatch(
		fread(f),
		error = function(e) NULL
	)
	if (is.null(dt) || nrow(dt) == 0L) next

	req <- c("outcome", "term", "estimate", "p_two_sided")
	missing <- setdiff(req, names(dt))
	if (length(missing) > 0L) {
		stop(
			sprintf(
				"File %s missing required columns: %s",
				f,
				paste(missing, collapse = ", ")
			),
			call. = FALSE
		)
	}

	row <- dt[outcome == outcome_target & term == term_target]

	if (nrow(row) == 0L) {
		out_rows[[length(out_rows) + 1L]] <- data.table(
			chr = chr,
			coeff = NA_real_,
			log_odds = NA_real_,
			p_val = NA_real_
		)
		next
	}

	coeff <- as.numeric(row$estimate[[1L]])
	p_val <- as.numeric(row$p_two_sided[[1L]])

	out_rows[[length(out_rows) + 1L]] <- data.table(
		chr = chr,
		coeff = coeff,
		log_odds = exp(coeff),
		p_val = p_val
	)
}

out_dt <- rbindlist(out_rows, use.names = TRUE, fill = TRUE)

# -----------------------------
# Chromosome ordering: 1–22
# -----------------------------
out_dt[, chr_num := suppressWarnings(as.integer(chr))]
setorder(out_dt, chr_num)
out_dt[, chr_num := NULL]

# -----------------------------
# Write output
# -----------------------------
out_fpath <- file.path(glm_dir, "glm_results_all.tsv")
fwrite(out_dt, out_fpath, sep = "\t", quote = FALSE, na = "NA")

message(sprintf("Wrote: %s", out_fpath))
