#!/usr/bin/env Rscript

suppressPackageStartupMessages({
	library(argparse)
	library(data.table)
	library(Matrix)
	library(rtracklayer)
	library(GenomicRanges)
	library(Signac)
})

# ----------
# Arguments
# ----------
parser <- ArgumentParser()

parser$add_argument("--renv_dir", required = TRUE,
	help = "Path to renv project root")
parser$add_argument("--sig_pgl_fpath", required = TRUE,
	help = "Significant PGL metadata TSV path (must include columns 'peak' and 'gene')")
parser$add_argument("--output_dir", required = TRUE,
	help = "Output directory for probability analysis outputs")
parser$add_argument("--max_d", default = 1000000,
	help = "Max genomic distance between peak pairs (bp)")
parser$add_argument("--skip_chrs", default = "",
	help = "Comma-separated list of chromosomes to skip (e.g. chr1,chr2)")
parser$add_argument("--only_chrs", default = "",
	help = "Comma-separated list of chromosomes to run (overrides skip_chrs if set)")
parser$add_argument("--overwrite", action = "store_true",
	help = "Overwrite per-chr outputs if they exist")
parser$add_argument("--verbose", action = "store_true",
	help = "Enable verbose logging (prints intermediate tables and diagnostics)")

args <- parser$parse_args()

renv::load(args$renv_dir)

sig_pgl_fpath <- args$sig_pgl_fpath
output_dir <- args$output_dir
max_d <- as.integer(args$max_d)
overwrite <- isTRUE(args$overwrite)
verbose <- isTRUE(args$verbose)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

linkedPeaks_dir <- file.path(output_dir, "linkedPeaks_metadata")
pairedPeaks_dir <- file.path(output_dir, "pairedPeaks_metadata")
glmResults_dir <- file.path(output_dir, "logistic_regression_results")

dir.create(linkedPeaks_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pairedPeaks_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(glmResults_dir, recursive = TRUE, showWarnings = FALSE)

# ----------
# Helpers
# ----------
vcat <- function(...) {
	if (isTRUE(verbose)) message(...)
}

vprint_head <- function(x, n = 6L, label = NULL) {
	if (!isTRUE(verbose)) return(invisible(NULL))
	if (!is.null(label)) message(label)
	print(utils::head(x, n = n))
	invisible(NULL)
}

parse_chr_list <- function(x) {
	x <- trimws(x)
	if (nchar(x) == 0L) return(character(0))
	trimws(unlist(strsplit(x, ",", fixed = TRUE)))
}

get_linkedPeaks_dt <- function(gr, pgl_dt) {
	# ---- Map peaks to linked genes ----
	gr_peak_names <- GRangesToString(gr)

	pgl_small <- pgl_dt[peak %chin% unique(gr_peak_names)]
	vprint_head(pgl_small, label = "Head(pgl_small):")

	if (isTRUE(verbose)) {
		print(all.equal(pgl_dt, pgl_small))
	}

	genes_per_peak <- pgl_small[
		,
		.(linkedGenes_list = list(unique(gene))),
		by = peak
	]

	vprint_head(genes_per_peak, label = "Head(genes_per_peak):")

	setnames(genes_per_peak, "peak", "peak_name")
	setkey(genes_per_peak, peak_name)

	# ---- Build peak metadata dt ----
	if (isTRUE(verbose)) {
		message("Head(gr):")
		print(utils::head(gr))
	}

	dt <- data.table(
		seqname   = as.character(seqnames(gr)),
		start     = start(gr),
		end       = end(gr),
		names     = as.character(gr$name),
		peak_name = gr_peak_names
	)

	vprint_head(dt, label = "Head(dt peak metadata):")

	dt[, linkedGenes_list := vector("list", .N)]
	setkey(dt, peak_name)
	dt[genes_per_peak, linkedGenes_list := i.linkedGenes_list]

	dt[, linkedGenes := vapply(
		linkedGenes_list,
		function(x) if (length(x) == 0L) "" else paste(x, collapse = ","),
		FUN.VALUE = character(1L)
	)]

	setorder(dt, seqname, start, end)
	dt[]
}

get_pairedPeaks_dt <- function(dt, max_d = 100L, with_stats = TRUE) {
	n <- nrow(dt)

	fmt_count <- function(x) {
		format(
			as.numeric(x),
			big.mark   = ",",
			scientific = (abs(x) >= 1e6)
		)
	}

	# ---- Handle n <= 1 ----
	if (n <= 1L) {
		pair_dt <- data.table(
			peak_names = character(),
			log_d      = numeric(),
			same_remo  = logical(),
			same_gene  = logical()
		)

		if (!with_stats) return(pair_dt)

		worst_case_pairs <- as.double(n) * (n - 1) / 2

		stats_dt <- data.table(
			metric    = c("n", "pair_count", "comparison_count_success", "worst_case_pairs", "fraction_of_worst"),
			value     = c(n, 0L, 0L, worst_case_pairs, NA_real_),
			formatted = c(
				fmt_count(n),
				fmt_count(0L),
				fmt_count(0L),
				fmt_count(worst_case_pairs),
				NA_character_
			)
		)

		return(list(pairs = pair_dt, stats = stats_dt))
	}

	# ---- Build pairs with sliding window ----
	capacity <- 1e6L
	idx_i <- integer(capacity)
	idx_j <- integer(capacity)
	k <- 0L

	comparison_count <- 0L
	pair_count <- 0L

	grow <- function() {
		new_cap <- as.integer(capacity * 2L)
		if (new_cap <= capacity) {
			stop("Cannot grow buffer further; capacity overflow.")
		}
		idx_i <<- c(idx_i, integer(new_cap - capacity))
		idx_j <<- c(idx_j, integer(new_cap - capacity))
		capacity <<- new_cap
	}

	for (i in seq_len(n - 1L)) {
		chr_i <- dt$seqname[i]
		end_i <- dt$end[i]
		j <- i + 1L

		while (
			j <= n &&
			dt$seqname[j] == chr_i &&
			dt$start[j] <= end_i + max_d
		) {
			comparison_count <- comparison_count + 1L
			k <- k + 1L
			if (k > capacity) grow()
			idx_i[k] <- i
			idx_j[k] <- j
			pair_count <- pair_count + 1L
			j <- j + 1L
		}
	}

	# ---- Construct output dt ----
	if (k == 0L) {
		pair_dt <- data.table(
			peak_names = character(),
			log_d      = numeric(),
			same_remo  = logical(),
			same_gene  = logical()
		)
	} else {
		idx_i <- idx_i[seq_len(k)]
		idx_j <- idx_j[seq_len(k)]

		peak_names <- paste(dt$peak_name[idx_i], dt$peak_name[idx_j], sep = "_")
		d_raw <- dt$start[idx_j] - dt$end[idx_i]
		log_d <- log(d_raw)

		same_remo <- dt$names[idx_i] == dt$names[idx_j]

		has_gene_overlap <- function(a, b) {
			if (length(a) == 0L || length(b) == 0L) return(FALSE)
			length(intersect(a, b)) > 0L
		}

		same_gene <- mapply(
			has_gene_overlap,
			dt$linkedGenes_list[idx_i],
			dt$linkedGenes_list[idx_j]
		)

		pair_dt <- data.table(
			peak_names = peak_names,
			log_d      = log_d,
			same_remo  = same_remo,
			same_gene  = same_gene
		)
	}

	if (!with_stats) return(pair_dt)

	worst_case_pairs <- as.double(n) * (n - 1) / 2
	frac_of_worst <- if (worst_case_pairs > 0) comparison_count / worst_case_pairs else NA_real_

	stats_dt <- data.table(
		metric    = c("n", "pair_count", "comparison_count_success", "worst_case_pairs", "fraction_of_worst"),
		value     = c(n, pair_count, comparison_count, worst_case_pairs, frac_of_worst),
		formatted = c(
			fmt_count(n),
			fmt_count(pair_count),
			fmt_count(comparison_count),
			fmt_count(worst_case_pairs),
			format(frac_of_worst, digits = 3, scientific = TRUE)
		)
	)

	if (isTRUE(verbose)) {
		vprint_head(pair_dt, label = "Head(pair_dt):")
		vprint_head(stats_dt, label = "stats_dt:")
	}

	list(pairs = pair_dt, stats = stats_dt)
}

p_gt0_scinot_from_z <- function(z, digits = 3) {
	log_p <- pnorm(z, lower.tail = FALSE, log.p = TRUE)
	log10_p <- log_p / log(10)
	exponent <- floor(log10_p)
	mantissa <- 10^(log10_p - exponent)
	sprintf(paste0("%.", digits, "fE%+d"), mantissa, exponent)
}

summarise_glm <- function(fit, model_name, outcome) {
	sm <- tryCatch(
		coef(summary(fit)),
		error = function(e) NULL
	)

	if (is.null(sm) || !is.matrix(sm) || nrow(sm) == 0L) {
		return(data.table(
			term = character(),
			estimate = numeric(),
			std_error = numeric(),
			z_value = numeric(),
			p_two_sided = numeric(),
			p_one_sided_gt0 = numeric(),
			model_name = character(),
			outcome = character()
		))
	}

	dt <- as.data.table(sm, keep.rownames = "term")

	old_cols <- intersect(names(dt), c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))

	rename_map <- c(
		"Estimate"   = "estimate",
		"Std. Error" = "std_error",
		"z value"    = "z_value",
		"Pr(>|z|)"   = "p_two_sided"
	)

	setnames(dt, old = old_cols, new = rename_map[old_cols])

	dt[, p_one_sided_gt0 := NA_real_]
	dt[is.finite(z_value), p_one_sided_gt0 := p_gt0_scinot_from_z(z_value)]

	dt[, `:=`(model_name = model_name, outcome = outcome)]

	setcolorder(dt, c(
		"model_name", "outcome", "term",
		"estimate", "std_error", "z_value",
		"p_two_sided", "p_one_sided_gt0"
	))

	dt[]
}

# ----------
# Load inputs
# ----------
skip_chrs <- parse_chr_list(args$skip_chrs)
only_chrs <- parse_chr_list(args$only_chrs)

sig_links <- fread(sig_pgl_fpath)

rel_levels <- seqlevels(REMO.v1.GRCh38)
if (length(only_chrs) > 0L) {
	rel_levels <- intersect(rel_levels, only_chrs)
} else if (length(skip_chrs) > 0L) {
	rel_levels <- setdiff(rel_levels, skip_chrs)
}

vcat(sprintf("max_d: %d", max_d))
vcat(sprintf("Chromosomes to run: %s", paste(rel_levels, collapse = ",")))

# ----------
# Subset to linked peaks
# ----------
linked_peak_strs <- unique(sig_links$peak)
gr_peaks_strs <- GRangesToString(REMO.v1.GRCh38)
gr_sub <- REMO.v1.GRCh38[gr_peaks_strs %in% linked_peak_strs]

# ----------
# Per-chromosome processing
# ----------
for (chr in rel_levels) {
	vcat(sprintf("Processing %s", chr))

	gr_sub_chr <- gr_sub[seqnames(gr_sub) == chr]
	if (length(gr_sub_chr) == 0L) {
		message(sprintf("Skipping %s: no linked peaks", chr))
		next
	}

	if (isTRUE(verbose)) {
		message(sprintf("n_peaks (linked) in %s: %d", chr, length(gr_sub_chr)))
	}

	dt <- get_linkedPeaks_dt(gr_sub_chr, sig_links)
	res <- get_pairedPeaks_dt(dt, max_d = max_d)

	# ---- Write linked peaks ----
	linkedPeaks_fpath <- file.path(linkedPeaks_dir, paste0("linkedPeaks_metadata_", chr, ".tsv"))
	if (overwrite || !file.exists(linkedPeaks_fpath)) {
		fwrite(dt[, linkedGenes_list := NULL], linkedPeaks_fpath, sep = "\t", quote = FALSE)
	}

	# ---- Write paired peaks ----
	pairedPeaks_fpath <- file.path(pairedPeaks_dir, paste0("pairedPeaks_metadata_", chr, ".tsv.gz"))
	if (overwrite || !file.exists(pairedPeaks_fpath)) {
		fwrite(res$pairs, pairedPeaks_fpath, sep = "\t", quote = FALSE)
	}

	# ---- Fit GLMs ----
	model_dt <- res$pairs[
		is.finite(log_d),
		.(same_remo = as.integer(same_remo), log_d = log_d, same_gene = as.integer(same_gene))
	]

	if (nrow(model_dt) == 0L) {
		message(sprintf("Skipping %s: no finite log_d rows", chr))
		next
	}

	if (isTRUE(verbose)) {
		message(sprintf("n_rows(model_dt) in %s: %d", chr, nrow(model_dt)))
	}

	fit_A <- glm(same_remo ~ log_d + same_gene, data = model_dt, family = binomial)
	fit_B <- glm(same_gene ~ log_d + same_remo, data = model_dt, family = binomial(link = "logit"))

	glm_res_A <- summarise_glm(fit_A, "Model_A", "same_remo")
	glm_res_B <- summarise_glm(fit_B, "Model_B", "same_gene")
	glm_results <- rbind(glm_res_A, glm_res_B, fill = TRUE)

	# ---- Write GLM summaries ----
	glmSummary_fpath <- file.path(glmResults_dir, paste0("glm_summary_", chr, ".txt"))
	summary_A <- capture.output(summary(fit_A))
	summary_B <- capture.output(summary(fit_B))

	writeLines(
		c(
			paste0("===== Model A (", chr, ") ====="),
			summary_A,
			"",
			paste0("===== Model B (", chr, ") ====="),
			summary_B,
			"\n"
		),
		con = glmSummary_fpath,
		sep = "\n",
		useBytes = TRUE
	)

	# ---- Write GLM table ----
	glmResults_fpath <- file.path(glmResults_dir, paste0("glm_results_", chr, ".tsv"))
	fwrite(glm_results, glmResults_fpath, sep = "\t", quote = FALSE)

	message(sprintf("Done: %s", chr))
}
