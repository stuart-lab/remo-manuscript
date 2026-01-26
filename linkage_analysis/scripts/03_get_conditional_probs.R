#!/usr/bin/env Rscript

suppressPackageStartupMessages({
	library(argparse)
	library(data.table)
	library(rtracklayer)
	library(GenomicRanges)
	library(Signac)
})

parser <- ArgumentParser()

parser$add_argument("--renv_dir", required = TRUE,
	help = "Path to renv project root")
parser$add_argument("--remo_bed", required = TRUE,
	help = "REMO BED(.gz) path")
parser$add_argument("--sig_pgl_fpath", required = TRUE,
	help = "Significant PGL metadata TSV path (must include column 'peak')")
parser$add_argument("--output_dir", required = TRUE,
	help = "Output directory")
parser$add_argument("--output_name", default = "conditional_prob_unlinked_singleCRE.tsv",
	help = "Output TSV filename")
parser$add_argument("--overwrite", action = "store_true",
	help = "Overwrite output file if exists")

args <- parser$parse_args()

renv::load(args$renv_dir)

remo_bed <- args$remo_bed
sig_pgl_fpath <- args$sig_pgl_fpath
output_dir <- args$output_dir
output_fpath <- file.path(output_dir, args$output_name)
overwrite <- isTRUE(args$overwrite)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (file.exists(output_fpath) && !overwrite) {
	message(sprintf("Output exists and --overwrite not set, skipping: %s", output_fpath))
	quit(save = "no", status = 0)
}

remo <- import(remo_bed, format = "BED")
start(remo) <- start(remo) - 1

sig_links <- fread(sig_pgl_fpath)

gr_dt <- as.data.table(remo)
gr_dt[, peak := GRangesToString(remo)]
gr_dt[, remo_name := remo$name]

linked_peaks <- unique(sig_links$peak)
gr_dt[, linked := peak %in% linked_peaks]

peaks_per_remo <- gr_dt[, .N, by = remo_name]
setnames(peaks_per_remo, "N", "n_peaks")

gr_dt <- merge(gr_dt, peaks_per_remo, by = "remo_name", all.x = TRUE)
gr_dt[, singleCRE := (n_peaks == 1L)]

freq <- gr_dt[, .N, by = .(linked, singleCRE)]

N_ls <- freq[linked == TRUE  & singleCRE == TRUE,  N]
N_lm <- freq[linked == TRUE  & singleCRE == FALSE, N]
N_us <- freq[linked == FALSE & singleCRE == TRUE,  N]
N_um <- freq[linked == FALSE & singleCRE == FALSE, N]

N_ls <- ifelse(length(N_ls) == 0L, 0L, N_ls)
N_lm <- ifelse(length(N_lm) == 0L, 0L, N_lm)
N_us <- ifelse(length(N_us) == 0L, 0L, N_us)
N_um <- ifelse(length(N_um) == 0L, 0L, N_um)

N_total <- N_ls + N_lm + N_us + N_um
N_linked <- N_ls + N_lm
N_unlinked <- N_us + N_um
N_singleCRE <- N_ls + N_us
N_multiCRE <- N_lm + N_um

P_unlinked <- N_unlinked / N_total
P_unlinked_given_singleCRE <- N_us / N_singleCRE
P_singleCRE_given_unlinked <- N_us / N_unlinked
P_linked_given_singleCRE <- N_ls / N_singleCRE
P_singleCRE_given_linked <- N_ls / N_linked
P_singleCRE <- N_singleCRE / N_total

cond_summary <- data.table(
	metric = c(
		"P(unlinked)",
		"P(unlinked | singleCRE)",
		"P(singleCRE | unlinked)",
		"P(linked | singleCRE)",
		"P(singleCRE | linked)",
		"P(singleCRE)"
	),
	value = c(
		P_unlinked,
		P_unlinked_given_singleCRE,
		P_singleCRE_given_unlinked,
		P_linked_given_singleCRE,
		P_singleCRE_given_linked,
		P_singleCRE
	)
)

print(cond_summary)
fwrite(cond_summary, output_fpath, sep = "\t", quote = FALSE)
message(sprintf("Wrote: %s", output_fpath))
