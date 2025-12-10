library(AnnotationHub)
library(Signac)
ah <- AnnotationHub()


ensdb_v98 <- ah[["AH75011"]]
annotations <- GetGRangesFromEnsDb(ensdb = ensdb_v98)
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"

outf <- snakemake@output[[1]]
message(outf)

saveRDS(annotations, outf)