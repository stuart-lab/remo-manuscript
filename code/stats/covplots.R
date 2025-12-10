library(Seurat)
library(Signac)
library(ggplot2)
library(GenomicRanges)

obj <- readRDS("../../../signac/vignette_data/pbmc_multiomic.rds")
ccre <- read.table("../../data/REMOv1_GRCh38.bed.gz", sep = "\t", header = FALSE)
colnames(ccre) <- c("chr", "start", "end", "REMO")

ccre <- makeGRangesFromDataFrame(ccre, seqnames.field = "chr", start.field = "start", end.field = "end", keep.extra.columns = TRUE)
ccre <- sort(ccre)

Idents(obj) <- "predicted.id"

ct.freq <- table(Idents(obj))
ct.use <- names(ct.freq)[ct.freq > 200]

obj <- SortIdents(obj)

cells.plot <- c("CD14 Mono", "CD16 Mono", "CD8 Naive", "CD4 Naive", "NK")

region = "chr19-54769391-54779916"
region.highlight = c(GRanges(seqnames = "chr19",ranges = IRanges(start = 54769695, end = 54770170)),
                   GRanges(seqnames = "chr19",ranges = IRanges(start = 54770390, end = 54770889)),
                   GRanges(seqnames = "chr19",ranges = IRanges(start = 54770891, end = 54771916)),
                   GRanges(seqnames = "chr19",ranges = IRanges(start = 54778981, end = 54779228)))

p1 <- CoveragePlot(obj, region = region, idents = cells.plot, region.highlight = region.highlight, peaks = FALSE, ranges = ccre, ranges.group.by = "REMO", annotation = FALSE, links = FALSE)

region = "chr8-92511123-92512247"
extend.upstream = 600
extend.downstream = 2900
region.highlight = c(GRanges(seqnames = "chr8",ranges = IRanges(start = 92511123, end = 92512247)),
                   GRanges(seqnames = "chr8",ranges = IRanges(start = 92513325, end = 92514056)),
                    GRanges(seqnames = "chr8",ranges = IRanges(start = 92514118, end = 92514914)))

p2 <- CoveragePlot(obj, region = region, idents = cells.plot, region.highlight = region.highlight, peaks = FALSE, ranges = ccre, extend.downstream = extend.downstream, extend.upstream = extend.upstream,  ranges.group.by = "REMO", annotation = FALSE, links = FALSE)

region = "chr8-23538452-23542267"
extend.upstream = 1000
extend.downstream = 1000
region.highlight = c(GRanges(seqnames = "chr8",ranges = IRanges(start = 23538452, end = 23539651)),
                   GRanges(seqnames = "chr8",ranges = IRanges(start = 23540417, end = 23542267)))

p3 <- CoveragePlot(obj, region = region, idents = cells.plot, region.highlight = region.highlight, peaks = FALSE, ranges = ccre, extend.downstream = extend.downstream, extend.upstream = extend.upstream,  ranges.group.by = "REMO", annotation = FALSE, links = FALSE)

remo <- readRDS("../../objects/pbmc_multiome_remo.rds")
labels <- readRDS("../../data/pbmc/multiome_labels.rds")
remo <- AddMetaData(remo, labels)

Idents(remo) <- "predicted.id"

mk <- FindMarkers(remo, ident.1 = "CD14+ Monocytes", only.pos = TRUE)
avg <- AverageExpression(remo, features = rownames(mk), layer = 'data', assay = 'RNA')[[1]]
mk$avg <- avg[rownames(mk), 'CD14+ Monocytes']
mk <- mk[mk$p_val_adj < 0.01 & mk$avg_log2FC > 1, ]
mk <- mk[order(mk$avg_log2FC, decreasing = TRUE), ]

FeaturePlot(remo, rownames(mk)[1], order = TRUE)

gr <- ccre[ccre$REMO == rownames(mk)[6]]

p4 <- CoveragePlot(obj, region = gr[2], idents = cells.plot, region.highlight = gr[2], peaks = FALSE, ranges = ccre, extend.upstream = 1000, extend.downstream = 1000,  ranges.group.by = "REMO", annotation = FALSE, links = FALSE)

p5 <- CoveragePlot(obj, region = gr[length(gr)], idents = cells.plot, region.highlight = gr[length(gr)], peaks = FALSE, ranges = ccre, extend.upstream = 1000, extend.downstream = 1000,  ranges.group.by = "REMO", annotation = FALSE, links = FALSE)

start(gr[length(gr)]) - end(gr[2])

ggsave("../../plots/covplot1.pdf", p1, height = 3, width = 6)
ggsave("../../plots/covplot2.pdf", p2, height = 3, width = 6)
ggsave("../../plots/covplot3.pdf", p3, height = 3, width = 6)
ggsave("../../plots/covplot4.pdf", p4, height = 3, width = 4)
ggsave("../../plots/covplot5.pdf", p5, height = 3, width = 4)


# find two peaks that are markers for different cell types that are within 10 kb of each other
cre_obj <- readRDS("../../objects/pbmc_multiome_ccre.rds")
mk0 <- FindMarkers(cre_obj, ident.1 = 0, ident.2 = 3, only.pos = TRUE)
mk4 <- FindMarkers(cre_obj, ident.1 = 3, ident.2 = 0, only.pos = TRUE)

mk0 <- mk0[mk0$avg_log2FC > 2 & mk0$p_val_adj < 0.01, ]
mk4 <- mk4[mk4$avg_log2FC > 2 & mk4$p_val_adj < 0.01, ]

# for each peak, find closest in other list
mk0_gr <- StringToGRanges(rownames(mk0))
mk4_gr <- StringToGRanges(rownames(mk4))

dist_04 <- distanceToNearest(mk0_gr, mk4_gr)
mk0$distance <- mcols(dist_04)$distance
mk0$nearest <- rownames(mk4)[subjectHits(dist_04)]

mk0 <- mk0[order(mk0$distance, decreasing = FALSE), ]
head(mk0, 10)

CoveragePlot(obj, "chr19-45676100-45676541", extend.upstream = 800, extend.downstream = 2000, idents = ct.use, ranges = ccre, ranges.group.by = "REMO", peaks = FALSE, links = FALSE)

cov1 <- CoveragePlot(obj, "chr9-88923059-88924425", extend.upstream = 1200, extend.downstream = 500, idents = ct.use, ranges = ccre, ranges.group.by = "REMO", peaks = FALSE, links = FALSE, annotation = FALSE)

remo_regions <- ccre$REMO
names(remo_regions) <- GRangesToString(ccre)

mk0$REMO <- remo_regions[rownames(mk0)]

remo_annot <- readr::read_tsv("../../data/REMOv1_GRCh38_annotations.tsv.gz")

mod_size <- remo_annot$CREs
names(mod_size) <- remo_annot$Module

mk0$size <- mod_size[mk0$REMO]

head(mk0)

regions <- ccre[ccre$REMO == "REMOv1.GRCh38-191600"]

cov2 <- CoveragePlot(obj, regions[25], extend.upstream = 300, extend.downstream = 300, idents = ct.use, ranges = ccre, ranges.group.by = "REMO", peaks = FALSE, links = FALSE, annotation = FALSE)

start(regions[25]) - 88924425

ggsave("../../plots/covplot_example_1.png", cov1, height = 8, width = 8, dpi=500)
ggsave("../../plots/covplot_example_2.png", cov2, height = 8, width = 6, dpi=500)
