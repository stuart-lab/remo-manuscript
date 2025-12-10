library(Signac)
library(Seurat)
library(BSgenome.Hsapiens.UCSC.hg38)

gex <- readRDS(snakemake@input[['gex']])
annot <- readRDS(snakemake@input[['annotation']])
ccre <-  Read10X(snakemake@input[['ccre']], gene.column=1)
frags <- snakemake@input[['frags']]

gex[['peaks']] <- CreateChromatinAssay(counts = ccre[, colnames(gex)], annotation = annot, fragments = frags)
DefaultAssay(gex) <- "peaks"
gex <- RegionStats(gex, genome = BSgenome.Hsapiens.UCSC.hg38)
ga <- GeneActivity(gex)
gex[['GA']] <- CreateAssayObject(counts = ga)

saveRDS(gex, snakemake@output[['obj']])