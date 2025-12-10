library(Seurat)

ref <- snakemake@input[['ref']]
query <- snakemake@input[['obj']]
output <- snakemake@output[[1]]

ref <- readRDS(ref)
ref <- UpdateSeuratObject(ref)

query <- readRDS(query)

transfer.anchors <- FindTransferAnchors(
  reference = ref,
  query = query,
  reduction = 'pcaproject'
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = ref$celltype,
  weight.reduction = query[['pca']],
  dims = 1:30
)

saveRDS(object = predicted.labels, file = output)