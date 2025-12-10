library(Signac)
library(Seurat)
library(Matrix)
library(ggplot2)
library(REMO.v1.GRCh38)
source("code/utilities.R")

dims <- 20
resolution <- 0.8

obj <- readRDS(snakemake@input[['multiome']])
remo_counts_dir <- snakemake@input[['remo']]

grouped <- Read10X(remo_counts_dir, gene.column=1)
grouped_filt <- grouped[, colnames(obj)]

remo_obj <- process_remo_obj(grouped_filt, dims = 1:dims, resolution = resolution)

# add multiome annotations
remo_obj <- AddMetaData(remo_obj, obj$predicted.id, col.name = "multiome_annotation")

remo_obj <- FindNeighbors(remo_obj, reduction = "pca", dims = 1:dims)
remo_obj <- FindClusters(remo_obj, algorithm=3, resolution=1, graph.name = "RNA_snn")

# map multiome annotations to remo clusters
cluster_id <- remo_obj$seurat_clusters
names(cluster_id) <- remo_obj$multiome_annotation

remo_obj$prediction_score <- obj$prediction.score.max
prediction_score <- remo_obj$prediction_score
names(prediction_score) <- remo_obj$multiome_annotation

ct_list <- c()
for (i in unique(cluster_id)) {
    idx <- which(cluster_id == i)
    cluster_ct <- cluster_id[idx]
    cluster_score <- prediction_score[idx]
    
    highscore <- which(cluster_score > 0.5)
    cluster_score <- cluster_score[highscore]
    
    ct_table <- table(names(cluster_score))
    majority_ct <- names(sort(ct_table, decreasing=TRUE))[1]
    ct_list[i] <- majority_ct
}

# add majority multiome celltype to remo seurat clusters
current.cluster.ids <- names(ct_list)
new.cluster.ids <- ct_list
remo_obj$cluster_celltype <- plyr::mapvalues(x = remo_obj$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)

# term enrichment whole body
Idents(remo_obj) <- remo_obj$cluster_celltype
predictions <- EnrichedTerms(remo_obj, terms = REMO.v1.GRCh38.terms)


# organ specific
cl <- REMO.v1.GRCh38.tissues$Blood
predictions_organ <- EnrichedTerms(remo_obj, terms = REMO.v1.GRCh38.terms[cl])

# add top enriched term as cluster labels
clusters <- unique(remo_obj$cluster_celltype)
top_term_list <- c()
top_term_list_organ <- c()
for (id in clusters) {
    fgsea_results <- predictions[[id]]
    fgsea_results_organ <- predictions_organ[[id]]
    top_term_list[[id]] <- fgsea_results$pathway[1]
    top_term_list_organ[[id]] <- fgsea_results_organ$pathway[1]
}

# if no enrichment, set as undefined
top_term_list <- lapply(top_term_list, function(x) {
    if (is.na(x)) "undefined" else x
})
top_term_list_organ <- lapply(top_term_list_organ, function(x) {
    if (is.na(x)) "undefined" else x
})

current.cluster.ids <- names(top_term_list)
new.cluster.ids <- sapply(top_term_list, function(x) x[1])
remo_obj$top_enriched_terms <- plyr::mapvalues(x = remo_obj$cluster_celltype, from = current.cluster.ids, to = new.cluster.ids)

new.cluster.ids <- sapply(top_term_list_organ, function(x) x[1])
remo_obj$top_enriched_terms_organ <- plyr::mapvalues(x = remo_obj$cluster_celltype, from = current.cluster.ids, to = new.cluster.ids)

saveRDS(remo_obj, snakemake@output[['remo']])
saveRDS(predictions, snakemake@output[['predictions']])
saveRDS(predictions_organ, snakemake@output[['organ']])