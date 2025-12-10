library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(patchwork)
library(stringr)
library(grid)
library(gridExtra)
library(ggrepel)

source("code/utilities.R")

pbmc_obj <- readRDS("data/term_enrichment/pbmc/pbmc_remo.rds")
pbmc_predictions <- readRDS("data/term_enrichment/pbmc/pbmc_predictions.rds")
pbmc_predictions_organ <- readRDS("data/term_enrichment/pbmc/pbmc_predictions_organ.rds")

islet_obj <- readRDS("data/term_enrichment/islet/islet_remo.rds")
islet_predictions <- readRDS("data/term_enrichment/islet/islet_predictions.rds")
islet_predictions_organ <- readRDS("data/term_enrichment/islet/islet_predictions_organ.rds")

brain_obj <- readRDS("data/term_enrichment/brain/brain_remo.rds")
brain_predictions <- readRDS("data/term_enrichment/brain/brain_predictions.rds")
brain_predictions_organ <- readRDS("data/term_enrichment/brain/brain_predictions_organ.rds")

label_match <- read.table("data/term_enrichment/label_match.csv", sep = ",", header = TRUE)

## calculate accuracy score
accuracy_score <- function(obj, predictions) {
    # get top1 and top3 enriched terms
    clusters <- unique(obj$cluster_celltype)
    top_term_list <- c()
    top3_terms_list <- list()
    for (id in clusters) {
        fgsea_results <- predictions[[id]]
        top_term_list[[id]] <- fgsea_results$pathway[1]
        top3_terms_list[[id]] <- fgsea_results$pathway[1:3]
    }

    # if no enrichment, set as undefined
    top_term_list <- lapply(top_term_list, function(x) {
        if (is.na(x)) "undefined" else x
    })

    # for each cell type, check if to enriched term is exact match
    labels <- names(top_term_list)
    exact_match <- sapply(labels, function(lb) {
      predicted <- top_term_list[[lb]]
      any(label_match$gene_ontology_name[label_match$label_transfer_name == lb] == predicted)
    })

    # for each cell type, check if top 3 enriched term is exact match 
    labels <- names(top3_terms_list)
    top3_match <- sapply(labels, function(lb) {
      predicted <- top3_terms_list[[lb]]  
      truth <- label_match$gene_ontology_name[label_match$label_transfer_name == lb]
      any(predicted %in% truth)
    })

    accuracy_1 <- mean(exact_match)
    accuracy_3 <- mean(top3_match)

    accuracy_metric <- c(accuracy_1, accuracy_3)
    names(accuracy_metric) <- c("accuracy_top1", "accuracy_top3")

    return(accuracy_metric)
}

pbmc_score <- accuracy_score(pbmc_obj, pbmc_predictions)
brain_score <- accuracy_score(brain_obj, brain_predictions)
islet_score <- accuracy_score(islet_obj, islet_predictions)

pbmc_score_organ <- accuracy_score(pbmc_obj, pbmc_predictions_organ)
brain_score_organ <- accuracy_score(brain_obj, brain_predictions_organ)
islet_score_organ <- accuracy_score(islet_obj, islet_predictions_organ)

df <- data.frame(
  tissue = rep(c("pbmc", "brain", "islet"), each = 2),
  metric = rep(c("Top term", "Within top 3"), times = 3),
  accuracy = c(
    pbmc_score[['accuracy_top1']], pbmc_score[['accuracy_top3']],
    brain_score[['accuracy_top1']], brain_score[['accuracy_top3']],
    islet_score[['accuracy_top1']], islet_score[['accuracy_top3']]
  )
)
df$mapping <- "Whole body"

df2 <- data.frame(
  tissue = rep(c("pbmc", "brain", "islet"), each = 2),
  metric = rep(c("Top term", "Within top 3"), times = 3),
  accuracy = c(
    pbmc_score_organ[['accuracy_top1']], pbmc_score_organ[['accuracy_top3']],
    brain_score_organ[['accuracy_top1']], brain_score_organ[['accuracy_top3']],
    islet_score_organ[['accuracy_top1']], islet_score_organ[['accuracy_top3']]
  )
)
df2$mapping <- "Tissue"

df <- rbind(df, df2)

p <- ggplot(df, aes(x = tissue, y = accuracy, fill = mapping)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~metric, ncol = 1) +
  scale_fill_manual(values = c("#FFA500", "#FFDFA3")) +
  ylim(0, 1) +
  labs(
    title = "REMO cell type\nannotation accuracy",
    x = "Dataset",
    y = "Accuracy",
    fill = "Mapping"
  ) +
  theme_minimal() + theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))


ggsave(file = "plots/term_enrichment_metrics.pdf", p, width = 3, height = 10)


plot_cell_annotation <- function(object, label_match) {
    
    df <- as.data.frame(Embeddings(object[['umap']]))
    
    df$celltype <- object$cluster_celltype
    df$term_wholebody <- object$top_enriched_terms
    df$term_tissue <- object$top_enriched_terms_organ

    label_df <- df[, c("celltype", "term_wholebody", "term_tissue")] %>% unique()

    label_df$wb_match <- NA
    label_df$tissue_match <- NA
    for (i in seq_len(nrow(label_df))) {
        row <- label_df[i, ]
        match_labels <- label_match[label_match$label_transfer_name == row$celltype, ]
        wb_match <- row$term_wholebody %in% match_labels$gene_ontology_name
        tissue_match <- row$term_tissue %in% match_labels$gene_ontology_name
        label_df[i, 'wb_match'] <- wb_match
        label_df[i, 'tissue_match'] <- tissue_match
    }

    label_data_gt <- df %>%
      group_by(celltype) %>%
      summarize(
        umap_1 = mean(umap_1),
        umap_2 = mean(umap_2),
        .groups = "drop"
      ) %>%
      mutate(wrapped_label = str_wrap(celltype, width = 20))

    label_data_wb <- df %>%
      group_by(term_wholebody) %>%
      summarize(
        umap_1 = mean(umap_1),
        umap_2 = mean(umap_2),
        .groups = "drop"
      ) %>%
      left_join(label_df, by = "term_wholebody") %>%
      mutate(wrapped_label = str_wrap(term_wholebody, width = 20))

    label_data_tissue <- df %>%
      group_by(term_tissue) %>%
      summarize(
        umap_1 = mean(umap_1),
        umap_2 = mean(umap_2),
        .groups = "drop"
      ) %>%
      left_join(label_df, by = "term_tissue") %>%
      mutate(wrapped_label = str_wrap(term_tissue, width = 20))
    
    transparent_cols <- c("TRUE" = "#3498DB4D", "FALSE" = "#E74C3C4D")

    ct <- ggplot(df, aes(x = umap_1, y = umap_2)) +
      geom_point(aes(color = celltype), size = 0.5) +
      geom_label_repel(
        data = label_data_gt,
        aes(label = wrapped_label),
        color = "black",
        size = 3.5,
        box.padding = 0.5,
        point.padding = 0.5,
        min.segment.length = 0,
        max.overlaps = Inf
      ) +
      scale_color_discrete(guide = "none") +
      theme_bw() +
      theme(legend.position = "none")

    wb <- ggplot(df, aes(x = umap_1, y = umap_2)) +
      geom_point(aes(color = term_wholebody), size = 0.5) +
      geom_label_repel(
        data = label_data_wb,
        aes(
          label = wrapped_label,
          fill = wb_match
        ),
        color = "black",
        size = 3.5,
        box.padding = 0.5,
        point.padding = 0.5,
        min.segment.length = 0,
        max.overlaps = Inf
      ) +
      scale_fill_manual(values = transparent_cols) +
      scale_color_discrete(guide = "none") +
      theme_bw() +
      theme(legend.position = "none")

    tissue <- ggplot(df, aes(x = umap_1, y = umap_2)) +
      geom_point(aes(color = term_tissue), size = 0.5) +
      geom_label_repel(
        data = label_data_tissue,
        aes(
          label = wrapped_label,
          fill = tissue_match
        ),
        color = "black",
        size = 3.5,
        box.padding = 0.5,
        point.padding = 0.5,
        min.segment.length = 0,
        max.overlaps = Inf
      ) +
      scale_fill_manual(values = transparent_cols) +
      scale_color_discrete(guide = "none") +
      theme_bw() +
      theme(legend.position = "none")
    return(list(ct, wb, tissue))
}

## plot umaps
brain_plot <- plot_cell_annotation(brain_obj, label_match)
pbmc_plot <- plot_cell_annotation(pbmc_obj, label_match)
islet_plot <- plot_cell_annotation(islet_obj, label_match)

pp <- (brain_plot[[1]] | brain_plot[[2]] | brain_plot[[3]]) / (pbmc_plot[[1]] | pbmc_plot[[2]] | pbmc_plot[[3]]) / (islet_plot[[1]] | islet_plot[[2]] | islet_plot[[3]])

ggsave(file = "plots/term_enrichment_umap.png", plot = pp, height = 13, width = 13, dpi = 500)

## plot top 5 enriched terms
fgsea_table <- function(predictions, ct_name = "pDC") {
    df <- head(predictions[[ct_name]], 5)
    
    # subset first 3 leadingEdge
    df$leadingEdge <- sapply(df$leadingEdge, function(y) {
        y <- y[1:min(3, length(y))]   # keep first 3
        paste(y, collapse = ", ")     # collapse
    })

    colnames(df) <- c("pathway", "pval", "padj", "log2err", "ES", "NES", "size", "first 3 leadingEdge")
    df <- as(df, 'data.frame')

    # adjust max width
    df$pathway <- str_wrap(df$pathway, width = 25)
    df$'first 3 leadingEdge' <- str_wrap(df$'first 3 leadingEdge', width = 10)
    return(df)
}

pDC <- fgsea_table(predictions = pbmc_predictions, ct_name = 'pDC')
beta <- fgsea_table(predictions = islet_predictions, ct_name = 'beta')

pdf("plots/fgsea_table.pdf", width = 15, height = 5)
grid.newpage()
tg <- tableGrob(pDC)
grid.draw(tg)

grid.newpage()
tg <- tableGrob(beta)
grid.draw(tg)
dev.off()

## save table top 5 term enrichment per celltype
get_top_n <- function(x, tissue, n=5) {
    celltypes <- names(x)
    fg <- data.frame()
    for (i in seq_along(celltypes)) {
        df <- head(x[[celltypes[i]]], 5)
        df$'label transfer cell type' <- celltypes[i]
        df$rank <- 1:min(5,nrow(df))
        fg <- rbind(fg, df)
    }
    fg$tissue <- tissue
    fg <- fg[,c(11, 9, 10, 1:8)]
    return(fg)
}

brain_fgsea <- get_top_n(brain_predictions, 'brain')
pbmc_fgsea <- get_top_n(pbmc_predictions, 'pbmc')
islet_fgsea <- get_top_n(islet_predictions, 'islet')
all_fgsea <- rbind(brain_fgsea, pbmc_fgsea, islet_fgsea)

# show only first 5 leading edge
all_fgsea$leadingEdge <- sapply(
    all_fgsea$leadingEdge,
    function(x) {
        n <- length(x)
        if (n > 5) {
            x <- c(x[1:5], "...")     # add ellipsis
        }
        paste(x, collapse = ", ")
    }
)

write.table(x = all_fgsea, file = "data/term_enrichment/fgsea_wholebody.csv", sep=",", row.names=FALSE)


brain_fgsea <- get_top_n(brain_predictions_organ, 'brain')
pbmc_fgsea <- get_top_n(pbmc_predictions_organ, 'pbmc')
islet_fgsea <- get_top_n(islet_predictions_organ, 'islet')
all_fgsea <- rbind(brain_fgsea, pbmc_fgsea, islet_fgsea)

# show only first 5 leading edge
all_fgsea$leadingEdge <- sapply(
    all_fgsea$leadingEdge,
    function(x) {
        n <- length(x)
        if (n > 5) {
            x <- c(x[1:5], "...")     # add ellipsis
        }
        paste(x, collapse = ", ")
    }
)

write.table(x = all_fgsea, file = "data/term_enrichment/fgsea_tissue.csv", sep=",", row.names=FALSE)
