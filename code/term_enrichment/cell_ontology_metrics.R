library(Seurat)
library(Signac)
library(igraph)
library(GenomicRanges)
library(REMO.v1.GRCh38)
library(ontoProc)
library(ggplot2)
set.seed(1234)

# load data
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

cl <- get_ontology("data/term_enrichment/cl.obo", extract_tags = "everything")

# convert cl.obo to igraph
parents <- cl$parents
self <- rep(names(parents), lengths(parents))
g <- make_graph(rbind(unlist(parents), self))

# remove non cell ontology nodes
cl_nodes <- V(g)[grepl("^CL:", V(g)$name)]
g_filtered <- induced_subgraph(g, cl_nodes)

# make lookup vector for remo terms and CL_IDs
idx <- names(cl$name)[which(names(cl$name) %in% names(REMO.v1.GRCh38.CL_ID))]
remo_cl <- cl$name[idx]

# get a dataframe of wrong terms
wrong_terms <- function(obj, predictions, lookup, label_match) {
    # get top1 enriched terms
    clusters <- unique(obj$cluster_celltype)
    top_term_list <- c()
    for (id in clusters) {
        fgsea_results <- predictions[[id]]
        top_term_list[[id]] <- fgsea_results$pathway[1]
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

    # get wrong terms
    wrong_terms <- top_term_list[!exact_match]

    if (length(wrong_terms) == 0) {
        message("no wrong terms, returning NA")
        return(NA)
    }

    wrong_term_df <- data.frame(
    cluster = names(wrong_terms),
    enriched_terms = as.character(wrong_terms)
    )

    # add cluster CL ID, from label match table
    wrong_term_df$cluster_CL.ID <- sapply(wrong_term_df$cluster, function(x) {
        paste(label_match$CL_ID[label_match$label_transfer_name == x], collapse = ";")
    })

    # add enriched term CL ID, from REMO.v1.GRCh8.CL_ID
    wrong_term_df$enriched_terms_CL.ID <- sapply(wrong_term_df$enriched_terms, function(x) {
        paste(lookup[x])
    })

    return(wrong_term_df)
}

# calculate distance between enriched term and cluster term
node_dist <- function(wt_df_row, g) {
    wrong_term_cl.id <- wt_df_row$enriched_terms_CL.ID
    cluster_cl.id <- wt_df_row$cluster_CL.ID
    if (grepl(";", cluster_cl.id)) {
      cluster_cl.id <- unlist(strsplit(unlist(cluster_cl.id), ";"))
    }

    dist_n <- 100 # initialize high number
    for (j in 1:length(cluster_cl.id)) {
        dist <- shortest_paths(
            graph = g,
            from = wrong_term_cl.id,
            to = cluster_cl.id[j],
            mode = c("all"),
            output = c("both"),
            algorithm = c("automatic"))
        dist_n <- min(length(dist$vpath[[1]]), dist_n)
    }
    return(dist_n)
}

# calculate lcad
calculate_shortest_dist <- function(wt_df, g, cl) {
    cluster_dist <- c()
    
    unique_clusters <- unique(wt_df$cluster)
    for (i in seq_along(unique_clusters)) {
        cluster_name <- unique_clusters[i]
    
        # get term ids
        row <- wt_df[wt_df$cluster == cluster_name, c("cluster_CL.ID", "enriched_terms_CL.ID")]
        
        # calculate distance between enriched term and cluster term
        d <- node_dist(row, g)
        cluster_dist <- c(cluster_dist, d)
    }
    
    wt_df$shortest_dist <- cluster_dist
    return(wt_df)
}

pbmc_wt_wholebody <- wrong_terms(obj = pbmc_obj, 
                                 predictions = pbmc_predictions, 
                                 lookup = setNames(names(remo_cl), remo_cl), 
                                 label_match = label_match)
pbmc_wt_organ <- wrong_terms(obj = pbmc_obj, 
                             predictions = pbmc_predictions_organ, 
                             lookup = setNames(names(remo_cl), remo_cl), 
                             label_match = label_match)

islet_wt_wholebody <- wrong_terms(obj = islet_obj, 
                                  predictions = islet_predictions, 
                                  lookup = setNames(names(remo_cl), remo_cl), 
                                  label_match = label_match)
islet_wt_organ <- wrong_terms(obj = islet_obj, 
                              predictions = islet_predictions_organ, 
                              lookup = setNames(names(remo_cl), remo_cl), 
                              label_match = label_match)

brain_wt_wholebody <- wrong_terms(obj = brain_obj, 
                                  predictions = brain_predictions, 
                                  lookup = setNames(names(remo_cl), remo_cl), 
                                  label_match = label_match)
brain_wt_organ <- wrong_terms(obj = brain_obj, 
                              predictions = brain_predictions_organ, 
                              lookup = setNames(names(remo_cl), remo_cl), 
                              label_match = label_match)
# note: brain got no wrong enrichment

# randomize assignment
pbmc_wt_random_wholebody <- pbmc_wt_wholebody
pbmc_wt_random_wholebody$enriched_terms_CL.ID <- sample(names(remo_cl), nrow(pbmc_wt_random_wholebody), replace = TRUE)
pbmc_wt_random_organ <- pbmc_wt_organ
pbmc_wt_random_organ$enriched_terms_CL.ID <- sample(names(remo_cl), nrow(pbmc_wt_random_organ), replace = TRUE)

islet_wt_random_wholebody <- islet_wt_wholebody
islet_wt_random_wholebody$enriched_terms_CL.ID <- sample(names(remo_cl), nrow(islet_wt_random_wholebody), replace = TRUE)
islet_wt_random_organ <- islet_wt_organ
islet_wt_random_organ$enriched_terms_CL.ID <- sample(names(remo_cl), nrow(islet_wt_random_organ), replace = TRUE)

# calculate LCAD
pbmc_wt_wholebody <- calculate_shortest_dist(pbmc_wt_wholebody, g_filtered, cl)
pbmc_wt_organ <- calculate_shortest_dist(pbmc_wt_organ, g_filtered, cl)
islet_wt_wholebody <- calculate_shortest_dist(islet_wt_wholebody, g_filtered, cl)
islet_wt_organ <- calculate_shortest_dist(islet_wt_organ, g_filtered, cl)

pbmc_wt_random_wholebody <- calculate_shortest_dist(pbmc_wt_random_wholebody, g_filtered, cl)
pbmc_wt_random_organ <- calculate_shortest_dist(pbmc_wt_random_organ, g_filtered, cl)
islet_wt_random_wholebody <- calculate_shortest_dist(islet_wt_random_wholebody, g_filtered, cl)
islet_wt_random_organ <- calculate_shortest_dist(islet_wt_random_organ, g_filtered, cl)

# dataframes
pbmc_wt_wholebody$tissue <- "pbmc"
pbmc_wt_wholebody$enrichment <- "wholebody"
pbmc_wt_organ$tissue <- "pbmc"
pbmc_wt_organ$enrichment <- "organ"

islet_wt_wholebody$tissue <- "islet"
islet_wt_wholebody$enrichment <- "wholebody"
islet_wt_organ$tissue <- "islet"
islet_wt_organ$enrichment <- "organ"

pbmc_wt_random_wholebody$tissue <- "pbmc"
pbmc_wt_random_wholebody$enrichment <- "wholebody"
pbmc_wt_random_organ$tissue <- "pbmc"
pbmc_wt_random_organ$enrichment <- "organ"

islet_wt_random_wholebody$tissue <- "islet"
islet_wt_random_wholebody$enrichment <- "wholebody"
islet_wt_random_organ$tissue <- "islet"
islet_wt_random_organ$enrichment <- "organ"

all_lcad_organ <- rbind(pbmc_wt_organ, islet_wt_organ)
all_lcad_wholebody <- rbind(pbmc_wt_wholebody, islet_wt_wholebody)

all_lcad_random_organ <- rbind(pbmc_wt_random_organ, islet_wt_random_organ)
all_lcad_random_wholebody <- rbind(pbmc_wt_random_wholebody, islet_wt_random_wholebody)

true_enrichment_lcad <- rbind(all_lcad_organ, all_lcad_wholebody)
write.csv(true_enrichment_lcad, "data/term_enrichment/enrichment_shortestdist.csv", row.names = FALSE)

random_enrichment_lcad <- rbind(all_lcad_random_organ, all_lcad_random_wholebody)
for (i in 1:nrow(random_enrichment_lcad)) {
    random_cl <- random_enrichment_lcad[i,]$enriched_terms_CL.ID
    random_term_name <- remo_cl[random_cl]
    random_enrichment_lcad[i,]$enriched_terms <- random_term_name
}
write.csv(random_enrichment_lcad, "data/term_enrichment/random_enrichment_shortestdist.csv", row.names = FALSE)


# dist plot
mean_dist_organ <- mean(all_lcad_organ$shortest_dist)
mean_dist_wholebody <- mean(all_lcad_wholebody$shortest_dist)
mean_dist_random_organ <- mean(all_lcad_random_organ$shortest_dist)
mean_dist_random_wholebody <- mean(all_lcad_random_wholebody$shortest_dist)

df <- data.frame(
    group = c("organ", "wholebody", "organ", "wholebody"),
    mean_dist = c(mean_dist_organ, mean_dist_wholebody, mean_dist_random_organ, mean_dist_random_wholebody),
    enrichment = c("true enrichment", "true enrichment", "randomized", "randomized")
)
df$enrichment <- factor(df$enrichment, levels = c("true enrichment", "randomized"))

p <- ggplot(df, aes(x = enrichment, y = mean_dist, fill = group)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  labs(
    x = "CL term",
    y = "Mean shortest dist",
    fill = "term enrichment"
  ) +
  theme_minimal()

ggsave(file = "plots/term_enrichment_shortestdist.png", p, width = 5, height = 4, dpi = 500)