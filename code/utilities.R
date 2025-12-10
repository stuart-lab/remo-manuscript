PEAKS_COL="#E85C0D"
REMO_COL="#56B4E9"
GENE_COL="#1DCD9F"
CCRE_COL="#CE6DFF"
RUST_COL="#DEA584"
R_COL="#198CE7"

library(Matrix)
library(Signac)
library(Seurat)

set.seed(1234)


theme_dimplot <- function(..., legend = TRUE) {
  custom.theme <- theme_void() +
    theme(
        panel.border = element_rect(color = "grey", fill = NA, linewidth = 1)
    )
  if (!legend) {
    custom.theme <- custom.theme +
      theme(
        legend.position = "none"
      )
  }
  return(custom.theme)
}

get_nonzero <- function(x) {
    # read first line of mtx file to get nonzero elements
    matrix_dims <- read.table(x, sep = " ", comment.char = "%", nrows = 1)
    return(matrix_dims[1, 3])
}

get_knn_purity <- function(emb, idents, method) {
    knp <- data.frame(
        "Score" = knn_purity(emb, idents),
        "Celltype" = idents
    )
    knp <- knp |>
        dplyr::group_by(Celltype) |>
        dplyr::mutate(mn = mean(Score)) |>
        dplyr::ungroup() |>
        dplyr::select(Celltype, mn) |>
        unique()
    knp$Method <- method
    return(knp)
}

get_silhouette <- function(x, idents) {
    dist.matrix <- dist(x)
    sil <- cluster::silhouette(x = as.numeric(x = as.factor(x = idents)), dist = dist.matrix)
    return(sil[, 3])
}

get_mean_sil <- function(emb, idents, method) {
    sil <- data.frame(
        "Score" = get_silhouette(emb, idents),
        "Celltype" = idents
    )
    sil <- sil |>
        dplyr::group_by(Celltype) |>
        dplyr::mutate(mn = mean(Score)) |>
        dplyr::ungroup() |>
        dplyr::select(Celltype, mn) |>
        unique()
    sil$Method <- method
    return(sil)
}

get_metrics <- function(remo_obj, peaks_obj,
                        remo_emb, peaks_emb,
                        tissue_name) {

    # KNN purity
    knp_remo <- get_knn_purity(remo_emb, remo_obj$cluster, "REMO")
    knp_peaks <- get_knn_purity(peaks_emb, peaks_obj$cluster, "Peaks")
    knp <- rbind(knp_remo, knp_peaks)
    knp$Dataset <- tissue_name
    knp$Method <- factor(knp$Method, levels = c("REMO", "Peaks"))

    # silhouette
    sil_remo <- get_mean_sil(remo_emb, remo_obj$cluster, "REMO")
    sil_peaks <- get_mean_sil(peaks_emb, peaks_obj$cluster, "Peaks")
    sil <- rbind(sil_remo, sil_peaks)
    sil$Dataset <- tissue_name
    sil$Method <- factor(sil$Method, levels = c("REMO", "Peaks"))

    # CH
    ch <- data.frame(
        "Dataset" = tissue_name,
        "Score" = c(
            fpc::calinhara(remo_emb, as.numeric(remo_obj$cluster)),
            fpc::calinhara(peaks_emb, as.numeric(peaks_obj$cluster))),
        "Method" = c("REMO", "Peaks"),
        "Metric" = "CH"
    )
    
    # ARI
    ari <- data.frame(
        "Dataset" = tissue_name,
        "Score" = c(
            mclust::adjustedRandIndex(remo_obj$seurat_clusters, remo_obj$cluster),
            mclust::adjustedRandIndex(peaks_obj$seurat_clusters, peaks_obj$cluster)),
        "Method" = c("REMO", "Peaks"),
        "Metric" = "ARI"
    )

    # NMI
    nmi <- data.frame(
        "Dataset" = tissue_name,
        "Score" = c(
            aricode:::NMI(remo_obj$seurat_clusters, remo_obj$cluster),
            aricode:::NMI(peaks_obj$seurat_clusters, peaks_obj$cluster)),
        "Method" = c("REMO", "Peaks"),
        "Metric" = "NMI"
    )

    return(list("AN" = rbind(ari, nmi, ch), "KNN" = knp, "SIL" = sil))
}

knn_purity <- function(embeddings, clusters, k = 20) {
  nn <- RANN::nn2(data = embeddings, k = k + 1)$nn.idx[, 2:k]
  nn_purity <- vector(mode = "numeric", length = length(x = clusters))
  for (i in seq_len(length.out = nrow(x = nn))) {
    nn_purity[i] <- sum(clusters[nn[i, ]] == clusters[i]) / k
  }
  return(nn_purity)
}

# Function to extract runtime and max memory from the profiling output
extract_profile_info <- function(profile_file) {
  # Read the profiling file
  lines <- readLines(profile_file)
  
  # Initialize variables to store the results
  elapsed_time <- NA
  max_memory <- NA
  
  # Loop through each line and extract relevant information
  for (line in lines) {
    if (grepl("Command exited with non-zero status 1", line)) {
      df <- data.frame("Time" = NA, "Memory" = NA)
      return(df)
    }
    if (grepl("Elapsed \\(wall clock\\) time", line)) {
      # Extract the elapsed time
      elapsed_time <- convert_time_to_seconds(sub(".*: ", "", line))
    } else if (grepl("Maximum resident set size", line)) {
      # Extract the max memory usage
      max_memory <- as.numeric(sub(".*: ", "", line))
    }
  }

  df <- data.frame("Time" = elapsed_time, "Memory" = max_memory/1024)
  return(df)
}

convert_time_to_seconds <- function(time_str) {
  # Check if the time format is h:mm:ss or m:ss
  if (grepl(":", time_str)) {
    time_parts <- as.numeric(unlist(strsplit(time_str, ":")))
    if (length(time_parts) == 3) {
      # Format is h:mm:ss
      return(time_parts[1] * 3600 + time_parts[2] * 60 + time_parts[3])
    } else if (length(time_parts) == 2) {
      # Format is m:ss
      return(time_parts[1] * 60 + time_parts[2])
    }
  }
  return(as.numeric(time_str))
}

revcomp <- function(x) {
    as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(x)))
}

process_rna_obj <- function(
    counts,
    dims = 1:30,
    var_features = NULL,
    nfeatures = 2000,
    scale.factor=10000,
    resolution = 0.8
) {
  obj <- CreateSeuratObject(counts = counts)
  obj <- NormalizeData(obj, scale.factor = scale.factor)
  if (is.null(var_features)) {
    obj <- FindVariableFeatures(obj, nfeatures = nfeatures)          
  } else {
    VariableFeatures(obj) <- var_features
  }
  obj <- ScaleData(obj)
  obj <- RunPCA(obj, verbose = FALSE)
  obj <- RunUMAP(obj, reduction = 'pca', dims = dims, verbose = FALSE)
  obj <- FindNeighbors(obj, reduction = "pca", dims = dims)
  obj <- FindClusters(obj, algorithm=3, resolution=resolution)
  return(obj)
}

process_remo_obj <- function(
    counts,
    dims = 1:30,
    var_features = NULL,
    nfeatures = 20000,
    resolution = 0.8
) {
  obj <- CreateSeuratObject(counts = counts)
  obj <- NormalizeData(obj, scale.factor = median(colSums(counts)))
  if (is.null(var_features)) {
    obj <- FitMeanVar(obj, nfeatures = nfeatures)
  } else {
    VariableFeatures(obj) <- var_features
  }
  obj <- RunSVD(obj, pca = TRUE)
  obj <- RunUMAP(obj, reduction = 'pca', dims = dims, verbose = FALSE)
  obj <- FindNeighbors(obj, reduction = "pca", dims = dims)
  obj <- FindClusters(obj, algorithm=3, resolution=resolution)
  return(obj)
}

process_atac_obj <- function(counts, dims = 2:30, features = NULL, resolution=0.8) {
  obj <- CreateSeuratObject(counts = counts, assay = "ATAC")
  obj <- RunTFIDF(obj)
  if (is.null(x = features)) {
      obj <- FindTopFeatures(obj)
  } else {
      VariableFeatures(obj) <- features
  }
  obj <- RunSVD(obj)
  obj <- RunUMAP(obj, reduction = 'lsi', dims = dims, verbose = FALSE)
  obj <- FindNeighbors(obj, reduction = "lsi", dims = dims)
  obj <- FindClusters(obj, algorithm=3, resolution=resolution)
  return(obj)
}
