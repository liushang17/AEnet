library(Seurat)
library(ggplot2)
library(umap)

#' @title ASP Cluster Enrichment Score
#' @description Calculates enrichment scores of ASP clusters for each cell.
#'
#' @param psi_matrix A numeric matrix: rows = junctions, columns = cells (e.g., PSI or ratio).
#' @param ann_junc A data frame of ASP annotations, must contain columns 'symbol' and 'PC' (cluster ID).
#'
#' @return A data frame of enrichment scores: rows = cells, columns = cluster enrichment scores.
#' @export
#'
#' @examples
#' enrichment_scores <- asp_score(psi_matrix, ann_junc)
library(dplyr)
asp_score <- function(psi_matrix, ann_junc) {
  stopifnot("symbol" %in% colnames(ann_junc), "PC" %in% colnames(ann_junc))
  stopifnot(all(rownames(psi_matrix) != ""))
  
  # Remove 'Other' patterns
  ann_junc <- subset(ann_junc, PC != "Other")
  
  # Get unique clusters
  clusters <- unique(as.character(ann_junc$PC))
  all_scores <- list()
  
  for (clus in clusters) {
    # Subset junctions for this cluster
    juncs_in_cluster <- ann_junc$symbol[ann_junc$PC == clus]
    matched_juncs <- intersect(juncs_in_cluster, rownames(psi_matrix))
    
    if (length(matched_juncs) < 2) next  # Skip small clusters
    
    psi_subset <- psi_matrix[matched_juncs, , drop = FALSE]
    
    # Remove rows with zero variance or all NA
    row_sds <- apply(psi_subset, 1, sd, na.rm = TRUE)
    psi_subset <- psi_subset[row_sds > 0, , drop = FALSE]
    
    if (nrow(psi_subset) == 0) next
    
    # Compute per-cell mean for the cluster
    cluster_scores <- colMeans(psi_subset, na.rm = TRUE)
    
    # Store as data frame
    cluster_df <- data.frame(Cell = names(cluster_scores), Score = cluster_scores)
    colnames(cluster_df)[2] <- clus
    
    all_scores[[clus]] <- cluster_df
  }
  
  # Merge all cluster scores by Cell
  if (length(all_scores) == 0) {
    warning("No valid clusters found.")
    return(data.frame())
  }
  
  score_df <- Reduce(function(x, y) merge(x, y, by = "Cell", all = TRUE), all_scores)
  return(score_df)
}

#' @title Gene Cluster Enrichment Score
#' @description Calculates enrichment scores of gene clusters for each cell based on expression matrix.
#'
#' @param exp A numeric matrix: rows = genes, columns = cells (expression).
#' @param geneann A data frame with columns 'symbol' (gene) and 'PC' (cluster ID).
#'
#' @return A data frame of enrichment scores: rows = cells, columns = clusters.
#' @export
#'
#' @examples
#' enrichment_scores <- gene_score(exp, geneann)

gene_score <- function(exp, geneann) {
  stopifnot(is.matrix(exp))
  stopifnot(all(c("symbol", "PC") %in% colnames(geneann)))
  
  # Filter out 'Other' cluster
  geneann <- subset(geneann, PC != "Other")
  
  # Unique cluster IDs
  clusters <- unique(as.character(geneann$PC))
  score_list <- list()
  
  for (clus in clusters) {
    genes_in_cluster <- geneann$symbol[geneann$PC == clus]
    matched_genes <- intersect(genes_in_cluster, rownames(exp))
    
    if (length(matched_genes) < 2) next  # skip small clusters
    
    sub_exp <- exp[matched_genes, , drop = FALSE]
    cluster_score <- colMeans(sub_exp, na.rm = TRUE)
    
    score_list[[clus]] <- cluster_score
  }
  
  if (length(score_list) == 0) {
    warning("No valid clusters with matching genes found.")
    return(data.frame())
  }
  
  score_matrix <- do.call(cbind, score_list)
  score_df <- log2(as.data.frame(score_matrix)+1)
  score_df$Cell <- colnames(exp)
  score_df <- score_df[, c("Cell", setdiff(colnames(score_df), "Cell"))]
  rownames(score_df) <- NULL
  
  return(score_df)
}

#' @title Cell Clusters
#' @description Cluster cells based on ASP or gene expression scores using Seurat.
#'
#' @param asp_score A data frame with ASP enrichment scores
#' @param exp_score A data frame with gene expression scores
#' @param resolution Clustering resolution (default=0.5)
#' @param min.dist Minimum distance for UMAP (default=1)
#' @param scale Whether to scale data (default=TRUE)
#' @param log_transform Whether to log-transform exp_score (default=TRUE)
#'
#' @return A data frame with Cell, cluster, UMAP coordinates
#' @export
cell_clus <- function(asp_score = NULL, exp_score = NULL, resolution = 0.5, 
                      min.dist = 1, scale = TRUE) {
  
  # Validate input - 修复的if语句
  if (is.null(asp_score) && is.null(exp_score)) {
    stop("At least one of asp_score or exp_score must be provided")
  }
  
  # Process ASP scores
  if (!is.null(asp_score)) {
    if (!"Cell" %in% colnames(asp_score)) {
      stop("asp_score must contain a 'Cell' column")
    }
    rownames(asp_score) <- asp_score$Cell
    asp_mat <- asp_score[, setdiff(colnames(asp_score), "Cell"), drop = FALSE]
    if (scale) {asp_mat <- t(scale(t(scale(asp_mat))))}else{ asp_mat <- scale(asp_mat)}
  }
  
  # Process expression scores
  if (!is.null(exp_score)) {
    if (!"Cell" %in% colnames(exp_score)) {
      stop("exp_score must contain a 'Cell' column")
    }
    rownames(exp_score) <- exp_score$Cell
    exp_mat <- exp_score[, setdiff(colnames(exp_score), "Cell"), drop = FALSE]
    if (scale){ exp_mat <- t(scale(t(scale(exp_mat))))}else{exp_mat <- scale(exp_mat)}
  }
  
  # Merge data
  if (!is.null(asp_score) && !is.null(exp_score)) {
    common_cells <- intersect(rownames(asp_mat), rownames(exp_mat))
    if (length(common_cells) == 0) stop("No common cells between inputs")
    mat <- cbind(asp_mat[common_cells, ], exp_mat[common_cells, ])
  } else {
    mat <- if (!is.null(asp_score)) asp_mat else exp_mat
  }

  # NA
  na.omit(mat)
  
  # Clustering
  snn <- Seurat::FindNeighbors(mat)$snn
  clusters <- Seurat::FindClusters(snn, resolution = resolution)
  umap <- Seurat::RunUMAP(mat, min.dist = min.dist)
 
  # Return results
  data.frame(
    Cell = rownames(mat),
    cluster = clusters[, 1],
    umap_1 = umap@cell.embeddings[, 1],
    umap_2 = umap@cell.embeddings[, 2],
    row.names = NULL
  )
}

#' @title Identify Key Splicing Factors
#' @description Identifies key splicing factors (SFs) for specific clusters from AEN network data
#'
#' @param corm A data frame containing the AEN network correlations (must contain 'symbol', 'junction', and 'correlation' columns)
#' @param sf A character vector of splicing factor gene symbols to consider
#' @param ann_junc A data frame containing annotated junctions with cluster assignments (must contain 'PC' and 'symbol' columns)
#' @param cluster The specific cluster identifier to analyze (default = "C_1")
#' @param min_connections Minimum number of connections required to consider an SF (default = 10)
#' @param cor_threshold Threshold for considering dominant correlation direction (default = 0.75)
#'
#' @return A data frame with columns: 
#'   - SF: splicing factor symbol
#'   - importance: proportion of cluster junctions connected to the SF
#'   - order: ranking by number of connections
#'   - direction: dominant correlation direction (1 = positive, -1 = negative, 0 = mixed)
#' @export
#'
#' @examples 
#' key_sf <- key_sf(corm, sf_list, ann_junc, cluster = "C_1")
key_sf <- function(corm, sf, ann_junc, cluster = "C_1", min_connections = 10, cor_threshold = 0.75) {
  
  # Input validation
  if (!all(c("symbol", "junction", "correlation") %in% colnames(corm))) {
    stop("corm must contain 'symbol', 'junction', and 'correlation' columns")
  }
  if (!all(c("PC", "symbol") %in% colnames(ann_junc))) {
    stop("ann_junc must contain 'PC' and 'symbol' columns")
  }
  if (!is.character(sf)) {
    stop("sf must be a character vector")
  }
  
  # Filter junctions for the target cluster
  cluster_junc <- ann_junc[ann_junc$PC == cluster, ]
  if (nrow(cluster_junc) == 0) {
    stop("No junctions found for cluster ", cluster)
  }
  
  # Find SF-junction connections
  sf_connections <- corm[corm$symbol %in% sf & corm$junction %in% cluster_junc$symbol, ]
  
  # Check if we have enough connections
  if (nrow(sf_connections) < min_connections) {
    message("Insufficient connections (", nrow(sf_connections), ") for cluster ", cluster)
    return(NULL)
  }
  
  # Calculate SF importance metrics
  sf_stats <- sf_connections %>%
    dplyr::group_by(symbol) %>%
    dplyr::summarize(
      Freq = dplyr::n(),
      pro = Freq / nrow(cluster_junc),
      pos_cor = sum(correlation > 0) / Freq,
      neg_cor = sum(correlation < 0) / Freq
    ) %>%
    dplyr::arrange(dplyr::desc(Freq)) %>%
    dplyr::mutate(
      order = seq_len(dplyr::n()),
      direction = dplyr::case_when(
        pos_cor > cor_threshold ~ 1L,
        neg_cor > cor_threshold ~ -1L,
        TRUE ~ 0L
      )
    ) %>%
    dplyr::select(
      SF = symbol,
      importance = pro,
      order,
      direction
    )
  
  return(sf_stats)
}

