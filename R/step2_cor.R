#' @title Spearman Correlation Between Junction Ratio and Gene Expression
#' @description Computes Spearman correlation between exon junction ratios and gene expression across cells.
#'              Designed for single-sample processing with internal parallelization.
#' @param ratio_mat A numeric matrix of exon junction ratios (rows: junction pairs, columns: cells).
#' @param gene_mat A numeric matrix of gene expression (rows: genes, columns: cells).
#' @param min_obs Minimum number of valid observations (default = 21).
#' @param parallel Logical. Whether to use parallel computation (default = FALSE).
#' @param n_cores Number of cores to use in parallel mode. Ignored if parallel = FALSE.
#' @return A data frame of significant Spearman correlations (rows: junction-gene pairs).
#' @export
#' @importFrom Rfast correls
#' @import parallel
junction_gene_spearman <- function(ratio_mat, gene_mat, min_obs = 21, min_fea = 6, parallel = FALSE, n_cores = 2) {
  
  # Align cells
  common_cells <- intersect(colnames(ratio_mat), colnames(gene_mat))
  ratio_mat <- ratio_mat[, common_cells]
  gene_mat <- gene_mat[, common_cells]
  
  # Input checks
  stopifnot(
    identical(colnames(ratio_mat), colnames(gene_mat)),
    is.matrix(ratio_mat),
    is.matrix(gene_mat)
  )
  
  # Ensure row names
  if (is.null(rownames(ratio_mat))) rownames(ratio_mat) <- paste0("JP", seq_len(nrow(ratio_mat)))
  if (is.null(rownames(gene_mat))) rownames(gene_mat) <- paste0("G", seq_len(nrow(gene_mat)))
  
  compute_cor <- function(i) {
    x <- ratio_mat[i, ]
    valid_cells <- !is.na(x)
    if (sum(valid_cells) < min_obs) return(NULL)
    
    gene_mat_valid <- gene_mat[, valid_cells]
    gene_mat_binary <- gene_mat_valid
    gene_mat_binary[gene_mat_binary > 0] <- 1
    gene_freq <- rowSums(gene_mat_binary)
    valid_genes <- names(gene_freq[gene_freq >= min_fea])
    if (length(valid_genes) == 0) return(NULL)
    
    corr <- Rfast::correls(
      y = as.numeric(x[valid_cells]),
      x = t(gene_mat[valid_genes, valid_cells]),
      type = "spearman"
    )
    
    sig_idx <- which(corr[, 5] < 0.01)
    if (length(sig_idx) >= 2) {
      result <- data.frame(corr[sig_idx, ])
      result$geneID <- rownames(corr)[sig_idx]
      result$junctionID <- rownames(ratio_mat)[i]
      return(result)
    }
    return(NULL)
  }
  
  if (parallel) {
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl))
    parallel::clusterExport(cl, varlist = c("ratio_mat", "gene_mat", "min_obs", "compute_cor"), envir = environment())
    parallel::clusterEvalQ(cl, library(Rfast))
    cor_list <- parallel::parLapply(cl, seq_len(nrow(ratio_mat)), compute_cor)
  } else {
    cor_list <- lapply(seq_len(nrow(ratio_mat)), compute_cor)
  }
  
  cor_df <- do.call(rbind, cor_list)
  return(cor_df)
}

#' @title Multi-Sample ASP-Gene Correlation Analysis (Sequential)
#' @description Computes ASP-gene correlation across one or multiple samples using Spearman correlation.
#'              Designed to run sequentially across samples, using internal parallelization within each sample.
#'
#' @param mat A numeric matrix of junction PSI values (rows: junctions, columns: cells).
#' @param exp A numeric matrix of gene expression (rows: genes, columns: cells).
#' @param met A data.frame with columns 'cell' and 'Patient', indicating sample identity.
#' @param cell_cutoff Minimum number of cells required (default = 21).
#' @param inner_cores Number of cores to use inside each sample's computation (default = 2).
#'
#' @return A named list of data.frames per sample with ASP-gene correlations.
#' @export
multi <- function(mat, exp, met, outdir, cell_cutoff = 21 , gene_cutoff = 6, parallel = TRUE, inner_cores = 2) {
  stopifnot(is.matrix(exp), is.matrix(mat), is.data.frame(met))
  stopifnot(all(c("cell", "Patient") %in% colnames(met)))
  
  patient_ids <- unique(met$Patient)
  res_list <- list()
  
  for (pid in patient_ids) {
    met_sub <- subset(met, Patient == pid)
    cells_use <- intersect(met_sub$cell, intersect(colnames(exp), colnames(mat)))
    if (length(cells_use) < cell_cutoff) {
      warning(paste("Skipping", pid, ": too few cells"))
      res_list[[pid]] <- NULL
      next
    }
    exp_sub <- exp[, cells_use, drop = FALSE]
    mat_sub <- mat[, cells_use, drop = FALSE]
    message("Running sample: ", pid)
    res <- junction_gene_spearman(mat_sub, exp_sub, min_obs = cell_cutoff, min_fea = gene_cutoff, parallel = TRUE, n_cores = inner_cores)
    saveRDS(res,file = paste0(outdir,"/step2.",pid,".rds"))
  }
}
