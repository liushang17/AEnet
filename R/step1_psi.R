#' @title Fast Junction Usage Ratio Calculator
#' @description Computes junction1/(junction1+junction2) ratios with thresholding (vectorized implementation)
#' @param junc_mat Sparse or dense matrix (cells x junctions)
#' @param asp_result ASP result from asp() function
#' @param min_total_count Minimum total reads required (default: 5)
#' @return Sparse matrix of ratios (pairs x cells)
#' @export
#' @import Matrix
library(Matrix)
calculate_psi <- function(junc_mat, asp_result, outdir, min_total_count = 5) {
  
  # Convert to sparse matrix if not already
  if (!inherits(junc_mat, "sparseMatrix")) {
    junc_mat <- Matrix(junc_mat, sparse = TRUE)
  }
  
  # Create junction name to index mapping
  junc_index <- setNames(1:nrow(junc_mat), rownames(junc_mat))
  
  # Get indices for all junctions
  j1_idx <- junc_index[asp_result$junction1]
  j2_idx <- junc_index[asp_result$junction2]
  valid_pairs <- !is.na(j1_idx) & !is.na(j2_idx)
  
  # Extract counts in one operation
  j1_counts <- junc_mat[j1_idx[valid_pairs], , drop = FALSE]
  j2_counts <- junc_mat[j2_idx[valid_pairs], , drop = FALSE]
  
  # Vectorized calculations
  total_counts <- j1_counts + j2_counts
  threshold_mask <- total_counts > min_total_count
  
  # Compute ratios (preserving sparsity)
  ratios <- Matrix(0, nrow = nrow(asp_result), ncol = ncol(junc_mat),
                   dimnames = list(asp_result$all, colnames(junc_mat)),
                   sparse = TRUE)
  
  ratios[valid_pairs, ] <- j1_counts / (total_counts + (total_counts == 0)) # Avoid div/0
  ratios[!threshold_mask] <- NA
  
  saveRDS(ratios,file = paste0(outdir,"/step1.psi.rds"))
  return(ratios)
}
