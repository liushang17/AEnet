#' @title Detect Alternative Splicing Patterns (ASP)
#'
#' @description
#' Identifies Alternative Splicing Patterns by analyzing splice junction combinations
#' that share either common start or end positions.
#'
#' @param junc_mat A matrix or data frame with junction names as row names in
#'        "chr_start_end_site" format
#' @param min_freq Minimum occurrence frequency for a site to be considered (default: 2)
#' @param n_cores Number of CPU cores for parallel processing (default: total cores - 1)
#'
#' @return A data.table with junction pairs sharing either:
#'         - Same chromosome, start position and splice site (chr_st), or
#'         - Same chromosome, end position and splice site (chr_en)
#' @export
#'
#' @examples
#' # Create example junction matrix
#' junc_mat <- matrix(rnorm(100), nrow = 4,
#'                   dimnames = list(c("chr1_100_200_A", "chr1_100_300_A",
#'                                    "chr1_150_200_A", "chr2_50_100_T"),
#'                                  paste0("sample", 1:25)))
#' asp_results <- asp(junc_mat)

asp <- function(junc_mat, outdir, min_freq = 2, n_cores = 2) {
  # Convert input to data.table format
  desj <- data.frame(V1 = rownames(junc_mat), type = "junction")
  setDT(desj)
  setnames(desj, 1, "V1")
  annj <- unique(desj)

  # Parse junction information
  # Format: "chr_start_end_site"
  annj[, c("chr", "st", "en", "site") := tstrsplit(V1, "_", fixed = TRUE)]

  # Create composite keys for start/end positions
  annj[, chr_st := paste(chr, st, site, sep = "_")]  # chr_start_site
  annj[, chr_en := paste(chr, en, site, sep = "_")]  # chr_end_site

  # Function to find junction pairs sharing common features
  get_pairs_parallel <- function(col) {
    # Filter for sites meeting frequency threshold
    sui <- annj[, .N, by = col][N >= min_freq]

    if (nrow(sui) == 0) return(NULL)

    # Parallel processing setup
    cl <- makeCluster(n_cores)
    clusterExport(cl, c("annj", "col"), envir = environment())
    clusterEvalQ(cl, library(data.table))

    # Generate all possible junction pairs for each qualifying site
    pairs_list <- parLapply(cl, sui[[col]], function(x) {
      junctions <- annj[get(col) == x, V1]
      if (length(junctions) >= 2) {
        combn(junctions, 2, simplify = FALSE)
      } else NULL
    })

    stopCluster(cl)

    # Format results
    pairs <- unlist(pairs_list, recursive = FALSE)
    if (length(pairs) > 0) {
      data.table(
        all = sapply(pairs, function(p) paste(p, collapse = ".")),
        junction1 = sapply(pairs, `[`, 1),
        junction2 = sapply(pairs, `[`, 2)
      )
    } else NULL
  }

  # Find pairs sharing either start or end positions
  mit <- rbindlist(list(
    get_pairs_parallel("chr_st"),  # Pairs with common starts
    get_pairs_parallel("chr_en")   # Pairs with common ends
  ), fill = TRUE)
  saveRDS(mit,file = paste0(outdir,"/step0.asp.rds"))
  return(mit)
}
