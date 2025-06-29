#' @title Process Correlation Results
#' @description Process and filter correlation results from multiple correlation result files
#' @param indir Input directory containing RDS files
#' @param min_freq Minimum frequency threshold for filtering (default = 2)
#' @return A processed data frame with aggregated and filtered correlation results
#' @export
#' @importFrom dplyr group_by summarise mutate select %>%
process_cor_results <- function(indir, min_freq = 2) {
  if (!dir.exists(indir)) stop("Input directory does not exist.")
  
  files <- list.files(indir, pattern = "step2", full.names = TRUE)
  if (length(files) == 0) stop("No matching 'step2' files found.")
  
  # Define a function to read and preprocess a single file
  process_file <- function(f) {
    df <- readRDS(f)
    df <- df[!is.na(df$geneID), ]
    if (nrow(df) == 0) return(NULL)
    df$alljunction <- df$junctionID
    df$type1 <- paste0(df$alljunction, ":", df$geneID)
    sign <- ifelse(df$correlation < 0, "-", "+")
    df$all_type <- paste0(df$type1, ":", sign)
    return(df)
  }
  
  if (length(files) == 1) {
    result <- process_file(files[1])
    if (is.null(result)) stop("File has no valid correlation entries.")
    result$junction <- result$alljunction
    result$symbol <- result$geneID
    result$Freq <- 1
    result$type1 <- paste0(result$alljunction, ":", result$geneID)
    final_result <- result[, c("junction", "symbol", "correlation", "type1", "Freq")]
    final_result$p.value <- result$p.value
    
  } else {
    # Read and combine all results
    cor_list <- lapply(files, process_file)
    cor_list <- Filter(Negate(is.null), cor_list)
    if (length(cor_list) == 0) stop("No valid correlation data found.")
    
    all_cor <- do.call(rbind, cor_list)
    
    # Filter by min_freq on all_type
    all_type_freq <- table(all_cor$all_type)
    keep_all_type <- names(all_type_freq)[all_type_freq >= min_freq]
    filtered <- all_cor[all_cor$all_type %in% keep_all_type, ]
    
    # Remove conflicting direction signs (+/-) per ASP-gene pair (type1)
    filtered1 <- unique(filtered[,c("all_type","type1")])
    direction_counts <- table(filtered1$type1)
    duplicated_pairs <- names(direction_counts[direction_counts > 1])
    filtered <- filtered[!filtered$type1 %in% duplicated_pairs, ]
    
    # Aggregate results by ASP-gene pair
    summary_df <- filtered %>%
      group_by(type1) %>%
      summarise(
        correlation = mean(correlation, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        junction = sub(":.*", "", type1),
        symbol = sub(".*:", "", type1)
      ) %>%
      select(junction, symbol, correlation, type1)
    
    # Add frequency info
    freq_df <- as.data.frame(table(filtered$type1))
    colnames(freq_df) <- c("type1", "Freq")
    
    final_result <- merge(summary_df, freq_df, by = "type1")
    final_result <- final_result[, c("junction", "symbol", "correlation", "type1", "Freq")]
    final_result$p.value <- 0
  }
  
  return(final_result)
}
