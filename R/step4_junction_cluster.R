#' @title AEN filter
#' @description Filters non-key links in the AEN network based on frequency and directional bias.
#'
#' @param corm Data frame containing ASP-gene relationships with columns: junction, correlation, Freq, etc.
#' @param sample_cutoff Minimum number of samples supporting each link (default = 2)
#' @param gene_cutoff Number of top junctions to keep based on number of links (default = 1500)
#' @param link_cutoff Minimum number of biased (directionally consistent) links per junction (default = 10)
#'
#' @return Filtered high-confidence ASP-gene correlation network (AEN)
#' @export
asp_selection <- function(corm, sample_cutoff = 2, gene_cutoff = 1500, link_cutoff = 10) {
  stopifnot(all(c("junction", "correlation", "Freq") %in% colnames(corm)))
  
  # Step 1: Filter by minimum sample support
  corm_filtered <- corm[corm$Freq > sample_cutoff, ]
  
  # Step 2: Calculate + / - correlation direction counts per junction
  all_counts <- as.data.frame(table(corm_filtered$junction))
  neg_counts <- as.data.frame(table(corm_filtered$junction[corm_filtered$correlation < 0]))
  names(all_counts) <- c("junction", "total")
  names(neg_counts) <- c("junction", "neg")
  
  junction_stat <- merge(all_counts, neg_counts, by = "junction", all.x = TRUE)
  junction_stat$neg[is.na(junction_stat$neg)] <- 0
  junction_stat$pos <- junction_stat$total - junction_stat$neg
  
  # Step 3: Keep junctions with enough biased links (either + or - direction dominates)
  keep_junctions <- junction_stat$junction[
    (junction_stat$pos > link_cutoff & junction_stat$neg > link_cutoff) 
  ]
  
  corm_biased <- corm_filtered[corm_filtered$junction %in% keep_junctions, ]
  
  # Step 4: Keep top N junctions with at least link_cutoff links
  junction_freq <- sort(table(corm_biased$junction), decreasing = TRUE)
  top_junctions <- names(junction_freq)[
    seq_len(min(gene_cutoff, sum(junction_freq >= link_cutoff)))
  ]
  top_junctions <- top_junctions[junction_freq[top_junctions] >= link_cutoff]
  
  corm_final <- corm_biased[corm_biased$junction %in% top_junctions, ]
  return(corm_final)
}

#' @title Gene Selection Based on Correlation and Annotation
#' @description Filters and selects genes based on correlation metrics and annotation data
#' @param corm Data frame containing correlation results with required columns: junction, correlation, Freq
#' @param aspann Data frame containing junction-gene annotations with columns: junction1, geneud, genename
#' @param sample_cutoff Minimum sample frequency threshold (default=2)
#' @param gene_cutoff Maximum number of genes to consider (default=1500) - currently not used
#' @param link_cutoff Minimum number of gene links required per junction (default=10)
#' @return Filtered data frame containing selected gene correlations with annotation
#' @export
#' @importFrom dplyr rename select left_join distinct count filter pull
library(dplyr)
gene_selection <- function(corm, aspann, sample_cutoff = 2, gene_cutoff = 1500, link_cutoff = 30) {
  # Validate input structure
  required_cols <- c("junction", "correlation", "Freq")
  if (!all(required_cols %in% colnames(corm))) {
    stop("Input dataframe 'corm' must contain columns: ", paste(required_cols, collapse = ", "))
  }
  
  # Step 1: Filter by sample frequency threshold
  corm_filtered <- corm[corm$Freq > sample_cutoff, ]
  
  # Swap junction and symbol columns using safe renaming
  corm_filtered$tmp <- corm_filtered$junction
  corm_filtered <- corm_filtered %>%
    rename(
      tmp = junction,
      junction = symbol,
      symbol = tmp
    ) %>%
    select(-tmp)
  
  # Step 2: Prepare and merge annotation data
  colnames(aspann) <- c("junction1", "geneud", "genename")
  
  # Perform left join to preserve row order
  corm_filtered$junction1 <- sub("\\..*", "", corm_filtered$symbol)
  merged <- merge(corm_filtered, aspann, by = "junction1")
  
  gene_asp_table <- unique(merged[, c("junction", "genename")])
  gene_freq <- as.data.frame(table(gene_asp_table$junction))
  colnames(gene_freq) <- c("gene", "n_asp")
  
  # Step 5: Keep genes with sufficient ASP links
  key_genes <- gene_freq$gene[gene_freq$n_asp > link_cutoff]
  corm_final <- merged[merged$junction %in% key_genes, ]
  
  return(corm_final)
}


#' @title AEN Clustering
#' @description The function to detect asp clusters based on the AEN network
#'
#' @param corm2 The data frame contating the high quality AEN network
#' @param cluster_num The number of clusters
#' @param asp_num The least number of ASPs within ASP Clusters
#'
#' @return a list containing a data frame of ASP clusters and a similarity matrix
#' @export
junction_clustering <- function(corm2, cluster_num = 25, asp_num = 10) {
  corm2$FC <- ifelse(corm2$correlation < 0, -1, 1)
  
  gene <- unique(corm2$symbol)
  asp <- unique(corm2$junction)
  
  gene_idx <- setNames(seq_along(gene), gene)
  asp_idx <- setNames(seq_along(asp), asp)
  
  mat <- sparseMatrix(i = gene_idx[corm2$symbol], j = asp_idx[corm2$junction], x = corm2$FC)
  mat <- as.matrix(mat)
  rownames(mat) <- gene
  colnames(mat) <- asp
  
  # Expand positive and negative profiles
  pos_mat <- t(mat)
  neg_mat <- t(mat)
  pos_mat[pos_mat == -1] <- 0
  neg_mat[neg_mat == 1] <- 0
  neg_mat[neg_mat == -1] <- 1
  colnames(neg_mat) <- paste0(colnames(neg_mat), "_neg")
  combined_mat <- cbind(pos_mat, neg_mat)
  
  # Similarity
  dist_mat <- proxy::dist(combined_mat, method = "Jaccard")
  sim_mat <- 1 - as.matrix(dist_mat)
  
  cluster_mat = function(mat, distance, method){
    if(!(method %in% c("ward.D", "ward.D2", "ward", "single", "complete", "average", "mcquitty", "median", "centroid"))){
      stop("clustering method has to one form the list: 'ward', 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median' or 'centroid'.")
    }
    if(!(distance[1] %in% c("correlation", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")) & !(inherits(distance, "dist"))){
      stop("distance has to be a dissimilarity structure as produced by dist or one measure  form the list: 'correlation', 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary', 'minkowski'")
    }
    if(distance[1] == "correlation"){
      d = as.dist(1 - cor(t(mat)))
    }
    else{
      if(inherits(distance, "dist")){
        d = distance
      }
      else{
        d = dist(mat, method = distance)
      }
    }
    
    return(hclust(d, method = method))
  }
  
  row_order <- cluster_mat(sim_mat,distance = "euclidean",method = "complete")
  clusters <- cutree(row_order, k = cluster_num)
  
  rowinfo <- data.frame(symbol = rownames(sim_mat), final_cluster = as.character(clusters))
  
  tesm1 <- sim_mat
  
  rowinfo <- rowinfo[order(as.numeric(rowinfo$final_cluster)),]
  rowinfo1 <- data.frame(clus = as.character(rowinfo$final_cluster))
  rownames(rowinfo1) <- rowinfo$symbol
  tesm2 <- tesm1[rownames(rowinfo1),rownames(rowinfo1)]
  
  ############### filter ##########################################
  clus <- unique(as.character(rowinfo$final_cluster))
  mit <- matrix(nrow = length(clus),ncol = 4)
  for(i in 1:length(clus)){
    pos <- which(rowinfo$final_cluster %in% clus[i])
    rowinfo1 <- rowinfo[pos,]
    tesm3 <- tesm1[rowinfo1$symbol,rowinfo1$symbol]
    tesm3[tesm3 > 0.1] <- 1
    tesm3[tesm3 < 0.1] <- 0
    tmp <- (sum(rowSums(tesm3))  - nrow(rowinfo1))
    mit[i,1] <- clus[i]
    mit[i,2] <- tmp / (nrow(rowinfo1) ) / (nrow(rowinfo1) - 1)
    
    tesm3 <- tesm1[rowinfo1$symbol,rowinfo1$symbol]
    tesm3[tesm3 > 0.05] <- 1
    tesm3[tesm3 < 0.05] <- 0
    tmp <- (sum(rowSums(tesm3)) - nrow(rowinfo1))
    mit[i,3] <- tmp / (nrow(rowinfo1) ) / (nrow(rowinfo1) - 1)
    
    tesm3 <- tesm1[rowinfo1$symbol,rowinfo1$symbol]
    tesm3[tesm3 > 0.01] <- 1
    tesm3[tesm3 < 0.01] <- 0
    tmp <- (sum(rowSums(tesm3)) - nrow(rowinfo1))
    mit[i,4] <- tmp / (nrow(rowinfo1) ) / (nrow(rowinfo1) - 1)
  }
  pos <- which(as.numeric(as.character(mit[,2])) > 0.1)
  mit1 <- mit[pos,]
  pos <- which(rowinfo$final_cluster %in% mit1[,1])
  rowinfo1 <- rowinfo[pos,]
  
  rowinfo1 <- rowinfo1[order(as.numeric(rowinfo1$final_cluster)),]
  rowinfo2 <- data.frame(clus = as.character(rowinfo1$final_cluster))
  rownames(rowinfo2) <- rowinfo1$symbol
  tesm2 <- tesm1[rownames(rowinfo2),rownames(rowinfo2)]
  
  ############### merge ##########################################
  sui <- data.frame(table(rowinfo1$final_cluster))
  pos <- which(sui$Freq < asp_num)
  sui1 <- sui[pos,]
  pos <- which(sui$Freq >= asp_num)
  sui2 <- sui[pos,]
  sui1$clu <- 0
  sui1$pro <- 0
  for(i in 1:nrow(sui1)){
    pos <- which(rowinfo1$final_cluster %in% sui1$Var1[i])
    rowinfo3 <- rowinfo1[pos,]
    sui2$siz <- 0
    for(j in 1:nrow(sui2)){
      pos <- which(rowinfo1$final_cluster %in% sui2$Var1[j])
      rowinfo4 <- rowinfo1[pos,]
      tesm3 <- tesm1[rowinfo3$symbol,rowinfo4$symbol]
      tesm3[tesm3 > 0.1] <- 1
      tesm3[tesm3 < 0.1] <- 0
      tmp <- (sum(rowSums(tesm3)))
      sui2$siz[j] <- tmp / nrow(rowinfo3) / nrow(rowinfo4)
    }
    sui3 <- sui2
    pos <- which(sui3$Var1 == sui1$Var1[i])
    if(length(pos) > 0){sui3 <- sui3[-pos,]}
    pos <- which(as.numeric(as.character(sui3$siz)) == max(as.numeric(as.character(sui3$siz))))
    sui1$clu[i] <- as.character(sui3$Var1)[pos[1]]
    sui1$pro[i] <- as.numeric(as.character(sui3$siz))[pos[1]]
  }
  
  rowinfo1$pattern <- rowinfo1$final_cluster
  for(i in 1:nrow(sui1)){
    pos <- which(rowinfo1$final_cluster %in% sui1$Var1[i])
    if(sui1$pro[i] > 0.1){
      rowinfo1$pattern[pos] <- sui1$clu[i]
    }
  }
  
  rowinfo1$final_cluster <- rowinfo1$pattern
  
  sui <- data.frame(table(rowinfo1$final_cluster))
  pos <- which(sui$Freq <= asp_num)
  sui1 <- data.frame(sui[pos,])
  pos <- which(sui$Freq > asp_num)
  sui2 <- sui[pos,]
  if(nrow(sui1) > 0){
    sui1$clu <- 0
    sui1$pro <- 0
    for(i in 1:nrow(sui1)){
      pos <- which(rowinfo1$final_cluster %in% sui1$Var1[i])
      rowinfo3 <- rowinfo1[pos,]
      sui2$siz <- 0
      for(j in 1:nrow(sui2)){
        pos <- which(rowinfo1$final_cluster %in% sui2$Var1[j])
        rowinfo4 <- rowinfo1[pos,]
        tesm3 <- tesm1[rowinfo3$symbol,rowinfo4$symbol]
        tesm3[tesm3 > 0.1] <- 1
        tesm3[tesm3 < 0.1] <- 0
        tmp <- (sum(rowSums(tesm3)))
        sui2$siz[j] <- tmp / nrow(rowinfo3) / nrow(rowinfo4)
      }
      sui3 <- sui2
      pos <- which(sui3$Var1 == sui1$Var1[i])
      if(length(pos) > 0){sui3 <- sui3[-pos,]}
      pos <- which(as.numeric(as.character(sui3$siz)) == max(as.numeric(as.character(sui3$siz))))
      sui1$clu[i] <- as.character(sui3$Var1)[pos[1]]
      sui1$pro[i] <- as.numeric(as.character(sui3$siz))[pos[1]]
    }
    
    rowinfo1$pattern <- rowinfo1$final_cluster
    for(i in 1:nrow(sui1)){
      pos <- which(rowinfo1$final_cluster %in% sui1$Var1[i])
      if(sui1$pro[i] > 0.1){
        rowinfo1$pattern[pos] <- sui1$clu[i]
      }
    }
  }else{
    rowinfo1$pattern <- rowinfo1$final_cluster
  }
  
  
  sui <- data.frame(table(rowinfo1$pattern))
  pos <- which(sui$Freq >= asp_num)
  sui2 <- sui[pos,]
  
  rowinfo$PC <- "Other"
  for(i in 1:nrow(sui2)){
    pos <- which(rowinfo1$pattern %in% sui2$Var1[i])
    rowinfo2 <- rowinfo1[pos,]
    pos <- which(rowinfo$symbol %in% rowinfo2$symbol)
    rowinfo$PC[pos] <- paste0("C_",i)
  }
  resall <- list(asp_clusters = rowinfo,asp_simm = tesm2)
  return(resall)
}

#' @title AEN Clustering
#' @description The function to detect asp clusters based on the AEN network
#'
#' @param corm2 The data frame contating the high quality AEN network with junction as symbol
#' @param cluster_num The number of clusters
#' @param asp_num The least number of ASPs within ASP Clusters
#'
#' @return a list containing a data frame of ASP clusters and a similarity matrix
#' @export
gene_clustering <- function(corm2, cluster_num = 25, asp_num = 10) {
  corm2$FC <- ifelse(corm2$correlation < 0, -1, 1)
  
  gene <- unique(corm2$symbol)
  asp <- unique(corm2$junction)
  
  gene_idx <- setNames(seq_along(gene), gene)
  asp_idx <- setNames(seq_along(asp), asp)
  
  mat <- sparseMatrix(i = gene_idx[corm2$symbol], j = asp_idx[corm2$junction], x = corm2$FC)
  mat <- as.matrix(mat)
  rownames(mat) <- gene
  colnames(mat) <- asp
  
  # Expand positive and negative profiles
  pos_mat <- t(mat)
  neg_mat <- t(mat)
  pos_mat[pos_mat == -1] <- 0
  neg_mat[neg_mat == 1] <- 0
  neg_mat[neg_mat == -1] <- 1
  colnames(neg_mat) <- paste0(colnames(neg_mat), "_neg")
  combined_mat <- cbind(pos_mat, neg_mat)
  
  # Similarity
  dist_mat <- proxy::dist(combined_mat, method = "Jaccard")
  sim_mat <- 1 - as.matrix(dist_mat)
  
  cluster_mat = function(mat, distance, method){
    if(!(method %in% c("ward.D", "ward.D2", "ward", "single", "complete", "average", "mcquitty", "median", "centroid"))){
      stop("clustering method has to one form the list: 'ward', 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median' or 'centroid'.")
    }
    if(!(distance[1] %in% c("correlation", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")) & !(inherits(distance, "dist"))){
      stop("distance has to be a dissimilarity structure as produced by dist or one measure  form the list: 'correlation', 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary', 'minkowski'")
    }
    if(distance[1] == "correlation"){
      d = as.dist(1 - cor(t(mat)))
    }
    else{
      if(inherits(distance, "dist")){
        d = distance
      }
      else{
        d = dist(mat, method = distance)
      }
    }
    
    return(hclust(d, method = method))
  }
  
  row_order <- cluster_mat(sim_mat,distance = "euclidean",method = "complete")
  clusters <- cutree(row_order, k = cluster_num)
  
  rowinfo <- data.frame(symbol = rownames(sim_mat), final_cluster = as.character(clusters))
  
  tesm1 <- sim_mat
  
  rowinfo <- rowinfo[order(as.numeric(rowinfo$final_cluster)),]
  rowinfo1 <- data.frame(clus = as.character(rowinfo$final_cluster))
  rownames(rowinfo1) <- rowinfo$symbol
  tesm2 <- tesm1[rownames(rowinfo1),rownames(rowinfo1)]
  
  ############### filter ##########################################
  clus <- unique(as.character(rowinfo$final_cluster))
  mit <- matrix(nrow = length(clus),ncol = 4)
  for(i in 1:length(clus)){
    pos <- which(rowinfo$final_cluster %in% clus[i])
    rowinfo1 <- rowinfo[pos,]
    tesm3 <- tesm1[rowinfo1$symbol,rowinfo1$symbol]
    tesm3[tesm3 > 0.1] <- 1
    tesm3[tesm3 < 0.1] <- 0
    tmp <- (sum(rowSums(tesm3))  - nrow(rowinfo1))
    mit[i,1] <- clus[i]
    mit[i,2] <- tmp / (nrow(rowinfo1) ) / (nrow(rowinfo1) - 1)
    
    tesm3 <- tesm1[rowinfo1$symbol,rowinfo1$symbol]
    tesm3[tesm3 > 0.05] <- 1
    tesm3[tesm3 < 0.05] <- 0
    tmp <- (sum(rowSums(tesm3)) - nrow(rowinfo1))
    mit[i,3] <- tmp / (nrow(rowinfo1) ) / (nrow(rowinfo1) - 1)
    
    tesm3 <- tesm1[rowinfo1$symbol,rowinfo1$symbol]
    tesm3[tesm3 > 0.01] <- 1
    tesm3[tesm3 < 0.01] <- 0
    tmp <- (sum(rowSums(tesm3)) - nrow(rowinfo1))
    mit[i,4] <- tmp / (nrow(rowinfo1) ) / (nrow(rowinfo1) - 1)
  }
  pos <- which(as.numeric(as.character(mit[,2])) > 0.1)
  mit1 <- mit[pos,]
  pos <- which(rowinfo$final_cluster %in% mit1[,1])
  rowinfo1 <- rowinfo[pos,]
  
  rowinfo1 <- rowinfo1[order(as.numeric(rowinfo1$final_cluster)),]
  rowinfo2 <- data.frame(clus = as.character(rowinfo1$final_cluster))
  rownames(rowinfo2) <- rowinfo1$symbol
  tesm2 <- tesm1[rownames(rowinfo2),rownames(rowinfo2)]
  
  ############### merge ##########################################
  sui <- data.frame(table(rowinfo1$final_cluster))
  pos <- which(sui$Freq == 1)
  if(length(pos) > 0){sui <- sui[-pos,]}
  
  pos <- which(sui$Freq < asp_num)
  sui1 <- sui[pos,]
  pos <- which(sui$Freq >= asp_num)
  sui2 <- sui[pos,]
  if(nrow(sui1) > 0){
    sui1$clu <- 0
    sui1$pro <- 0
    for(i in 1:nrow(sui1)){
      pos <- which(rowinfo1$final_cluster %in% sui1$Var1[i])
      rowinfo3 <- rowinfo1[pos,]
      sui2$siz <- 0
      for(j in 1:nrow(sui2)){
        pos <- which(rowinfo1$final_cluster %in% sui2$Var1[j])
        rowinfo4 <- rowinfo1[pos,]
        tesm3 <- tesm1[rowinfo3$symbol,rowinfo4$symbol]
        tesm3[tesm3 > 0.1] <- 1
        tesm3[tesm3 < 0.1] <- 0
        tmp <- (sum(rowSums(tesm3)))
        sui2$siz[j] <- tmp / nrow(rowinfo3) / nrow(rowinfo4)
      }
      sui3 <- sui2
      pos <- which(sui3$Var1 == sui1$Var1[i])
      if(length(pos) > 0){sui3 <- sui3[-pos,]}
      pos <- which(as.numeric(as.character(sui3$siz)) == max(as.numeric(as.character(sui3$siz))))
      sui1$clu[i] <- as.character(sui3$Var1)[pos[1]]
      sui1$pro[i] <- as.numeric(as.character(sui3$siz))[pos[1]]
    }
    
    rowinfo1$pattern <- rowinfo1$final_cluster
    for(i in 1:nrow(sui1)){
      pos <- which(rowinfo1$final_cluster %in% sui1$Var1[i])
      if(sui1$pro[i] > 0.1){
        rowinfo1$pattern[pos] <- sui1$clu[i]
      }
    }
    rowinfo1$final_cluster <- rowinfo1$pattern
  }
  
  
  sui <- data.frame(table(rowinfo1$final_cluster))
  pos <- which(sui$Freq <= asp_num)
  sui1 <- data.frame(sui[pos,])
  pos <- which(sui$Freq > asp_num)
  sui2 <- sui[pos,]
  if(nrow(sui1) > 0){
    sui1$clu <- 0
    sui1$pro <- 0
    for(i in 1:nrow(sui1)){
      pos <- which(rowinfo1$final_cluster %in% sui1$Var1[i])
      rowinfo3 <- rowinfo1[pos,]
      sui2$siz <- 0
      for(j in 1:nrow(sui2)){
        pos <- which(rowinfo1$final_cluster %in% sui2$Var1[j])
        rowinfo4 <- rowinfo1[pos,]
        tesm3 <- tesm1[rowinfo3$symbol,rowinfo4$symbol]
        tesm3[tesm3 > 0.1] <- 1
        tesm3[tesm3 < 0.1] <- 0
        tmp <- (sum(rowSums(tesm3)))
        sui2$siz[j] <- tmp / nrow(rowinfo3) / nrow(rowinfo4)
      }
      sui3 <- sui2
      pos <- which(sui3$Var1 == sui1$Var1[i])
      if(length(pos) > 0){sui3 <- sui3[-pos,]}
      pos <- which(as.numeric(as.character(sui3$siz)) == max(as.numeric(as.character(sui3$siz))))
      sui1$clu[i] <- as.character(sui3$Var1)[pos[1]]
      sui1$pro[i] <- as.numeric(as.character(sui3$siz))[pos[1]]
    }
    
    rowinfo1$pattern <- rowinfo1$final_cluster
    for(i in 1:nrow(sui1)){
      pos <- which(rowinfo1$final_cluster %in% sui1$Var1[i])
      if(sui1$pro[i] > 0.1){
        rowinfo1$pattern[pos] <- sui1$clu[i]
      }
    }
  }else{
    rowinfo1$pattern <- rowinfo1$final_cluster
  }
  
  
  sui <- data.frame(table(rowinfo1$pattern))
  pos <- which(sui$Freq >= asp_num)
  sui2 <- sui[pos,]
  
  rowinfo$PC <- "Other"
  for(i in 1:nrow(sui2)){
    pos <- which(rowinfo1$pattern %in% sui2$Var1[i])
    rowinfo2 <- rowinfo1[pos,]
    pos <- which(rowinfo$symbol %in% rowinfo2$symbol)
    rowinfo$PC[pos] <- paste0("C_",i)
  }
  resall <- list(asp_clusters = rowinfo,asp_simm = tesm2)
  return(resall)
}
