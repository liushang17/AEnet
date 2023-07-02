library(pheatmap)
library(amap)

pattern_analysis <- function(datatmp,cell_cutoff,cor_cutoff){
  cluster_info <- list()
  mat_used <- datatmp
  
  ##filter
  junction_pro <- matrix(nrow = length(rownames(mat_used)),ncol =2)
  for(i in 1:length(rownames(mat_used))){
    tmp <- as.numeric(as.character(as.matrix(mat_used[i,])))
    pos <- which(tmp > 1)
    junction_pro[i,1] <- rownames(mat_used)[i]
    junction_pro[i,2] <- length(pos)  
  }
  
  maxnum <- max(as.numeric(as.character(junction_pro[,2])))
  junction_pro <- data.frame(junction_pro)
  junction_pro$X3 <- as.numeric(as.character(junction_pro$X2)) / maxnum
  pos <- which(junction_pro$X3 > 0.05)
  junction_pro_filter <- junction_pro[pos,]
  pos <- which(rownames(mat_used) %in% junction_pro_filter$X1)
  mat_used <- mat_used[pos,]
  
  if(length(rownames(mat_used)) <= 1){
    return(cluster_info)
  }else{
    
    ##zero to nonzero in above matrix
    for(i in 1:length(rownames(mat_used))){
      tmp <- 1 / (i+2) / 1000
      mat_used[i,] <- mat_used[i,] * i / 2 + tmp
    }
    
    ### Clustering cells based on correlation with Hierarchical Clustering
    cell_dist <- Dist(t(log2(as.matrix(mat_used)+1)),method = "correlation")
    cell_clus <- hclust(cell_dist)
    
    ###split the cells in extremely high similarity
    cor_cutoff <- cor_cutoff
    cell_clus_info <- data.frame(cutree(cell_clus,h=cor_cutoff))
    cell_clus_info <- data.frame(cell = rownames(cell_clus_info),cluster = cell_clus_info[,1])
    
    ###find the junction pattern of the above cluster
    clus <- unique(as.character(cell_clus_info$cluster))
    junction_pattern_cluster <- matrix(nrow = length(clus),ncol = 2)
    
    for(i in 1:length(clus)){
      juncs <- "chr"
      pos <- which(cell_clus_info$cluster %in% clus[i])
      if(length(pos) == 1){
        cellname <- as.character(cell_clus_info$cell)[pos]
        pos <- which(colnames(mat_used) %in% cellname)
        tmp <- as.numeric(as.character(as.matrix(mat_used[,pos])))
        pos <- which(tmp > 1)
        juns <- rownames(mat_used)[pos]
        juns <- paste(juns,collapse = ",")
        juncs <- paste(juncs,juns,sep = ",")
      }else{
      
        celltmp <- cell_clus_info$cell[pos]
        pos <- which(colnames(mat_used) %in% celltmp)
        mattmp <- mat_used[,pos]
      
        cutofftmp <- length(colnames(mattmp)) * 0.5
        for(j in 1:length(rownames(mattmp))){
          rowtmp <- as.numeric(as.character(as.matrix(mattmp[j,])))
          pos <- which(rowtmp > 1 * j /2)
          if(length(pos) > cutofftmp){juncs <- paste(juncs,rownames(mattmp)[j],sep = ",")}
        }
      }  
        junction_pattern_cluster[i,1] <- clus[i]
        junction_pattern_cluster[i,2] <- juncs
    }
    junction_pattern_cluster <- data.frame(junction_pattern_cluster)
    
    
    ### merge the clusters with the same pattern
    junction_pattern_clusters <- unique(as.character(junction_pattern_cluster$X2))
    cell_clus_info$cluster_merge <- 0
    junction_pattern_number <- matrix(nrow = length(junction_pattern_clusters),ncol = 3)
    for(i in 1:length(junction_pattern_clusters)){
      pos <- which(as.character(junction_pattern_cluster$X2) %in% junction_pattern_clusters[i])
      junction_pattern_tmp <- junction_pattern_cluster[pos,]
      pos <- which(cell_clus_info$cluster %in% as.character(junction_pattern_tmp$X1))
      cell_clus_info$cluster_merge[pos] <- i
      junction_pattern_number[i,1] <- junction_pattern_clusters[i]
      junction_pattern_number[i,2] <- length(pos)
      junction_pattern_number[i,3] <- i
    }
    
    ##find the final pattern
    junction_pattern_number <- data.frame(junction_pattern_number)
    pos <- which(as.numeric(as.character(junction_pattern_number$X2)) > cell_cutoff )
    junction_pattern_number_filter <- junction_pattern_number[pos,]
    pos <- which(junction_pattern_number_filter$X1 %in% "chr")
    if(length(pos) > 0){junction_pattern_number_filter <- junction_pattern_number_filter[-pos,]}
    if(length(rownames(junction_pattern_number_filter)) == 0){
      return (cluster_info)
      if(length(rownames(junction_pattern_number)) >= 2){
        print("not enough cells to detect the junction pattern of each cluster")
      }else{
        print("not enough reads to detect the junction pattern of each cluster")
      }
    }else{
      junction_patterns <- NULL
      for(i in 1:length(rownames(junction_pattern_number_filter))){
        tmp <- strsplit(as.character(junction_pattern_number_filter$X1)[i],split = ",")[[1]]
        junction_patterns <- c(junction_patterns,tmp)
      }
      juns <- unique(junction_patterns)
      pos <- which(juns %in% "chr")
      juns <- juns[-pos]
    
      if(length(juns) == 0){
        return (cluster_info)
      }else{
        if(length(juns) == 1){
          final_pattern <- data.frame(clu = 1,jun = juns)
          final_pattern <- list(final_pattern)
          cell_clus_info <- list(cell_clus_info)
          cluster_info <- c(final_pattern,cell_clus_info)
          return(cluster_info)
        }else{
          juns_stat <- matrix(nrow = length(juns),ncol = length(rownames(junction_pattern_number_filter)))
          rownames(juns_stat) <- juns
          colnames(juns_stat) <- junction_pattern_number_filter$X3
          for(i in 1:length(rownames(junction_pattern_number_filter))){
            tmp <- strsplit(as.character(junction_pattern_number_filter$X1)[i],split = ",")[[1]]
            pos <- which(rownames(juns_stat) %in% tmp)
            juns_stat[pos,i] <- 1
            juns_stat[-pos,i] <- 0
          }
        
          #pheatmap(juns_stat)
          juns_stat_dis <- dist(juns_stat,method = "euclidean")
          juns_stat_hier <- hclust(juns_stat_dis)
          #plclust(sui6)
          juns_stat_cluster <- data.frame(cutree(juns_stat_hier,h=0))
          final_pattern <- data.frame(clu = juns_stat_cluster[,1],jun = rownames(juns_stat_cluster))
          final_pattern <- list(final_pattern)
          cell_clus_info <- list(cell_clus_info)
          cluster_info <- c(final_pattern,cell_clus_info)
          return(cluster_info)
        }
      }
    }
  }
}
