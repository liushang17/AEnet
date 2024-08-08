library(Seurat)
library(ggplot2)
library(umap)

#' @title ASP Clusters enrichment Score
#' @description The enrichment score of ASP Clusters for each cell
#'
#' @param asp Cell-junction Count Matrix with Cell as columns and junction as rows
#' @param ann_junc The ASP Clusters
#'
#' @return Cell-ASP Clusters enrichment Score matrix
#' @export
#'
#' @examples asp_mat <- asp_score(asp,ann_junc)


asp_score <- function(asp,ann_junc){
  ann_junc$pattern <- ann_junc$PC
  pos <- which(ann_junc$pattern %in% "Other")
  if(length(pos) > 0){ann_junc <- ann_junc[-pos,]}

  mitall <- NULL
  clus <- unique(ann_junc$pattern)

  for(n in 1:length(clus)){
    pos <- which(ann_junc$pattern %in% clus[n])
    ann_junc1 <- ann_junc[pos,]

    mit <- NULL
    for(i in 1:nrow(ann_junc1)){
      juncs <- strsplit(ann_junc1$symbol[i],split = "\\.")[[1]]
      pos <- which(rownames(asp) %in% juncs)
      if(length(pos) == 2){
        asp1 <- asp[as.character(juncs),]
        tmp <- colSums(asp1)
        pos <- which(tmp > 10)
        asp1 <- asp1[,pos]
        tmp <- colSums(asp1)
        psi1 <- as.numeric((as.matrix(asp1[1,]))) / tmp
        psi2 <- -as.numeric((as.matrix(asp1[2,]))) / tmp
        pos <- which(psi1 == 0)
        psi1[pos] <- psi2[pos]
        asp2 <- data.frame(cell = colnames(asp1),junction = ann_junc1$symbol[i],value = psi1)
        mit <- rbind(mit,asp2)
      }

    }
    colsToKeep <- "value"
    med=mit %>%
      group_by(cell) %>%
      summarise_at(vars(colsToKeep), mean)
    med <- data.frame(med)
    colnames(med) <- c("Cell",clus[n])

    if(n == 1){
      mitall <- med
    }else{
      mitall <- merge(mitall,med,by = "Cell")
    }
  }

  return(mitall)
}

#' @title Cell Clusters
#' @description The Cell Clusters determined by ASP Clusters
#'
#' @param asp_score Cell-ASP Clusters enrichment Score matrix
#' @param resolution The resolution for the function FindClusters
#' @param min.dist The min distance for the function RunUMAP
#'
#' @return Cell clusters and the umap determined by ASP
#' @export
#'
#' @examples cell_info <- cell_clus(asp_score,resolution = 0.5,min.dist = 1)

cell_clus <- function(asp_score,resolution = 0.5,min.dist = 1){
  med4 <- asp_score[,2:ncol(asp_score)]
  rownames(med4) <- asp_score$Cell
  aspmat <- med4

  aspmat <- scale(aspmat)
  tms <- FindNeighbors(aspmat)
  tms1 <- FindClusters(tms$snn,resolution = resolution)
  colnames(tms1) <- "cluster"

  tsne_out <- RunUMAP((aspmat),min.dist = min.dist)
  #tsne_out <- RunTSNE((aspmat))
  umapinfo <- tsne_out@cell.embeddings
  tms1$umap_1 <- umapinfo[,1]
  tms1$umap_2 <- umapinfo[,2]

  return(tms1)
}


#' @title The key SF
#' @description The key splicing factors for the specific clusters
#'
#' @param corm AEN network
#' @param ann_junc The ASP Clusters
#' @param cluster The specific clusters
#'
#' @return The key splicing factors
#' @export
#'
#' @examples key_sf <- key_sf(corm,sf,ann_junc,cluster = "C_1")

key_sf <- function(corm,sf,ann_junc,cluster = "C_1"){
  pos <- which(ann_junc$PC %in% cluster)
  ann_junc1 <- ann_junc[pos,]
  pos <- which(corm$symbol %in% sf & corm$junction %in% ann_junc1$symbol)
  corm1 <- corm[pos,]
  sui <- data.frame(table(corm1$symbol))
  sui$pro <- sui$Freq / nrow(ann_junc1)
  sui <- sui[order(-sui$Freq),]
  sui$order <- 1:nrow(sui)

  sui$cor <- 0
  for(i in 1:nrow(sui)){
    pos <- which(corm$symbol %in% sui$Var1[i] & corm$junction %in% ann_junc1$symbol)
    sui$cor[i] <- mean(corm$correlation[pos])
  }

  return(sui)
}

