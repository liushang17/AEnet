library(pheatmap)
library("reshape2")
library("tidyr")
library(data.table)
library(Matrix)
library(umap)
library(ggsci)
library(ggplot2)

#' @title AEN filter
#' @description The function to filter the non-key links with AEN network
#'
#' @param corm The data frame contating the relationship (AEN) between ASP and Gene from the above step 2 or 1.
#' @param sample_cutoff The number cutoff of samples for the links within the AEN
#' @param gene_cutoff The top key ASPs used to the downstream analysis
#' @param link_cutoff the link cutoff for the ASPs
#'
#' @return The high quality AEN
#' @export
#'
#' @examples AEN_h <- asp_selction(corm,sample_filter = T,sample_cutoff = 2,gene_cutoff = 1500,link_cutoff = 10)


asp_selection <- function(corm,sample_filter = T,sample_cutoff = 2,gene_cutoff = 1500,link_cutoff = 10){

  ############# filter link ##############
  pos <- which(corm$Freq > sample_cutoff)
  corm2 <- corm[pos,]

  sui <- data.frame(table(corm2$junction))
  sui <- sui[order(-sui$Freq),]
  sui1 <- sui[1:gene_cutoff,]
  pos <- which(sui1$Freq >= link_cutoff)
  sui1 <- sui1[pos,]
  pos <- which(corm2$junction %in% sui1$Var1)
  corm2 <- corm2[pos,]

  sui$order <- 1:nrow(sui)
  sui$type <- "blue"
  pos <- which(sui$Var1 %in% sui1$Var1)
  sui$type[pos] <- "red"
  sui <- sui[1:5000,]
  sui$type <- factor(sui$type,levels = c("red","blue"))
  #p1 <- ggplot(sui,aes(x=order,y=Freq,color = type))+geom_point()+theme_classic()
  #print(p1)

  return(corm2)
}

#' @title AEN Clustering
#' @description The function to detect asp clusters based on the AEN network
#'
#' @param corm2 The data frame contating the high quality AEN network
#' @param cluster_num The number of clusters
#' @param asp_num The least number of ASPs within ASP Clusters
#'
#' @return a data frame containing the information of ASP clusters
#' @export
#'
#' @examples junction_clusters <- junction_clustering(corm2,cluster_num = 25,asp_num = 10)


junction_clustering <- function(corm2,cluster_num = 25,asp_num = 10){
  ############# matrix transformation #########
  corm2$FC <- 1
  pos <- which(corm2$correlation < 0)
  corm2$FC[pos] <- -1

  corm2$cellID<-corm2$junction
  corm2$geneID <- corm2$symbol
  gene=unique(corm2$geneID)
  cell=unique(corm2$cellID)
  gene_idx=c(1:length(gene))
  cell_idx=c(1:length(cell))
  names(gene_idx)=gene
  names(cell_idx)=cell
  mat=sparseMatrix(i=gene_idx[corm2$geneID],j=cell_idx[corm2$cellID],x=as.numeric(corm2$FC))
  mat <- as.matrix(mat)
  rownames(mat)=gene
  colnames(mat)=cell

  ########## feature clustering : step1 ########################
  mat1 <- t(mat)
  mat1[mat1 == -1] <- 0
  mat2 <- t(mat)
  mat2[mat2 == 1] <- 0
  mat2[mat2 == -1] <- 1
  colnames(mat2) <- paste0(colnames(mat2),"_neg")
  mat3 <- cbind(mat1,mat2)

  tesm <- proxy::dist(mat3, by_rows = TRUE, method = "Jaccard")
  tesm <- as.matrix(tesm)
  tesm1 <- 1 - tesm
  p1 <- pheatmap::pheatmap(tesm1,cluster_rows = T,cluster_cols = T)

  num <- cluster_num
  rowinfo <- data.frame(cutree(p1$tree_row,k = num))
  rowinfo$name <- rownames(rowinfo)
  colnames(rowinfo) <- c("cluster","symbol")
  rowinfo$final_cluster <- rowinfo$cluster

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
  pos <- which(sui$Freq >= 10)
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
