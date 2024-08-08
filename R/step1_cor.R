library(Rfast)

#' @title Correlation
#' @description The function to detect the correlation between ASP (Alternative Splicing Pattern) and Gene Expression Network (AEN) within one sample or multiple samples
#'
#' @param exp Cell-Gene Expression Matrix with Cell as columns and Gene as rows
#' @param mat Cell-junction Count Matrix with Cell as columns and junction as rows
#' @param annj The ASP annotation files from Step 0
#' @param met The Cell annotation files with the name of columns : cell and Patient
#' @param pa The Name of Patient or Sample chose to detect the correlation between ASP and gene
#' @param coverage_cutoff The cutoff of coverage for ASPs
#' @param dcell_cutoff The number cutoff of cells with coverage > coverage_cutoff
#' @param gene_cutoff The number cutoff for genes with the expression within cells > gene_cutoff
#'
#' @return a data frame containing the information of the relationship between ASP and gene with columns: ASP, Gene, correlation, P value, and Others.
#' @export
#'
#' @examples cor <- step1_cor(exp,mat,annj,met,pa = "ALL",coverage_cutoff = 5,cell_cutoff = 20,gene_cutoff = 5)


step1_cor <- function(exp,mat,annj,met,pa = "ALL",coverage_cutoff = 5,cell_cutoff = 20,gene_cutoff = 5){
  ############ cell selection ###############
  if(pa == "ALL"){
    met <- met
  }else{
    pos <- which(met$Patient %in% pa)
    met <- met[pos,]
  }

  pos <- which(colnames(exp) %in% met$cell)
  exp <- exp[,pos]
  pos <- which(colnames(mat) %in% met$cell)
  mat <- mat[,pos]
  pos <- which(colnames(mat) %in% colnames(exp))
  mat <- mat[,pos]
  exp <- exp[,colnames(mat)]

  ############ correlation ###############
  mitall <- NULL
  for(i in 1:nrow(annj)){
    pos <- which(rownames(mat) %in% c(annj$junction1[i],annj$junction2[i]))
    mat1 <- mat[pos,]
    tmp <- colSums(mat1)
    pos <- which(tmp > coverage_cutoff)
    len <- length(pos) / length(tmp)
    if(length(pos) > cell_cutoff & len > 0.01){
      mat1 <- mat1[,pos]
      tmp <- tmp[pos]

      psi <- as.numeric((as.matrix(mat1[1,]))) / tmp
      exp1 <- exp[,colnames(mat1)]
      exp2 <- exp1
      exp2[exp2 > 0] <- 1
      gene.freq <- rowSums(exp2)
      pos <- which(gene.freq > gene_cutoff)
      exp1 <- exp1[pos,]

      tmp2 <- correls(psi,t(exp1),type = "spearman")
      tmp2 <- data.frame(tmp2)
      tmp2$geneID <- rownames(tmp2)
      tmp2$alljunction <- annj$all[i]
      pos <- which(tmp2$p.value < 0.01)
      mit <- tmp2[pos,]

      mit <- data.frame(mit)
      if(nrow(mit) > 2){
#        mit1 <- merge(mit,geneinfo,by = "geneID")
        mit1 <- mit
        mitall <- rbind(mitall,mit1)
      }
    }
  }
  return(mitall)
}

#' @title multi
#' @description The function to detect the correlation between ASP (Alternative Splicing Pattern) and Gene Expression Network (AEN) within one sample or multiple samples
#'
#' @param exp Cell-Gene Expression Matrix with Cell as columns and Gene as rows
#' @param mat Cell-junction Count Matrix with Cell as columns and junction as rows
#' @param annj The ASP annotation files from Step 0
#' @param met The Cell annotation files with the name of columns : cell and Patient
#' @param pa The Name of Patient or Sample chose to detect the correlation between ASP and gene
#' @param coverage_cutoff The cutoff of coverage for ASPs
#' @param dcell_cutoff The number cutoff of cells with coverage > coverage_cutoff
#' @param gene_cutoff The number cutoff for genes with the expression within cells > gene_cutoff
#'
#' @return a list with multiple data frame containing the information of the relationship between ASP and gene with columns: ASP, Gene, correlation, P value, and Others.
#' @export
#'
#' @examples cor <- multi(exp,mat,annj,met,pa = "ALL",coverage_cutoff = 5,cell_cutoff = 20,gene_cutoff = 5)

mutli <- function(exp,mat,annj,met,pa = "ALL",coverage_cutoff = 5,cell_cutoff = 20,gene_cutoff = 5){
  pas <- unique(met$Patient)
  resall <- list()
  if(length(pas) == 1){
    restmp <- step1_cor(exp,mat,annj,met,pa = "ALL",coverage_cutoff,cell_cutoff,gene_cutoff)
    restmp <- list(restmp)
    names(restmp) <- "ALL"
    resall <- restmp
  }else{
    for(n in 1:length(pas)){
      restmp <- step1_cor(exp,mat,annj,met,pa = pas[n],coverage_cutoff,cell_cutoff,gene_cutoff)
      restmp <- list(restmp)
      names(restmp) <- pas[n]
      resall <- c(resall,restmp)
    }
  }
  return(resall)
}
