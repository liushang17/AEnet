library(dplyr)

#' @title Merge Correlation
#' @description The function to merge the correlation across samples from the step1
#'
#' @param cor_list the list containing data.frame from the above steps for each sample
#'
#' @return a data frame containing the correlation link (AEN linke) from multiple samples.
#' @export
#'
#' @examples merge_cor <- merge_cor(cor_list)


#### The output form 
merge_cor <- function(cor_list){
  if(length(cor_list) == 1){
    med2 <- cor_list[[1]]
    med2$junction <- med2$alljunction
    med2$symbol <- med2$geneID
    return(med2)
  }else{
    mit <- NULL
    for(i in 1:length(cor_list)){
      mat <- cor_list[[i]]
      
      if(class(mat) == "data.frame"){
        pos <- which(mat$symbol %in% NA)
        if(length(pos) > 0){mat <- mat[-pos,]}
        mat$junction <- mat$alljunction
        mat$symbol <- mat$geneID
        mit <- rbind(mit,mat)
      }
    }
    
    mit$type <- "+"
    pos <- which(mit$correlation < 0)
    mit$type[pos] <- "-"
    
    mit$all_type <- paste0(mit$alljunction,":",mit$symbol,":",mit$type)
#    saveRDS(mit,file = paste0(oudir,"/","All.cor.rds"))
    
    sui <- data.frame(table(mit$all_type))
    pos <- which(sui$Freq >=2)
    sui1 <- sui[pos,]
    
    pos <- which(mit$all_type %in% sui1$Var1)
    mit1 <- mit[pos,]
#    saveRDS(mit1,file = paste0(oudir,"/","All.cor.filter.rds"))
    
    mit1$type1 <- paste0(mit1$alljunction,":",mit1$symbol)
    mit2 <- mit1[,c("all_type","type1")]
    mit2 <- unique(mit2)
    
    sui2 <- data.frame(table(mit2$type1))
    pos <- which(sui2$Freq >= 2)
    if(length(pos) > 0){
      sui3 <- sui2[pos,]
      pos <- which(mit1$type1 %in% sui3$Var1)
      mit1 <- mit1[-pos,]
    }
    
#    saveRDS(mit1,file = paste0(oudir,"/","All.cor.filter.v2.rds"))
    
    colsToKeep <- "correlation"
    med=mit1 %>%
      group_by(type1) %>%
      summarise_at(vars(colsToKeep), mean)
    med <- data.frame(med)
    med$asp <- gsub(":.*","",med$type1)
    med$symbol <- gsub(".*:","",med$type1)
    med1 <- med[,c("asp","symbol","correlation","type1")]
    colnames(med1)[1] <- "junction"
#    saveRDS(med1,file = paste0(oudir,"/","All.cor.final.rds"))
    sui <- data.frame(table(mit1$type1))
    colnames(sui) <- c("type1","Freq")
    med2 <- merge(med1,sui)
 #   saveRDS(med2,file = paste0(oudir,"/","All.cor.final1.rds"))
    return(med2)
    
  }
}